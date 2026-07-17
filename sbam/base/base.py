import math
import pysam
import numpy as np
from collections import defaultdict
from scipy.stats import binomtest
from Bio import SeqIO

class BaseAccuracy:
    def __init__(self, bam_path, fasta_path, threads=4, k_list=None, min_concordance=0.75, depth_variance_limit=2.5):
        self.bam_path = bam_path
        self.fasta_path = fasta_path
        self.k_list = k_list if k_list else [4, 5, 6]
        self.min_concordance = min_concordance
        self.depth_variance_limit = depth_variance_limit
        
        # Load the longest contig
        self.primary_contig = max(SeqIO.parse(self.fasta_path, "fasta"), key=lambda r: len(r.seq))
        self.seq_len = len(self.primary_contig.seq)
        self.original_contig_id = self.primary_contig.id
        self.bam_contig_id = f"{self.original_contig_id}_cyclic"
        self.ref_seq = str(self.primary_contig.seq).upper()

    def build_fidelity_profile(self):
        """
        Generates genome-wide depth and concordance arrays
        """
        print(f" > Profiling base concordance for {self.original_contig_id} ({self.seq_len} bp)...")
        
        depth_arr = np.zeros(self.seq_len, dtype=int)
        concord_arr = np.ones(self.seq_len, dtype=float) # Default perfect concordance
        
        # Open BAM file
        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            # truncate=True ensures to only read exactly within [0, seq_len)
            for pileupcolumn in bam.pileup(self.bam_contig_id, 0, self.seq_len, truncate=True, min_base_quality=0):
                pos = pileupcolumn.reference_pos
                depth = pileupcolumn.nsegments
                
                if depth == 0 or pos >= self.seq_len:
                    continue
                    
                ref_base = self.ref_seq[pos]
                
                # fast path get all bases at this position instantly
                # This avoids looping through every single read object
                query_bases = pileupcolumn.get_query_sequences(add_indels=False)
                
                # Count how many uppercase query bases match the reference
                matches = sum(1 for qb in query_bases if qb.upper() == ref_base)
                
                depth_arr[pos] = depth
                concord_arr[pos] = matches / depth if depth > 0 else 0.0

        # Masking structural noise: identify bases with depth outside the expected range
        median_depth = np.median(depth_arr[depth_arr > 0])
        if np.isnan(median_depth):
            median_depth = 1  
            
        upper_depth = median_depth * self.depth_variance_limit
        lower_depth = median_depth / self.depth_variance_limit
        
        # Boolean mask: True for bases to be masked (either too high/low depth or low concordance)
        mask_arr = (depth_arr > upper_depth) | (depth_arr < lower_depth) | (concord_arr < self.min_concordance)
        
        masked_count = np.sum(mask_arr)

        self.masked_bases = int(masked_count)
        self.masked_pct = (self.masked_bases / self.seq_len) * 100 if self.seq_len > 0 else 0.0
        print(f"   - Median Depth: {median_depth:.1f}x")
        print(f"   - Masked {masked_count} bases ({(masked_count/self.seq_len)*100:.2f}%) due to structural noise/discordance.")
        
        return concord_arr, mask_arr

    def analyze_motifs(self):
        """
        Scans unmasked regions to find sequence motifs with statistically significant
        basecalling errors (discordance).
        """
        concord_arr, mask_arr = self.build_fidelity_profile()
        print(f" > [Fidelity] Testing systematic errors in {self.k_list}-mer motifs...")
        
        # mask cumsum allows O(1) checks for whether a k-mer window overlaps any masked bases
        mask_cumsum = np.zeros(self.seq_len + 1, dtype=int)
        mask_cumsum[1:] = np.cumsum(mask_arr)
        
        results = []
        
        for k in self.k_list:
            motif_totals = defaultdict(int)
            motif_discordant = defaultdict(int)
            
            for i in range(self.seq_len - k + 1):
                # check if this k-mer overlaps any masked bases
                if (mask_cumsum[i+k] - mask_cumsum[i]) > 0:
                    continue
                    
                kmer = self.ref_seq[i:i+k]
                if 'N' in kmer:
                    continue
                    
                center_idx = i + (k // 2)
                
                motif_totals[kmer] += 1
                if concord_arr[center_idx] < self.min_concordance:
                    motif_discordant[kmer] += 1
                    
            # statistical testing for each k-mer
            total_bases_checked = sum(motif_totals.values())
            total_discordant_bases = sum(motif_discordant.values())
            
            if total_bases_checked == 0:
                continue
                
            global_error_rate = total_discordant_bases / total_bases_checked
            
            for kmer, total in motif_totals.items():
                if total < 500:
                    continue
                    
                discordant_count = motif_discordant[kmer]
                
                # Binomial test
                test = binomtest(k=discordant_count, n=total, p=global_error_rate, alternative='greater')
                
                if test.pvalue < 0.01:
                    error_rate = discordant_count / total
                    fold_increase = error_rate / global_error_rate
                    
                    if fold_increase > 2.0:
                        results.append({
                            "kmer_size": k,
                            "motif": kmer,
                            "occurrences": total,
                            "discordant_count": discordant_count,
                            "error_rate": error_rate,
                            "fold_increase": fold_increase,
                            "p_value": test.pvalue
                        })
                        
        results.sort(key=lambda x: x['fold_increase'], reverse=True)
        return results