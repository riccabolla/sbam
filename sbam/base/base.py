import math
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pysam
from Bio import SeqIO
from scipy.stats import binomtest


def _multi_thread(args):

    bam_path, contig, start, end, ref_seq_chunk = args
    chunk_len = end - start
    
    depths = np.zeros(chunk_len, dtype=np.int32)
    matches = np.zeros(chunk_len, dtype=np.int32)
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for pileupcolumn in bam.pileup(contig, start, end, truncate=True, min_base_quality=0):
            pos = pileupcolumn.reference_pos
            local_idx = pos - start
            
            if local_idx < 0 or local_idx >= chunk_len:
                continue
                
            depth = pileupcolumn.nsegments
            depths[local_idx] = depth
            
            if depth == 0:
                continue
                
            ref_base = ref_seq_chunk[local_idx]
            
            query_bases = pileupcolumn.get_query_sequences(add_indels=False)
            match_count = sum(1 for qb in query_bases if qb.upper() == ref_base)
            
            matches[local_idx] = match_count
            
    return start, end, depths, matches


class BaseAccuracy:
    def __init__(self, bam_path, fasta_path, threads=4, k_list=None, min_concordance=0.75, depth_variance_limit=2.5):
        self.bam_path = bam_path
        self.fasta_path = fasta_path
        self.threads = threads
        self.k_list = k_list if k_list else [4, 5, 6]
        self.min_concordance = min_concordance
        self.depth_variance_limit = depth_variance_limit
        
        self.primary_contig = max(SeqIO.parse(self.fasta_path, "fasta"), key=lambda r: len(r.seq))
        self.seq_len = len(self.primary_contig.seq)
        
        self.original_contig_id = self.primary_contig.id
        self.bam_contig_id = f"{self.original_contig_id}_cyclic"
        self.ref_seq = str(self.primary_contig.seq).upper()
        
        self.masked_bases = 0
        self.masked_pct = 0.0

    def build_fidelity_profile(self):
        """
        Chunks the genome and delegates pileup to multiprocessing workers.
        """
        print(f" > [Fidelity] Profiling base concordance for {self.original_contig_id} ({self.seq_len} bp) using {self.threads} threads...")
        
        chunk_size = math.ceil(self.seq_len / self.threads)
        tasks = []
        for i in range(self.threads):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, self.seq_len)
            if start >= end: 
                break
            ref_seq_chunk = self.ref_seq[start:end]
            tasks.append((self.bam_path, self.bam_contig_id, start, end, ref_seq_chunk))
            
        depth_arr = np.zeros(self.seq_len, dtype=int)
        concord_arr = np.ones(self.seq_len, dtype=float) 
        
        # Execute parallel workers
        with mp.Pool(processes=self.threads) as pool:
            results = pool.map(_multi_thread, tasks)
            
        # Reassemble arrays
        for start, end, d_chunk, m_chunk in results:
            depth_arr[start:end] = d_chunk
            # Calculate concordance safely
            valid_mask = d_chunk > 0
            concord_arr[start:end][valid_mask] = m_chunk[valid_mask] / d_chunk[valid_mask]

        median_depth = np.median(depth_arr[depth_arr > 0])
        if np.isnan(median_depth):
            median_depth = 1  
            
        upper_depth = median_depth * self.depth_variance_limit
        lower_depth = median_depth / self.depth_variance_limit
        
        mask_arr = (depth_arr > upper_depth) | (depth_arr < lower_depth) | (concord_arr < self.min_concordance)
        
        self.masked_bases = int(np.sum(mask_arr))
        self.masked_pct = (self.masked_bases / self.seq_len) * 100 if self.seq_len > 0 else 0.0
        
        print(f"   - Median Depth: {median_depth:.1f}x")
        print(f"   - Masked {self.masked_bases} bases ({self.masked_pct:.2f}%) due to structural noise/discordance.")

        intervals = []
        # Fast numpy edge detection to find start/end of true blocks
        edges = np.diff(np.concatenate(([0], mask_arr.view(np.int8), [0])))
        starts = np.where(edges == 1)[0]
        ends = np.where(edges == -1)[0]
        
        for s, e in zip(starts, ends):
            if e - s >= 50:  # Filter out noisy drops (keep blocks >50bp)
                mean_depth = np.mean(depth_arr[s:e])
                if mean_depth > upper_depth: 
                    reason = "High Depth (Collapsed Repeat / Multi-copy Operon)"
                elif mean_depth < lower_depth: 
                    reason = "Low Depth (Gap / Poor Coverage)"
                else: 
                    reason = "Low Concordance (Severe Basecalling Discordance)"
                
                intervals.append({"start": int(s), "end": int(e), "length": int(e-s), "reason": reason})
                
        intervals.sort(key=lambda x: x['length'], reverse=True)
        self.masked_regions = intervals[:50]  # Store top 50 largest regions
        
        return concord_arr, mask_arr

    def analyze_motifs(self):
        """
        Scans unmasked regions to find sequence motifs with statistically significant basecalling errors.
        """
        concord_arr, mask_arr = self.build_fidelity_profile()
        print(f" > [Fidelity] Testing systematic errors in {self.k_list}-mer motifs...")
        
        mask_cumsum = np.zeros(self.seq_len + 1, dtype=int)
        mask_cumsum[1:] = np.cumsum(mask_arr)
        
        results = []
        
        for k in self.k_list:
            motif_totals = defaultdict(int)
            motif_discordant = defaultdict(int)
            
            for i in range(self.seq_len - k + 1):
                # O(1) mask lookup
                if (mask_cumsum[i+k] - mask_cumsum[i]) > 0:
                    continue
                    
                kmer = self.ref_seq[i:i+k]
                if 'N' in kmer:
                    continue
                    
                center_idx = i + (k // 2)
                
                motif_totals[kmer] += 1
                if concord_arr[center_idx] < self.min_concordance:
                    motif_discordant[kmer] += 1
                    
            total_bases_checked = sum(motif_totals.values())
            total_discordant_bases = sum(motif_discordant.values())
            
            if total_bases_checked == 0:
                continue
                
            global_error_rate = total_discordant_bases / total_bases_checked
            
            for kmer, total in motif_totals.items():
                if total < 500:
                    continue
                    
                discordant_count = motif_discordant[kmer]
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