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
            if depth == 0:
                continue

            ref_base = ref_seq_chunk[local_idx]

            # add_indels=True ensures insertions, deletions, and substitutions
            # all register as discordant against ref_base.
            query_bases = pileupcolumn.get_query_sequences(add_indels=True)
            actual_depth = len(query_bases) if query_bases else depth
            depths[local_idx] = actual_depth

            match_count = sum(1 for qb in query_bases if qb.upper() == ref_base)
            matches[local_idx] = match_count

    return start, end, depths, matches


class BaseAccuracy:
    def __init__(
        self,
        bam_path,
        fasta_path,
        threads=4,
        k_list=None,
        min_concordance=0.75,
        depth_variance_limit=2.5,
        merge_gap_bp=500,
    ):
        self.bam_path = bam_path
        self.fasta_path = fasta_path
        self.threads = threads
        self.k_list = k_list if k_list else [4, 5, 6]
        self.min_concordance = min_concordance
        self.depth_variance_limit = depth_variance_limit
        self.merge_gap_bp = merge_gap_bp

        self.primary_contig = max(
            SeqIO.parse(self.fasta_path, "fasta"), key=lambda r: len(r.seq)
        )
        self.seq_len = len(self.primary_contig.seq)

        self.original_contig_id = self.primary_contig.id
        self.bam_contig_id = f"{self.original_contig_id}_cyclic"
        self.ref_seq = str(self.primary_contig.seq).upper()

        self.masked_bases = 0
        self.masked_pct = 0.0
        self.masked_regions = []

    def _get_buffer_length(self):
        try:
            with pysam.AlignmentFile(self.bam_path, "rb") as bam:
                if self.bam_contig_id in bam.references:
                    bam_ref_len = bam.get_reference_length(self.bam_contig_id)
                    if isinstance(bam_ref_len, int) and bam_ref_len > self.seq_len:
                        return min(bam_ref_len - self.seq_len, self.seq_len)
        except Exception:
            pass
        return 0

    def build_fidelity_profile(self):
        """
        Chunks the genome (including any cyclic buffer), delegates pileup to
        workers, folds buffer coverage back into [0, buffer_len), and merges
        nearby masked intervals (<500 bp apart).
        """
        print(
            f" > [Fidelity] Profiling base concordance for {self.original_contig_id} "
            f"({self.seq_len} bp) using {self.threads} threads..."
        )

        buffer_len = self._get_buffer_length()
        total_profile_len = self.seq_len + buffer_len
        cyclic_ref_seq = self.ref_seq + self.ref_seq[:buffer_len]

        chunk_size = math.ceil(total_profile_len / self.threads)
        tasks = []
        for i in range(self.threads):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, total_profile_len)
            if start >= end:
                break
            ref_seq_chunk = cyclic_ref_seq[start:end]
            tasks.append((self.bam_path, self.bam_contig_id, start, end, ref_seq_chunk))

        raw_depth_arr = np.zeros(total_profile_len, dtype=int)
        raw_match_arr = np.zeros(total_profile_len, dtype=int)

        with mp.Pool(processes=self.threads) as pool:
            results = pool.map(_multi_thread, tasks)

        for start, end, d_chunk, m_chunk in results:
            raw_depth_arr[start:end] = d_chunk
            raw_match_arr[start:end] = m_chunk

        # Fold cyclic buffer depths and matches back into the start of the contig
        depth_arr = raw_depth_arr[: self.seq_len].copy()
        match_arr = raw_match_arr[: self.seq_len].copy()
        if buffer_len > 0:
            depth_arr[:buffer_len] += raw_depth_arr[
                self.seq_len : self.seq_len + buffer_len
            ]
            match_arr[:buffer_len] += raw_match_arr[
                self.seq_len : self.seq_len + buffer_len
            ]

        concord_arr = np.ones(self.seq_len, dtype=float)
        valid_mask = depth_arr > 0
        concord_arr[valid_mask] = match_arr[valid_mask] / depth_arr[valid_mask]
        concord_arr[~valid_mask] = 0.0

        non_zero_depths = depth_arr[valid_mask]
        median_depth = (
            float(np.median(non_zero_depths)) if non_zero_depths.size > 0 else 1.0
        )
        if np.isnan(median_depth) or median_depth <= 0:
            median_depth = 1.0

        upper_depth = median_depth * self.depth_variance_limit
        lower_depth = median_depth / self.depth_variance_limit

        mask_arr = (
            (depth_arr > upper_depth)
            | (depth_arr < lower_depth)
            | (concord_arr < self.min_concordance)
        )

        self.masked_bases = int(np.sum(mask_arr))
        self.masked_pct = (
            (self.masked_bases / self.seq_len) * 100 if self.seq_len > 0 else 0.0
        )

        print(f"   - Median Depth: {median_depth:.1f}x")
        print(
            f"   - Masked {self.masked_bases} bases ({self.masked_pct:.2f}%) "
            f"due to structural noise/discordance."
        )

        edges = np.diff(np.concatenate(([0], mask_arr.view(np.int8), [0])))
        starts = np.where(edges == 1)[0]
        ends = np.where(edges == -1)[0]

        # Only keep masked intervals that are at least 20 bp long
        structural_blocks = [
            [int(s), int(e)] for s, e in zip(starts, ends) if (e - s) >= 20
        ]

        # Merge adjacent masked intervals separated by < self.merge_gap_bp
        merged_bounds = []
        for s, e in structural_blocks:
            if merged_bounds and (s - merged_bounds[-1][1]) < self.merge_gap_bp:
                merged_bounds[-1][1] = e
            else:
                merged_bounds.append([s,e])

        intervals = []
        for s, e in merged_bounds:
            if e - s >= 50:
                mean_depth = float(np.mean(depth_arr[s:e]))
                mean_concord = float(np.mean(concord_arr[s:e])) * 100.0
                if mean_depth > upper_depth:
                    reason = "High Depth (Collapsed Repeat / Multi-copy Operon)"
                elif mean_depth < lower_depth:
                    reason = "Low Depth (Gap / Poor Coverage)"
                else:
                    reason = "Low Concordance (Severe Basecalling Discordance)"

                intervals.append(
                    {
                        "start": int(s),
                        "end": int(e),
                        "length": int(e - s),
                        "mean_depth": round(mean_depth, 1),
                        "mean_concordance": round(mean_concord, 1),
                        "reason": reason,
                    }
                )

        intervals.sort(key=lambda x: x["length"], reverse=True)
        self.masked_regions = intervals[:50]

        return concord_arr, mask_arr

    def analyze_motifs(self):
        """
        Scans unmasked regions to find sequence motifs with statistically
        significant basecalling errors (including indels and substitutions).
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
                if (mask_cumsum[i + k] - mask_cumsum[i]) > 0:
                    continue

                kmer = self.ref_seq[i : i + k]
                if "N" in kmer:
                    continue

                center_idx = i + (k // 2)

                motif_totals[kmer] += 1
                if concord_arr[center_idx] < self.min_concordance:
                    motif_discordant[kmer] += 1

            total_bases_checked = sum(motif_totals.values())
            total_discordant_bases = sum(motif_discordant.values())

            if total_bases_checked == 0 or total_discordant_bases == 0:
                continue

            global_error_rate = total_discordant_bases / total_bases_checked

            eligible_motifs = {
                kmer: total for kmer, total in motif_totals.items() if total >= 500
            }
            if not eligible_motifs:
                continue

            corrected_alpha = 0.01 / len(eligible_motifs)

            for kmer, total in eligible_motifs.items():
                discordant_count = motif_discordant[kmer]
                if discordant_count == 0:
                    continue

                test = binomtest(
                    k=discordant_count,
                    n=total,
                    p=global_error_rate,
                    alternative="greater",
                )

                if test.pvalue < corrected_alpha:
                    error_rate = discordant_count / total
                    fold_increase = error_rate / global_error_rate

                    if fold_increase > 2.0:
                        results.append(
                            {
                                "kmer_size": k,
                                "motif": kmer,
                                "occurrences": total,
                                "discordant_count": discordant_count,
                                "error_rate": error_rate,
                                "fold_increase": fold_increase,
                                "p_value": test.pvalue,
                                "adj_p_value": min(
                                    1.0, test.pvalue * len(eligible_motifs)
                                ),
                            }
                        )

        results.sort(key=lambda x: x["fold_increase"], reverse=True)
        return results