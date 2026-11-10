import math
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pysam
from Bio import SeqIO
from scipy.stats import binomtest

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def wilson_lower_bound(k: int, n: int, z: float = 1.96) -> float:
    """Computes the lower bound of the 95% Wilson score confidence interval for a binomial proportion."""
    if n <= 0 or k <= 0:
        return 0.0
    p_hat = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = p_hat + z2 / (2.0 * n)
    rad = z * math.sqrt((p_hat * (1.0 - p_hat) + z2 / (4.0 * n)) / n)
    return max(0.0, (center - rad) / denom)


def compute_empirical_qv(error_count: int, evaluated_bases: int) -> tuple[float, bool]:
    """
    Computes an empirical read-concordance QV proxy: -10 * log10(error_count / evaluated_bases).
    If error_count == 0, returns the lower-bound QV using 0.5 / evaluated_bases and is_capped=True.
    """
    if evaluated_bases <= 0:
        return 0.0, False
    if error_count <= 0:
        return round(-10.0 * math.log10(0.5 / evaluated_bases), 1), True
    return round(-10.0 * math.log10(error_count / evaluated_bases), 1), False


def _share_core_or_substring(m1_pair: str, m2_pair: str) -> bool:
    """Checks if two canonical motif pairs ('FWD/REV') are nested substrings or share a >=5bp core."""
    s1_f, s1_r = m1_pair.split("/")
    s2_f, s2_r = m2_pair.split("/")
    for a in (s1_f, s1_r):
        for b in (s2_f, s2_r):
            if a in b or b in a:
                return True
            if len(a) >= 5 and len(b) >= 5:
                if (
                    a[1:] == b[:-1]
                    or b[1:] == a[:-1]
                    or a[1:] == b[1:]
                    or a[:-1] == b[:-1]
                ):
                    return True
    return False


def _multi_thread(args):
    bam_path, contig, start, end, ref_seq_chunk, max_depth = args
    chunk_len = end - start

    depths = np.zeros(chunk_len, dtype=np.int32)
    matches = np.zeros(chunk_len, dtype=np.int32)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for pileupcolumn in bam.pileup(
            contig,
            start,
            end,
            truncate=True,
            min_base_quality=0,
            max_depth=max_depth,
        ):
            pos = pileupcolumn.reference_pos
            local_idx = pos - start

            if local_idx < 0 or local_idx >= chunk_len:
                continue

            depth = pileupcolumn.nsegments
            if depth == 0:
                continue

            ref_base = ref_seq_chunk[local_idx]

            # Strip deletion lookaheads ('-1X') on the preceding base so deletions are
            # counted once at the deleted position ('*'), while insertions ('+1X') are
            # counted once at the preceding base.
            query_bases = pileupcolumn.get_query_sequences(add_indels=True)
            actual_depth = len(query_bases) if query_bases else depth
            depths[local_idx] = actual_depth

            match_count = sum(
                1 for qb in query_bases if qb.split("-")[0].upper() == ref_base
            )
            matches[local_idx] = match_count

    return start, end, depths, matches


class BaseAccuracy:
    def __init__(
        self,
        bam_path,
        fasta_path,
        threads=4,
        k_list=None,
        min_concordance=0.85,
        majority_threshold=0.50,
        depth_variance_limit=2.5,
        merge_gap_bp=500,
        min_structural_block_bp=20,
        max_depth=100000,
        canonical=True,
    ):
        self.bam_path = bam_path
        self.fasta_path = fasta_path
        self.threads = threads
        self.k_list = k_list if k_list else [4, 5, 6]
        self.min_concordance = min_concordance
        self.majority_threshold = majority_threshold
        self.depth_variance_limit = depth_variance_limit
        self.merge_gap_bp = merge_gap_bp
        self.min_structural_block_bp = min_structural_block_bp
        self.max_depth = max_depth
        self.canonical = canonical

        self.records = list(SeqIO.parse(self.fasta_path, "fasta"))
        self.primary_contig = max(self.records, key=lambda r: len(r.seq))
        self.seq_len = len(self.primary_contig.seq)

        self.original_contig_id = self.primary_contig.id
        self.bam_contig_id = f"{self.original_contig_id}_cyclic"
        self.ref_seq = str(self.primary_contig.seq).upper()

        self.masked_bases = 0
        self.masked_pct = 0.0
        self.masked_regions = []
        self.low_concordance_blocks = []

        self.evaluated_bases = 0
        self.discordant_bases = 0
        self.consensus_error_bases = 0
        self.subconsensus_bias_bases = 0
        self.baseline_discordance_rate = 0.0
        self.consensus_error_rate = 0.0
        self.consensus_errors_per_mb = 0.0
        self.subconsensus_bias_per_mb = 0.0
        self.empirical_qv = 0.0
        self.empirical_qv_capped = False

        self.contig_base_metrics = {}
        self.raw_motif_hits = []
        self.explained_discordant_bases = 0
        self.explained_discordant_pct = 0.0

    def _get_buffer_length_for(self, bam_contig_id: str, orig_len: int) -> int:
        try:
            with pysam.AlignmentFile(self.bam_path, "rb") as bam:
                if bam_contig_id in bam.references:
                    bam_ref_len = bam.get_reference_length(bam_contig_id)
                    if isinstance(bam_ref_len, int) and bam_ref_len > orig_len:
                        return min(bam_ref_len - orig_len, orig_len)
        except Exception:
            pass
        return 0

    def _get_buffer_length(self) -> int:
        return self._get_buffer_length_for(self.bam_contig_id, self.seq_len)

    def _profile_single_contig(self, contig_id: str, seq_str: str, pool=None):
        """Profiles base-level depth, structural masking, and consensus errors for one contig."""
        seq_len = len(seq_str)
        bam_contig_id = f"{contig_id}_cyclic"
        buffer_len = (
            self._get_buffer_length()
            if contig_id == self.original_contig_id
            else self._get_buffer_length_for(bam_contig_id, seq_len)
        )
        total_profile_len = seq_len + buffer_len
        cyclic_ref_seq = seq_str + seq_str[:buffer_len]

        chunk_size = max(1, math.ceil(total_profile_len / self.threads))
        tasks = []
        for i in range(self.threads):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, total_profile_len)
            if start >= end:
                break
            ref_seq_chunk = cyclic_ref_seq[start:end]
            tasks.append(
                (
                    self.bam_path,
                    bam_contig_id,
                    start,
                    end,
                    ref_seq_chunk,
                    self.max_depth,
                )
            )

        raw_depth_arr = np.zeros(total_profile_len, dtype=int)
        raw_match_arr = np.zeros(total_profile_len, dtype=int)

        if pool is not None:
            results = pool.map(_multi_thread, tasks)
        else:
            with mp.Pool(processes=self.threads) as local_pool:
                results = local_pool.map(_multi_thread, tasks)

        for start, end, d_chunk, m_chunk in results:
            raw_depth_arr[start:end] = d_chunk
            raw_match_arr[start:end] = m_chunk

        depth_arr = raw_depth_arr[:seq_len].copy()
        match_arr = raw_match_arr[:seq_len].copy()
        if buffer_len > 0:
            depth_arr[:buffer_len] += raw_depth_arr[seq_len : seq_len + buffer_len]
            match_arr[:buffer_len] += raw_match_arr[seq_len : seq_len + buffer_len]

        concord_arr = np.ones(seq_len, dtype=float)
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

        depth_outlier_mask = (depth_arr > upper_depth) | (depth_arr < lower_depth)

        raw_discord = concord_arr < self.min_concordance
        discord_edges = np.diff(np.concatenate(([0], raw_discord.view(np.int8), [0])))
        d_starts = np.where(discord_edges == 1)[0]
        d_ends = np.where(discord_edges == -1)[0]

        contiguous_discord_mask = np.zeros(seq_len, dtype=bool)
        for s, e in zip(d_starts, d_ends):
            if (e - s) >= self.min_structural_block_bp:
                contiguous_discord_mask[s:e] = True

        structural_mask = depth_outlier_mask | contiguous_discord_mask
        masked_bases = int(np.sum(structural_mask))
        masked_pct = (masked_bases / seq_len * 100.0) if seq_len > 0 else 0.0

        unmasked_flags = ~structural_mask
        evaluated_bases = int(np.sum(unmasked_flags))
        discordant_flags = unmasked_flags & (concord_arr < self.min_concordance)
        majority_err_flags = unmasked_flags & (concord_arr < self.majority_threshold)
        minority_bias_flags = discordant_flags & (~majority_err_flags)

        discordant_bases = int(np.sum(discordant_flags))
        consensus_error_bases = int(np.sum(majority_err_flags))
        subconsensus_bias_bases = int(np.sum(minority_bias_flags))

        s_edges = np.diff(np.concatenate(([0], structural_mask.view(np.int8), [0])))
        starts = np.where(s_edges == 1)[0]
        ends = np.where(s_edges == -1)[0]

        structural_blocks = [
            [int(s), int(e)]
            for s, e in zip(starts, ends)
            if (e - s) >= self.min_structural_block_bp
        ]

        merged_bounds = []
        for s, e in structural_blocks:
            if merged_bounds and (s - merged_bounds[-1][1]) < self.merge_gap_bp:
                merged_bounds[-1][1] = e
            else:
                merged_bounds.append([s, e])

        intervals = []
        low_concord_blocks = []
        for s, e in merged_bounds:
            if e - s >= 50:
                mean_depth = float(np.mean(depth_arr[s:e]))
                mean_concord = float(np.mean(concord_arr[s:e])) * 100.0
                rel_depth = mean_depth / median_depth if median_depth > 0 else 0.0

                if mean_depth > upper_depth:
                    reason = f"High Depth ({rel_depth:.1f}x of median; Collapsed Repeat / Multi-copy Operon)"
                elif mean_depth < 10.0:
                    reason = "Low Depth (<10x; Coverage Gap / Dropout)"
                elif mean_depth < lower_depth:
                    reason = (
                        f"Relative Depth Dip ({rel_depth:.2f}x of median; "
                        f"Structural Heterogeneity / Prophage / Mappability)"
                    )
                else:
                    reason = "Contiguous Low Concordance (Unpolished Block / Structural Variant)"

                block_info = {
                    "contig": contig_id,
                    "start": int(s),
                    "end": int(e),
                    "length": int(e - s),
                    "mean_depth": round(mean_depth, 1),
                    "mean_concordance": round(mean_concord, 1),
                    "reason": reason,
                }
                intervals.append(block_info)
                if mean_concord < 90.0:
                    low_concord_blocks.append(block_info)

        intervals.sort(key=lambda x: x["length"], reverse=True)
        qv_val, qv_capped = compute_empirical_qv(
            consensus_error_bases, evaluated_bases
        )
        cons_per_mb = (
            consensus_error_bases / (evaluated_bases / 1_000_000)
            if evaluated_bases > 0
            else 0.0
        )

        metrics = {
            "length": seq_len,
            "median_depth": round(median_depth, 1),
            "masked_bases": masked_bases,
            "masked_pct": round(masked_pct, 2),
            "evaluated_bases": evaluated_bases,
            "discordant_bases": discordant_bases,
            "consensus_error_bases": consensus_error_bases,
            "subconsensus_bias_bases": subconsensus_bias_bases,
            "consensus_errors_per_mb": round(cons_per_mb, 2),
            "empirical_qv": qv_val,
            "empirical_qv_capped": qv_capped,
            "low_concordance_blocks": low_concord_blocks,
            "masked_regions": intervals[:50],
        }
        return concord_arr, structural_mask, metrics

    def build_fidelity_profile(self):
        print(
            f" > [Fidelity] Profiling base concordance across {len(self.records)} contig(s) "
            f"(primary: {self.original_contig_id}, {self.seq_len} bp) using {self.threads} threads..."
        )

        primary_concord_arr = None
        primary_structural_mask = None

        with mp.Pool(processes=self.threads) as pool:
            for rec in self.records:
                cid = rec.id
                seq_str = str(rec.seq).upper()
                try:
                    c_arr, s_mask, m_dict = self._profile_single_contig(
                        cid, seq_str, pool=pool
                    )
                    self.contig_base_metrics[cid] = m_dict
                    if cid == self.original_contig_id:
                        primary_concord_arr = c_arr
                        primary_structural_mask = s_mask
                except ValueError:
                    # Secondary contig absent in mocked test BAM
                    if cid == self.original_contig_id:
                        raise

        p_metrics = self.contig_base_metrics[self.original_contig_id]
        self.masked_bases = p_metrics["masked_bases"]
        self.masked_pct = (
            (self.masked_bases / self.seq_len) * 100.0 if self.seq_len > 0 else 0.0
        )
        self.evaluated_bases = p_metrics["evaluated_bases"]
        self.discordant_bases = p_metrics["discordant_bases"]
        self.consensus_error_bases = p_metrics["consensus_error_bases"]
        self.subconsensus_bias_bases = p_metrics["subconsensus_bias_bases"]
        self.baseline_discordance_rate = (
            self.discordant_bases / self.evaluated_bases
            if self.evaluated_bases > 0
            else 0.0
        )
        self.consensus_error_rate = (
            self.consensus_error_bases / self.evaluated_bases
            if self.evaluated_bases > 0
            else 0.0
        )
        self.consensus_errors_per_mb = p_metrics["consensus_errors_per_mb"]
        self.subconsensus_bias_per_mb = (
            round(self.subconsensus_bias_bases / (self.evaluated_bases / 1_000_000), 1)
            if self.evaluated_bases > 0
            else 0.0
        )
        self.empirical_qv = p_metrics["empirical_qv"]
        self.empirical_qv_capped = p_metrics["empirical_qv_capped"]
        self.masked_regions = p_metrics["masked_regions"]
        self.low_concordance_blocks = p_metrics["low_concordance_blocks"]

        qv_str = (
            f">{self.empirical_qv:.1f}"
            if self.empirical_qv_capped
            else f"{self.empirical_qv:.1f}"
        )
        print(f"   - Primary Median Depth: {p_metrics['median_depth']:.1f}x")
        print(
            f"   - Structurally Masked: {self.masked_bases} bp ({self.masked_pct:.2f}%)"
        )
        print(
            f"   - Evaluated Positions: {self.evaluated_bases} bp | "
            f"Total Discordant (<{self.min_concordance:.0%}): {self.discordant_bases} bp "
            f"({self.baseline_discordance_rate:.4%})"
        )
        print(
            f"     * Consensus Errors (<{self.majority_threshold:.0%} support): "
            f"{self.consensus_error_bases} bp ({self.consensus_errors_per_mb:.2f}/Mb, "
            f"Concordance QV ≈ {qv_str}) | "
            f"Sub-consensus Read Bias ({self.majority_threshold:.0%}–{self.min_concordance:.0%}): "
            f"{self.subconsensus_bias_bases} bp ({self.subconsensus_bias_per_mb:.1f}/Mb)"
        )

        return primary_concord_arr, primary_structural_mask

    def analyze_motifs(self, min_occurrences=100, min_discordant_count=10):
        """
        Scans unmasked regions on the primary chromosome for k-mer contexts enriched
        for read-to-assembly discordance (< min_concordance). Uses a 3-state empirical
        orientation classification ('Orientation-Asymmetric', 'Orientation-Balanced',
        'Undetermined') and assigns discordant sites mutually exclusively across families
        so the table rows sum to the exact headline union.
        """
        concord_arr, structural_mask = self.build_fidelity_profile()
        print(f" > [Fidelity] Testing systematic errors in {self.k_list}-mer motifs...")

        if self.evaluated_bases == 0 or self.discordant_bases == 0:
            self.raw_motif_hits = []
            self.explained_discordant_bases = 0
            self.explained_discordant_pct = 0.0
            return []

        global_error_rate = self.baseline_discordance_rate

        mask_cumsum = np.zeros(self.seq_len + 1, dtype=int)
        mask_cumsum[1:] = np.cumsum(structural_mask)

        candidate_tests = []

        for k in self.k_list:
            motif_totals = defaultdict(int)
            motif_discordant = defaultdict(int)
            motif_fwd_totals = defaultdict(int)
            motif_rev_totals = defaultdict(int)
            motif_fwd_disc = defaultdict(int)
            motif_rev_disc = defaultdict(int)
            motif_majority_err = defaultdict(int)
            motif_is_palindrome = {}
            motif_positions = defaultdict(set)
            motif_majority_positions = defaultdict(set)

            for i in range(self.seq_len - k + 1):
                if (mask_cumsum[i + k] - mask_cumsum[i]) > 0:
                    continue

                raw_kmer = self.ref_seq[i : i + k]
                if "N" in raw_kmer:
                    continue

                if not self.canonical:
                    pos = i + (k // 2)
                    motif_is_palindrome[raw_kmer] = False
                    motif_totals[raw_kmer] += 1
                    motif_fwd_totals[raw_kmer] += 1
                    if concord_arr[pos] < self.min_concordance:
                        motif_discordant[raw_kmer] += 1
                        motif_fwd_disc[raw_kmer] += 1
                        motif_positions[raw_kmer].add(pos)
                        if concord_arr[pos] < self.majority_threshold:
                            motif_majority_err[raw_kmer] += 1
                            motif_majority_positions[raw_kmer].add(pos)
                    continue

                rc_kmer = revcomp(raw_kmer)
                is_pal = raw_kmer == rc_kmer
                is_fwd = raw_kmer <= rc_kmer
                canon_key = (
                    f"{raw_kmer}/{rc_kmer}" if is_fwd else f"{rc_kmer}/{raw_kmer}"
                )
                motif_is_palindrome[canon_key] = is_pal

                offset = (k // 2) if (is_pal or is_fwd) else (k - 1 - (k // 2))
                pos = i + offset

                motif_totals[canon_key] += 1
                if is_pal or is_fwd:
                    motif_fwd_totals[canon_key] += 1
                else:
                    motif_rev_totals[canon_key] += 1

                if concord_arr[pos] < self.min_concordance:
                    motif_discordant[canon_key] += 1
                    if is_pal or is_fwd:
                        motif_fwd_disc[canon_key] += 1
                    else:
                        motif_rev_disc[canon_key] += 1
                    motif_positions[canon_key].add(pos)
                    if concord_arr[pos] < self.majority_threshold:
                        motif_majority_err[canon_key] += 1
                        motif_majority_positions[canon_key].add(pos)

            for kmer_key, total in motif_totals.items():
                if total < min_occurrences:
                    continue
                discordant_count = motif_discordant[kmer_key]
                if discordant_count < min_discordant_count:
                    continue

                test = binomtest(
                    k=discordant_count,
                    n=total,
                    p=global_error_rate,
                    alternative="greater",
                )
                error_rate = discordant_count / total
                fold_increase = error_rate / global_error_rate
                lb_rate = wilson_lower_bound(discordant_count, total)
                ci_lower_fold = lb_rate / global_error_rate

                is_pal = motif_is_palindrome.get(kmer_key, False)
                f_disc = motif_fwd_disc[kmer_key]
                r_disc = motif_rev_disc[kmer_key]
                f_tot = motif_fwd_totals[kmer_key]
                r_tot = motif_rev_totals[kmer_key]

                # 3-State Empirical Orientation Classification (Asymmetric / Balanced / Undetermined)
                orientation_state = "Undetermined"
                strand_biased = False
                strand_pval = 1.0
                if is_pal:
                    orientation_state = "Palindromic"
                elif (f_tot + r_tot) > 0 and discordant_count > 0:
                    p_fwd = f_tot / (f_tot + r_tot)
                    strand_pval = binomtest(
                        k=f_disc, n=discordant_count, p=p_fwd, alternative="two-sided"
                    ).pvalue
                    max_share = max(f_disc, r_disc) / discordant_count
                    if strand_pval < 0.01 and max_share >= 0.80:
                        orientation_state = "Orientation-Asymmetric"
                        strand_biased = True
                    elif (
                        discordant_count >= 20
                        and strand_pval >= 0.05
                        and max_share <= 0.65
                    ):
                        orientation_state = "Orientation-Balanced"
                    else:
                        orientation_state = "Undetermined (Low Count / Intermediate)"

                candidate_tests.append(
                    {
                        "kmer_size": k,
                        "motif": kmer_key,
                        "is_palindrome": is_pal,
                        "occurrences": total,
                        "discordant_count": discordant_count,
                        "consensus_error_count": motif_majority_err[kmer_key],
                        "subconsensus_bias_count": discordant_count
                        - motif_majority_err[kmer_key],
                        "fwd_discordant": f_disc,
                        "rev_discordant": r_disc,
                        "fwd_occurrences": f_tot,
                        "rev_occurrences": r_tot,
                        "strand_biased": strand_biased,
                        "orientation_state": orientation_state,
                        "strand_pvalue": strand_pval,
                        "error_rate": error_rate,
                        "fold_increase": fold_increase,
                        "ci_lower_fold": ci_lower_fold,
                        "p_value": test.pvalue,
                        "_positions": motif_positions[kmer_key],
                        "_maj_positions": motif_majority_positions[kmer_key],
                    }
                )

        if not candidate_tests:
            self.raw_motif_hits = []
            self.explained_discordant_bases = 0
            self.explained_discordant_pct = 0.0
            return []

        num_hypotheses = len(candidate_tests)
        corrected_alpha = 0.01 / num_hypotheses

        significant_hits = []
        for item in candidate_tests:
            if item["p_value"] < corrected_alpha and item["ci_lower_fold"] >= 2.0:
                item["adj_p_value"] = min(1.0, item["p_value"] * num_hypotheses)
                significant_hits.append(item)

        if not significant_hits:
            self.raw_motif_hits = []
            self.explained_discordant_bases = 0
            self.explained_discordant_pct = 0.0
            return []

        families = self._cluster_into_families(significant_hits)

        for item in significant_hits:
            item.pop("_positions", None)
            item.pop("_maj_positions", None)

        self.raw_motif_hits = significant_hits

        print(
            f"   - Significant k-mers: {len(significant_hits)} -> Clustered into "
            f"{len(families)} Motif Family(ies) explaining {self.explained_discordant_bases}/"
            f"{self.discordant_bases} ({self.explained_discordant_pct:.1f}%) of discordant sites."
        )
        return families

    def _cluster_into_families(self, significant_hits):
        sorted_hits = sorted(
            significant_hits,
            key=lambda x: (x["p_value"], -x["ci_lower_fold"], -x["discordant_count"]),
        )

        families = []
        for hit in sorted_hits:
            hit_pos = hit["_positions"]
            hit_maj_pos = hit["_maj_positions"]
            placed = False

            for fam in families:
                lead_pos = fam["_lead_positions"]
                min_size = min(len(hit_pos), len(lead_pos))
                coord_overlap = (
                    (len(hit_pos & lead_pos) / min_size) if min_size > 0 else 0.0
                )

                if coord_overlap >= 0.25 or _share_core_or_substring(
                    hit["motif"], fam["motif"]
                ):
                    fam["member_motifs"].append(hit["motif"])
                    fam["family_size"] += 1
                    fam["_union_positions"].update(hit_pos)
                    fam["_union_maj_positions"].update(hit_maj_pos)
                    if hit["error_rate"] > fam["max_member_error_rate"]:
                        fam["max_member_error_rate"] = hit["error_rate"]
                        fam["max_member_fold"] = hit["fold_increase"]
                        fam["peak_member_motif"] = hit["motif"]
                    placed = True
                    break

            if not placed:
                families.append(
                    {
                        "kmer_size": hit["kmer_size"],
                        "motif": hit["motif"],
                        "is_palindrome": hit["is_palindrome"],
                        "member_motifs": [hit["motif"]],
                        "family_size": 1,
                        "occurrences": hit["occurrences"],
                        "discordant_count": hit["discordant_count"],
                        "consensus_error_count": hit["consensus_error_count"],
                        "subconsensus_bias_count": hit["subconsensus_bias_count"],
                        "fwd_discordant": hit["fwd_discordant"],
                        "rev_discordant": hit["rev_discordant"],
                        "fwd_occurrences": hit["fwd_occurrences"],
                        "rev_occurrences": hit["rev_occurrences"],
                        "strand_biased": hit["strand_biased"],
                        "orientation_state": hit["orientation_state"],
                        "error_rate": hit["error_rate"],
                        "fold_increase": hit["fold_increase"],
                        "ci_lower_fold": round(hit["ci_lower_fold"], 2),
                        "max_member_error_rate": hit["error_rate"],
                        "max_member_fold": hit["fold_increase"],
                        "peak_member_motif": hit["motif"],
                        "p_value": hit["p_value"],
                        "adj_p_value": hit["adj_p_value"],
                        "_lead_positions": set(hit_pos),
                        "_union_positions": set(hit_pos),
                        "_union_maj_positions": set(hit_maj_pos),
                    }
                )

        # Sort families by raw unique discordant sites first
        families.sort(
            key=lambda x: (len(x["_union_positions"]), x["ci_lower_fold"]),
            reverse=True,
        )

        # Partition discordant coordinates mutually exclusively across families
        # so that `assigned_unique_sites` sums to the exact global union (e.g. 1,483)
        claimed_coords = set()
        for fam in families:
            raw_union = fam.pop("_union_positions")
            maj_union = fam.pop("_union_maj_positions")
            fam.pop("_lead_positions", None)

            exclusive_coords = raw_union - claimed_coords
            claimed_coords.update(raw_union)

            raw_count = len(raw_union)
            assigned_count = len(exclusive_coords)
            maj_count = len(maj_union)

            fam["raw_family_sites"] = raw_count
            fam["assigned_unique_sites"] = assigned_count
            fam["family_unique_sites"] = assigned_count
            fam["pct_discordant_explained"] = (
                round(assigned_count / self.discordant_bases * 100.0, 2)
                if self.discordant_bases > 0
                else 0.0
            )
            fam["raw_pct_explained"] = (
                round(raw_count / self.discordant_bases * 100.0, 2)
                if self.discordant_bases > 0
                else 0.0
            )
            fam["consensus_error_count"] = max(
                fam["consensus_error_count"], maj_count
            )

            related = fam["member_motifs"][1:6]
            fam["related_motifs_str"] = ", ".join(related) if related else "-"

            # Empirical, non-causal signal classification
            if fam["consensus_error_count"] > 0:
                fam["artifact_class"] = (
                    f"Consensus Error ({fam['consensus_error_count']} site(s) <50%)"
                )
            else:
                fam["artifact_class"] = fam["orientation_state"]

        self.explained_discordant_bases = len(claimed_coords)
        self.explained_discordant_pct = (
            (self.explained_discordant_bases / self.discordant_bases * 100.0)
            if self.discordant_bases > 0
            else 0.0
        )

        # Re-sort by exclusively assigned sites so the table is ordered cleanly
        families.sort(
            key=lambda x: (
                x["assigned_unique_sites"],
                x["raw_family_sites"],
                x["ci_lower_fold"],
            ),
            reverse=True,
        )
        return families