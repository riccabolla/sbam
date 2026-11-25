"""Base-level read concordance, structural masking, strand/indel stratification, and motif analysis."""

from __future__ import annotations

import logging
import math
import multiprocessing as mp
from collections import defaultdict
from typing import Any

import numpy as np
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scipy.stats import binomtest, fisher_exact

from sbam.utils import compute_empirical_qv, revcomp, wilson_lower_bound

logger = logging.getLogger("sbam")

_BASE_LUT = np.full(256, 4, dtype=np.uint8)
for _ch, _code in ((ord("A"), 0), (ord("C"), 1), (ord("G"), 2), (ord("T"), 3)):
    _BASE_LUT[_ch] = _code
    _BASE_LUT[_ch + 32] = _code

_CODE_TO_BASE = ("A", "C", "G", "T")
_DEFAULT_CHUNK_BP = 250_000


def _decode_kmer(kid: int, k: int) -> str:
    chars = [""] * k
    val = int(kid)
    for idx in range(k - 1, -1, -1):
        chars[idx] = _CODE_TO_BASE[val & 3]
        val >>= 2
    return "".join(chars)


def _revcomp_kmer_id(kid: int, k: int) -> int:
    rc = 0
    val = int(kid)
    for _ in range(k):
        rc = (rc << 2) | ((val & 3) ^ 3)
        val >>= 2
    return rc


def _build_kmer_lut(k: int, canonical: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Precompute canonical mapping LUTs for all 4^k k-mers (k <= 7)."""
    n_kmers = 1 << (2 * k)
    canon_lut = np.empty(n_kmers, dtype=np.int32)
    is_fwd_lut = np.ones(n_kmers, dtype=bool)
    is_pal_lut = np.zeros(n_kmers, dtype=bool)

    if not canonical:
        labels = [_decode_kmer(kid, k) for kid in range(n_kmers)]
        canon_lut[:] = np.arange(n_kmers, dtype=np.int32)
        return canon_lut, is_fwd_lut, is_pal_lut, labels

    canon_to_idx: dict[int, int] = {}
    labels: list[str] = []

    for kid in range(n_kmers):
        rc_id = _revcomp_kmer_id(kid, k)
        is_pal = kid == rc_id
        is_fwd = kid <= rc_id
        canon_id = kid if is_fwd else rc_id
        idx = canon_to_idx.get(canon_id)
        if idx is None:
            idx = len(labels)
            canon_to_idx[canon_id] = idx
            f_str = _decode_kmer(canon_id, k)
            r_str = _decode_kmer(_revcomp_kmer_id(canon_id, k), k)
            labels.append(f"{f_str}/{r_str}")
        canon_lut[kid] = idx
        is_fwd_lut[kid] = is_pal or is_fwd
        is_pal_lut[kid] = is_pal

    return canon_lut, is_fwd_lut, is_pal_lut, labels


def _share_core_or_substring(m1_pair: str, m2_pair: str) -> bool:
    """Check whether two motifs (canonical 'FWD/REV' or raw) are nested or share a >=5 bp core."""
    parts1 = m1_pair.split("/")
    parts2 = m2_pair.split("/")
    for a in parts1:
        for b in parts2:
            if a in b or b in a:
                return True
            if len(a) >= 6 and len(b) >= 6:
                if (
                    a[1:] == b[:-1]
                    or b[1:] == a[:-1]
                    or a[1:] == b[1:]
                    or a[:-1] == b[:-1]
                ):
                    return True
    return False


def _profile_chunk(args: tuple[str, str, str, int, int, str, int]):
    """Read-centric CIGAR sweep over [start, end) using NumPy difference arrays."""
    contig_id, bam_path, bam_contig, start, end, ref_seq_chunk, _max_depth = args
    chunk_len = end - start
    diff_len = chunk_len + 1

    depth_diff = np.zeros(diff_len, dtype=np.int32)
    match_diff = np.zeros(diff_len, dtype=np.int32)
    mismatch_diff = np.zeros(diff_len, dtype=np.int32)
    del_diff = np.zeros(diff_len, dtype=np.int32)
    fwd_depth_diff = np.zeros(diff_len, dtype=np.int32)
    fwd_match_diff = np.zeros(diff_len, dtype=np.int32)
    rev_depth_diff = np.zeros(diff_len, dtype=np.int32)
    rev_match_diff = np.zeros(diff_len, dtype=np.int32)
    ins_arr = np.zeros(chunk_len, dtype=np.int32)

    direct_match: np.ndarray | None = None
    direct_mismatch: np.ndarray | None = None
    direct_fwd_match: np.ndarray | None = None
    direct_rev_match: np.ndarray | None = None
    ref_bytes: np.ndarray | None = None

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(bam_contig, start, end):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_qcfail
                or read.is_duplicate
            ):
                continue

            cigartuples = read.cigartuples
            rpos = read.reference_start
            if cigartuples is None or rpos is None:
                continue

            qpos = 0
            is_rev = read.is_reverse
            prev_matched_pos = -1

            for op, length in cigartuples:
                if op == 7:  # = (sequence match)
                    a = max(0, rpos - start)
                    b = min(chunk_len, rpos + length - start)
                    if a < b:
                        depth_diff[a] += 1
                        depth_diff[b] -= 1
                        match_diff[a] += 1
                        match_diff[b] -= 1
                        if is_rev:
                            rev_depth_diff[a] += 1
                            rev_depth_diff[b] -= 1
                            rev_match_diff[a] += 1
                            rev_match_diff[b] -= 1
                        else:
                            fwd_depth_diff[a] += 1
                            fwd_depth_diff[b] -= 1
                            fwd_match_diff[a] += 1
                            fwd_match_diff[b] -= 1
                    last_rpos = rpos + length - 1
                    prev_matched_pos = last_rpos if (start <= last_rpos < end) else -1
                    rpos += length
                    qpos += length
                elif op == 8:  # X (sequence mismatch)
                    a = max(0, rpos - start)
                    b = min(chunk_len, rpos + length - start)
                    if a < b:
                        depth_diff[a] += 1
                        depth_diff[b] -= 1
                        mismatch_diff[a] += 1
                        mismatch_diff[b] -= 1
                        if is_rev:
                            rev_depth_diff[a] += 1
                            rev_depth_diff[b] -= 1
                        else:
                            fwd_depth_diff[a] += 1
                            fwd_depth_diff[b] -= 1
                    prev_matched_pos = -1
                    rpos += length
                    qpos += length
                elif op == 2:  # D (deletion)
                    a = max(0, rpos - start)
                    b = min(chunk_len, rpos + length - start)
                    if a < b:
                        depth_diff[a] += 1
                        depth_diff[b] -= 1
                        del_diff[a] += 1
                        del_diff[b] -= 1
                        if is_rev:
                            rev_depth_diff[a] += 1
                            rev_depth_diff[b] -= 1
                        else:
                            fwd_depth_diff[a] += 1
                            fwd_depth_diff[b] -= 1
                    prev_matched_pos = -1
                    rpos += length
                elif op == 1:  # I (insertion)
                    anchor = rpos - 1
                    if start <= anchor < end:
                        l_idx = anchor - start
                        ins_arr[l_idx] += 1
                        if prev_matched_pos == anchor:
                            match_diff[l_idx] -= 1
                            match_diff[l_idx + 1] += 1
                            if is_rev:
                                rev_match_diff[l_idx] -= 1
                                rev_match_diff[l_idx + 1] += 1
                            else:
                                fwd_match_diff[l_idx] -= 1
                                fwd_match_diff[l_idx + 1] += 1
                            prev_matched_pos = -1
                    qpos += length
                elif op == 0:  # M (fallback for external BAMs without --eqx)
                    a = max(0, rpos - start)
                    b = min(chunk_len, rpos + length - start)
                    if a < b:
                        depth_diff[a] += 1
                        depth_diff[b] -= 1
                        if is_rev:
                            rev_depth_diff[a] += 1
                            rev_depth_diff[b] -= 1
                        else:
                            fwd_depth_diff[a] += 1
                            fwd_depth_diff[b] -= 1

                        q_seq = read.query_sequence
                        if q_seq is not None:
                            if direct_match is None:
                                direct_match = np.zeros(chunk_len, dtype=np.int32)
                                direct_mismatch = np.zeros(chunk_len, dtype=np.int32)
                                direct_fwd_match = np.zeros(chunk_len, dtype=np.int32)
                                direct_rev_match = np.zeros(chunk_len, dtype=np.int32)
                                ref_bytes = np.frombuffer(
                                    ref_seq_chunk.encode("ascii"), dtype=np.uint8
                                )
                            ov_len = b - a
                            q_off = qpos + (start + a - rpos)
                            q_arr = np.frombuffer(
                                q_seq[q_off : q_off + ov_len].upper().encode("ascii"),
                                dtype=np.uint8,
                            )
                            eq_mask = q_arr == ref_bytes[a:b]
                            direct_match[a:b] += eq_mask
                            direct_mismatch[a:b] += ~eq_mask
                            if is_rev:
                                direct_rev_match[a:b] += eq_mask
                            else:
                                direct_fwd_match[a:b] += eq_mask
                            last_rpos = rpos + length - 1
                            prev_matched_pos = (
                                last_rpos
                                if (start <= last_rpos < end and bool(eq_mask[-1]))
                                else -1
                            )
                        else:
                            prev_matched_pos = -1
                    else:
                        prev_matched_pos = -1
                    rpos += length
                    qpos += length
                elif op == 3:  # N (reference skip)
                    prev_matched_pos = -1
                    rpos += length
                elif op == 4:  # S (soft clip)
                    prev_matched_pos = -1
                    qpos += length
                else:  # H (5), P (6), etc.
                    prev_matched_pos = -1

                if rpos > end:
                    break

    depths = np.cumsum(depth_diff[:-1], dtype=np.int32)
    matches = np.cumsum(match_diff[:-1], dtype=np.int32)
    mismatches = np.cumsum(mismatch_diff[:-1], dtype=np.int32)
    deletions = np.cumsum(del_diff[:-1], dtype=np.int32)
    fwd_depths = np.cumsum(fwd_depth_diff[:-1], dtype=np.int32)
    fwd_matches = np.cumsum(fwd_match_diff[:-1], dtype=np.int32)
    rev_depths = np.cumsum(rev_depth_diff[:-1], dtype=np.int32)
    rev_matches = np.cumsum(rev_match_diff[:-1], dtype=np.int32)

    if direct_match is not None:
        matches += direct_match
        mismatches += direct_mismatch
        fwd_matches += direct_fwd_match
        rev_matches += direct_rev_match

    return (
        contig_id,
        start,
        end,
        depths,
        matches,
        mismatches,
        ins_arr,
        deletions,
        fwd_depths,
        fwd_matches,
        rev_depths,
        rev_matches,
    )


_multi_thread = _profile_chunk


class BaseAccuracy:
    """Profiles base-level read concordance, structural masks, and enriched discordant k-mers."""

    def __init__(
        self,
        bam_path: str,
        fasta_path: str,
        threads: int = 4,
        k_list: list[int] | None = None,
        min_concordance: float = 0.85,
        majority_threshold: float = 0.50,
        depth_variance_limit: float = 2.5,
        merge_gap_bp: int = 500,
        min_structural_block_bp: int = 20,
        max_depth: int = 100000,
        canonical: bool = True,
        records: list[SeqRecord] | None = None,
        linear_mode: bool = False,
        is_cyclic: bool | None = None,
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
        self.linear_mode = linear_mode
        self.is_cyclic = is_cyclic

        raw_records = records if records is not None else list(SeqIO.parse(self.fasta_path, "fasta"))
        self.records = [rec for rec in raw_records if len(rec.seq) > 0]
        if not self.records:
            raise ValueError(f"No non-empty FASTA records found in {self.fasta_path}")

        self.primary_contig = max(self.records, key=lambda r: len(r.seq))
        self.seq_len = len(self.primary_contig.seq)

        self.original_contig_id = self.primary_contig.id
        self.bam_contig_id = (
            self.original_contig_id if self.linear_mode else f"{self.original_contig_id}_cyclic"
        )
        self.ref_seq = str(self.primary_contig.seq).upper()

        self.masked_bases = 0
        self.masked_pct = 0.0
        self.masked_regions: list[dict[str, Any]] = []
        self.low_concordance_blocks: list[dict[str, Any]] = []
        self.all_masked_regions: list[dict[str, Any]] = []
        self.consensus_error_sites: list[dict[str, Any]] = []

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

        self.total_mismatches = 0
        self.total_insertions = 0
        self.total_deletions = 0
        self.snp_rate = 0.0
        self.indel_rate = 0.0
        self.consensus_error_snp_sites = 0
        self.consensus_error_indel_sites = 0
        self.read_strand_biased_sites = 0

        self.contig_base_metrics: dict[str, dict[str, Any]] = {}
        self.raw_motif_hits: list[dict[str, Any]] = []
        self.explained_discordant_bases = 0
        self.explained_discordant_pct = 0.0

        self._primary_arrays: dict[str, np.ndarray] = {}

    def _resolve_bam_contig_and_buffer(
        self, contig_id: str, orig_len: int, bam_refs: dict[str, int] | None = None
    ) -> tuple[str, int]:
        if self.linear_mode or self.is_cyclic is False:
            return contig_id, 0

        cyclic_name = f"{contig_id}_cyclic"
        if bam_refs is not None:
            if cyclic_name in bam_refs:
                bam_len = bam_refs[cyclic_name]
                buf = min(bam_len - orig_len, orig_len) if bam_len > orig_len else 0
                return cyclic_name, buf
            if contig_id in bam_refs:
                return contig_id, 0
            raise ValueError(f"Contig {cyclic_name} or {contig_id} not found in BAM")

        buf = (
            self._get_buffer_length()
            if contig_id == self.original_contig_id
            else self._get_buffer_length_for(cyclic_name, orig_len)
        )
        return cyclic_name, buf

    def _get_buffer_length_for(self, bam_contig_id: str, orig_len: int) -> int:
        if self.linear_mode or self.is_cyclic is False:
            return 0
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

    def _get_bam_ref_lengths(self) -> dict[str, int] | None:
        try:
            with pysam.AlignmentFile(self.bam_path, "rb") as bam:
                refs = getattr(bam, "references", None)
                lengths = getattr(bam, "lengths", None)
                if isinstance(refs, (list, tuple)) and isinstance(lengths, (list, tuple)):
                    return dict(zip(refs, lengths))
        except Exception:
            pass
        return None

    def _summarize_contig(
        self,
        contig_id: str,
        seq_str: str,
        buffer_len: int,
        raw_arrays: dict[str, np.ndarray],
    ) -> tuple[np.ndarray, np.ndarray, dict[str, Any], dict[str, np.ndarray]]:
        seq_len = len(seq_str)
        folded: dict[str, np.ndarray] = {}
        for key, arr in raw_arrays.items():
            base_arr = arr[:seq_len].copy()
            if buffer_len > 0:
                base_arr[:buffer_len] += arr[seq_len : seq_len + buffer_len]
            folded[key] = base_arr

        depth_arr = folded["depth"]
        match_arr = folded["match"]
        mismatch_arr = folded["mismatch"]
        ins_arr = folded["ins"]
        del_arr = folded["del"]
        fwd_depth_arr = folded["fwd_depth"]
        fwd_match_arr = folded["fwd_match"]
        rev_depth_arr = folded["rev_depth"]
        rev_match_arr = folded["rev_match"]

        concord_arr = np.ones(seq_len, dtype=float)
        valid_mask = depth_arr > 0
        concord_arr[valid_mask] = match_arr[valid_mask] / depth_arr[valid_mask]
        concord_arr[~valid_mask] = 0.0

        non_zero_depths = depth_arr[valid_mask]
        median_depth = float(np.median(non_zero_depths)) if non_zero_depths.size > 0 else 1.0
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

        # Classify masked bases into distinct categories so merging does not mix or dilute blocks
        cat_arr = np.zeros(seq_len, dtype=np.int8)
        cat_arr[contiguous_discord_mask] = 3
        cat_arr[depth_arr < lower_depth] = 2
        cat_arr[depth_arr > upper_depth] = 1

        raw_blocks: list[tuple[int, int, int]] = []
        for cat_code in (1, 2, 3):
            cat_mask = cat_arr == cat_code
            if not np.any(cat_mask):
                continue
            edges = np.diff(np.concatenate(([0], cat_mask.view(np.int8), [0])))
            starts = np.where(edges == 1)[0]
            ends = np.where(edges == -1)[0]
            for s, e in zip(starts, ends):
                if (e - s) >= self.min_structural_block_bp:
                    raw_blocks.append((int(s), int(e), cat_code))

        raw_blocks.sort(key=lambda x: x[0])
        merged_blocks: list[list[int]] = []
        for s, e, cat_code in raw_blocks:
            if (
                merged_blocks
                and merged_blocks[-1][2] == cat_code
                and (s - merged_blocks[-1][1]) < self.merge_gap_bp
            ):
                merged_blocks[-1][1] = e
            else:
                merged_blocks.append([s, e, cat_code])

        all_intervals: list[dict[str, Any]] = []
        low_concord_blocks: list[dict[str, Any]] = []

        for s, e, cat_code in merged_blocks:
            if (e - s) < self.min_structural_block_bp:
                continue
            masked_slice = cat_arr[s:e] == cat_code
            if not np.any(masked_slice):
                masked_slice = structural_mask[s:e]
            mean_depth = float(np.mean(depth_arr[s:e][masked_slice]))
            mean_concord = float(np.mean(concord_arr[s:e][masked_slice])) * 100.0
            rel_depth = mean_depth / median_depth if median_depth > 0 else 0.0

            if cat_code == 1 or mean_depth > upper_depth:
                category = "HIGH_DEPTH"
                reason = f"High Depth ({rel_depth:.1f}x of median; Collapsed Repeat / Multi-copy Operon)"
            elif cat_code == 2 and mean_depth < 10.0:
                category = "LOW_DEPTH"
                reason = "Low Depth (<10x; Coverage Gap / Dropout)"
            elif cat_code == 2:
                category = "LOW_DEPTH"
                reason = (
                    f"Relative Depth Dip ({rel_depth:.2f}x of median; "
                    f"Structural Heterogeneity / Prophage / Mappability)"
                )
            else:
                category = "LOW_CONCORDANCE"
                reason = "Contiguous Low Concordance (Unpolished Block / Structural Variant)"

            block_info = {
                "contig": contig_id,
                "start": int(s),
                "end": int(e),
                "length": int(e - s),
                "mean_depth": round(mean_depth, 1),
                "mean_concordance": round(mean_concord, 1),
                "reason": reason,
                "category": category,
            }
            all_intervals.append(block_info)
            if mean_concord < 90.0:
                low_concord_blocks.append(block_info)

        self.all_masked_regions.extend(all_intervals)

        # Keep ALL low-concordance blocks plus top 50 longest remaining regions for reporting
        low_concord_blocks.sort(key=lambda x: x["length"], reverse=True)
        other_blocks = [b for b in all_intervals if b["mean_concordance"] >= 90.0]
        other_blocks.sort(key=lambda x: x["length"], reverse=True)
        reported_regions = low_concord_blocks + other_blocks[:50]
        reported_regions.sort(key=lambda x: x["length"], reverse=True)

        # SNP vs Indel & Read-Strand Stratification across unmasked evaluated bases
        eval_read_bases = int(np.sum(depth_arr[unmasked_flags]))
        tot_mm = int(np.sum(mismatch_arr[unmasked_flags]))
        tot_ins = int(np.sum(ins_arr[unmasked_flags]))
        tot_del = int(np.sum(del_arr[unmasked_flags]))
        snp_rate = (tot_mm / eval_read_bases) if eval_read_bases > 0 else 0.0
        indel_rate = ((tot_ins + tot_del) / eval_read_bases) if eval_read_bases > 0 else 0.0

        cons_snp_sites = 0
        cons_indel_sites = 0
        strand_biased_sites = 0

        err_positions = np.flatnonzero(majority_err_flags)
        for pos_idx in err_positions:
            pos = int(pos_idx)
            d_val = int(depth_arr[pos])
            m_val = int(match_arr[pos])
            mm_val = int(mismatch_arr[pos])
            ins_val = int(ins_arr[pos])
            del_val = int(del_arr[pos])
            conc = float(concord_arr[pos])

            indel_val = ins_val + del_val
            if mm_val == 0 and indel_val == 0:
                dom_type = "OTHER"
            elif mm_val >= indel_val and mm_val > 0 and indel_val == 0:
                dom_type = "SNP"
                cons_snp_sites += 1
            elif indel_val > 0 and mm_val == 0:
                dom_type = "INS" if ins_val >= del_val else "DEL"
                cons_indel_sites += 1
            else:
                dom_type = "MIXED"
                if mm_val >= indel_val:
                    cons_snp_sites += 1
                else:
                    cons_indel_sites += 1

            fd = int(fwd_depth_arr[pos])
            fm = int(fwd_match_arr[pos])
            rd = int(rev_depth_arr[pos])
            rm = int(rev_match_arr[pos])
            f_conc = (fm / fd) if fd > 0 else 0.0
            r_conc = (rm / rd) if rd > 0 else 0.0

            if fd > 0 and rd > 0:
                table = [[fm, max(0, fd - fm)], [rm, max(0, rd - rm)]]
                _, s_pval = fisher_exact(table, alternative="two-sided")
                s_pval = float(s_pval)
            else:
                s_pval = 1.0

            if s_pval < 0.01:
                strand_biased_sites += 1

            ctx_start = max(0, pos - 5)
            ctx_end = min(seq_len, pos + 6)
            self.consensus_error_sites.append(
                {
                    "contig": contig_id,
                    "pos_1based": pos + 1,
                    "ref_base": seq_str[pos],
                    "depth": d_val,
                    "concordance": round(conc, 4),
                    "matches": m_val,
                    "mismatches": mm_val,
                    "insertions": ins_val,
                    "deletions": del_val,
                    "dominant_error_type": dom_type,
                    "fwd_depth": fd,
                    "fwd_concordance": round(f_conc, 4),
                    "rev_depth": rd,
                    "rev_concordance": round(r_conc, 4),
                    "strand_pvalue": round(s_pval, 6),
                    "context": seq_str[ctx_start:ctx_end],
                }
            )

        qv_val, qv_capped = compute_empirical_qv(consensus_error_bases, evaluated_bases)
        cons_per_mb = (
            consensus_error_bases / (evaluated_bases / 1_000_000) if evaluated_bases > 0 else 0.0
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
            "masked_regions": reported_regions,
            "total_mismatches": tot_mm,
            "total_insertions": tot_ins,
            "total_deletions": tot_del,
            "snp_rate": round(snp_rate, 6),
            "indel_rate": round(indel_rate, 6),
            "consensus_error_snp_sites": cons_snp_sites,
            "consensus_error_indel_sites": cons_indel_sites,
            "read_strand_biased_sites": strand_biased_sites,
        }
        return concord_arr, structural_mask, metrics, folded

    def build_fidelity_profile(self) -> tuple[np.ndarray, np.ndarray]:
        logger.info(
            " > [Fidelity] Profiling base concordance across %d contig(s) "
            "(primary: %s, %d bp) using %d thread(s)...",
            len(self.records),
            self.original_contig_id,
            self.seq_len,
            self.threads,
        )

        self.all_masked_regions = []
        self.consensus_error_sites = []
        self.contig_base_metrics = {}

        bam_refs = self._get_bam_ref_lengths()
        contig_plans: list[tuple[str, str, str, int, int]] = []
        all_tasks: list[tuple[str, str, str, int, int, str, int]] = []
        raw_buffers: dict[str, dict[str, np.ndarray]] = {}

        for rec in self.records:
            cid = rec.id
            seq_str = str(rec.seq).upper()
            seq_len = len(seq_str)
            try:
                bam_cid, buffer_len = self._resolve_bam_contig_and_buffer(
                    cid, seq_len, bam_refs=bam_refs
                )
            except ValueError:
                if cid == self.original_contig_id:
                    raise
                continue

            total_profile_len = seq_len + buffer_len
            cyclic_ref_seq = seq_str + seq_str[:buffer_len]
            contig_plans.append((cid, seq_str, bam_cid, buffer_len, total_profile_len))

            raw_buffers[cid] = {
                "depth": np.zeros(total_profile_len, dtype=np.int32),
                "match": np.zeros(total_profile_len, dtype=np.int32),
                "mismatch": np.zeros(total_profile_len, dtype=np.int32),
                "ins": np.zeros(total_profile_len, dtype=np.int32),
                "del": np.zeros(total_profile_len, dtype=np.int32),
                "fwd_depth": np.zeros(total_profile_len, dtype=np.int32),
                "fwd_match": np.zeros(total_profile_len, dtype=np.int32),
                "rev_depth": np.zeros(total_profile_len, dtype=np.int32),
                "rev_match": np.zeros(total_profile_len, dtype=np.int32),
            }

            if total_profile_len <= _DEFAULT_CHUNK_BP or self.threads == 1:
                n_chunks = 1
            else:
                n_chunks = max(
                    self.threads, math.ceil(total_profile_len / _DEFAULT_CHUNK_BP)
                )
            chunk_size = max(1, math.ceil(total_profile_len / n_chunks))
            for i in range(n_chunks):
                start = i * chunk_size
                end = min((i + 1) * chunk_size, total_profile_len)
                if start >= end:
                    break
                all_tasks.append(
                    (
                        cid,
                        self.bam_path,
                        bam_cid,
                        start,
                        end,
                        cyclic_ref_seq[start:end],
                        self.max_depth,
                    )
                )

        if not all_tasks:
            raise ValueError("No alignment tasks generated for assembly contigs.")

        if self.threads <= 1 or len(all_tasks) <= 1:
            results = [_profile_chunk(t) for t in all_tasks]
        else:
            with mp.Pool(processes=min(self.threads, len(all_tasks))) as pool:
                results = pool.map(_profile_chunk, all_tasks)

        for res in results:
            cid, start, end, d, m, mm, ins_c, del_c, fd, fm, rd, rm = res
            buf = raw_buffers[cid]
            buf["depth"][start:end] = d
            buf["match"][start:end] = m
            buf["mismatch"][start:end] = mm
            buf["ins"][start:end] = ins_c
            buf["del"][start:end] = del_c
            buf["fwd_depth"][start:end] = fd
            buf["fwd_match"][start:end] = fm
            buf["rev_depth"][start:end] = rd
            buf["rev_match"][start:end] = rm

        primary_concord_arr = None
        primary_structural_mask = None

        for cid, seq_str, _, buffer_len, _ in contig_plans:
            c_arr, s_mask, m_dict, folded = self._summarize_contig(
                cid, seq_str, buffer_len, raw_buffers[cid]
            )
            self.contig_base_metrics[cid] = m_dict
            if cid == self.original_contig_id:
                primary_concord_arr = c_arr
                primary_structural_mask = s_mask
                self._primary_arrays = folded

        assert primary_concord_arr is not None and primary_structural_mask is not None
        p_metrics = self.contig_base_metrics[self.original_contig_id]

        self.masked_bases = p_metrics["masked_bases"]
        self.masked_pct = (self.masked_bases / self.seq_len) * 100.0 if self.seq_len > 0 else 0.0
        self.evaluated_bases = p_metrics["evaluated_bases"]
        self.discordant_bases = p_metrics["discordant_bases"]
        self.consensus_error_bases = p_metrics["consensus_error_bases"]
        self.subconsensus_bias_bases = p_metrics["subconsensus_bias_bases"]
        self.baseline_discordance_rate = (
            self.discordant_bases / self.evaluated_bases if self.evaluated_bases > 0 else 0.0
        )
        self.consensus_error_rate = (
            self.consensus_error_bases / self.evaluated_bases if self.evaluated_bases > 0 else 0.0
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

        self.total_mismatches = p_metrics["total_mismatches"]
        self.total_insertions = p_metrics["total_insertions"]
        self.total_deletions = p_metrics["total_deletions"]
        self.snp_rate = p_metrics["snp_rate"]
        self.indel_rate = p_metrics["indel_rate"]
        self.consensus_error_snp_sites = p_metrics["consensus_error_snp_sites"]
        self.consensus_error_indel_sites = p_metrics["consensus_error_indel_sites"]
        self.read_strand_biased_sites = p_metrics["read_strand_biased_sites"]

        qv_str = f">{self.empirical_qv:.1f}" if self.empirical_qv_capped else f"{self.empirical_qv:.1f}"
        logger.info("   - Primary Median Depth: %.1fx", p_metrics["median_depth"])
        logger.info("   - Structurally Masked: %d bp (%.2f%%)", self.masked_bases, self.masked_pct)
        logger.info(
            "   - Evaluated Positions: %d bp | Total Discordant (<%.0f%%): %d bp (%.4f%%)",
            self.evaluated_bases,
            self.min_concordance * 100,
            self.discordant_bases,
            self.baseline_discordance_rate * 100,
        )
        logger.info(
            "     * Consensus Errors (<%.0f%% support): %d bp (%.2f/Mb, Concordance QV ≈ %s) | "
            "Sub-consensus Read Bias (%.0f%%–%.0f%%): %d bp (%.1f/Mb)",
            self.majority_threshold * 100,
            self.consensus_error_bases,
            self.consensus_errors_per_mb,
            qv_str,
            self.majority_threshold * 100,
            self.min_concordance * 100,
            self.subconsensus_bias_bases,
            self.subconsensus_bias_per_mb,
        )

        return primary_concord_arr, primary_structural_mask

    def analyze_motifs(
        self, min_occurrences: int = 100, min_discordant_count: int = 10
    ) -> list[dict[str, Any]]:
        """
        Vectorized scan of unmasked regions on the primary chromosome for k-mer contexts
        enriched for read-to-assembly discordance (< min_concordance).
        """
        concord_arr, structural_mask = self.build_fidelity_profile()
        logger.info(" > [Fidelity] Testing systematic errors in %s-mer motifs...", self.k_list)

        if self.evaluated_bases == 0 or self.discordant_bases == 0:
            self.raw_motif_hits = []
            self.explained_discordant_bases = 0
            self.explained_discordant_pct = 0.0
            return []

        global_error_rate = self.baseline_discordance_rate

        raw_bytes = np.frombuffer(self.ref_seq.encode("ascii", errors="replace"), dtype=np.uint8)
        codes = _BASE_LUT[raw_bytes]
        bad_pos = structural_mask | (codes >= 4)
        bad_cumsum = np.zeros(self.seq_len + 1, dtype=np.int32)
        np.cumsum(bad_pos.astype(np.int32), out=bad_cumsum[1:])
        c_clean = np.where(codes < 4, codes, 0).astype(np.uint32)

        disc_mask = (~structural_mask) & (concord_arr < self.min_concordance)
        maj_mask = (~structural_mask) & (concord_arr < self.majority_threshold)

        candidate_tests: list[dict[str, Any]] = []

        for k in self.k_list:
            if self.seq_len < k:
                continue

            num_windows = self.seq_len - k + 1
            valid_win = (bad_cumsum[k:] - bad_cumsum[:num_windows]) == 0
            if not np.any(valid_win):
                continue

            valid_starts = np.flatnonzero(valid_win)
            kmer_ids = np.zeros(num_windows, dtype=np.uint32)
            for j in range(k):
                kmer_ids = (kmer_ids << 2) | c_clean[j : j + num_windows]
            valid_kids = kmer_ids[valid_starts]

            if k <= 7:
                canon_lut, is_fwd_lut, is_pal_lut, labels = _build_kmer_lut(k, self.canonical)
                mapped_canon = canon_lut[valid_kids]
                mapped_fwd = is_fwd_lut[valid_kids]
                mapped_pal = is_pal_lut[valid_kids]
                n_canon = len(labels)
            else:
                # Sparse path for k >= 8
                uniq_kids, inv = np.unique(valid_kids, return_inverse=True)
                labels_map: dict[int, int] = {}
                labels = []
                u_canon = np.empty(len(uniq_kids), dtype=np.int32)
                u_fwd = np.empty(len(uniq_kids), dtype=bool)
                u_pal = np.empty(len(uniq_kids), dtype=bool)
                for u_i, kid in enumerate(uniq_kids):
                    kid_int = int(kid)
                    if not self.canonical:
                        labels_map[kid_int] = u_i
                        labels.append(_decode_kmer(kid_int, k))
                        u_canon[u_i] = u_i
                        u_fwd[u_i] = True
                        u_pal[u_i] = False
                    else:
                        rc_id = _revcomp_kmer_id(kid_int, k)
                        is_pal = kid_int == rc_id
                        is_fwd = kid_int <= rc_id
                        cid = kid_int if is_fwd else rc_id
                        c_idx = labels_map.get(cid)
                        if c_idx is None:
                            c_idx = len(labels)
                            labels_map[cid] = c_idx
                            labels.append(
                                f"{_decode_kmer(cid, k)}/{_decode_kmer(_revcomp_kmer_id(cid, k), k)}"
                            )
                        u_canon[u_i] = c_idx
                        u_fwd[u_i] = is_pal or is_fwd
                        u_pal[u_i] = is_pal
                mapped_canon = u_canon[inv]
                mapped_fwd = u_fwd[inv]
                mapped_pal = u_pal[inv]
                n_canon = len(labels)

            motif_totals = np.bincount(mapped_canon, minlength=n_canon)
            eligible_mask = motif_totals >= min_occurrences
            if not np.any(eligible_mask):
                continue

            fwd_totals = np.bincount(mapped_canon[mapped_fwd], minlength=n_canon)
            rev_totals = motif_totals - fwd_totals

            mid_hi = k // 2
            mid_lo = (k // 2) - 1
            even_pal_arr = mapped_pal & ((k % 2) == 0)

            if not self.canonical:
                primary_offsets = np.full(len(valid_starts), mid_hi, dtype=np.int32)
            else:
                primary_offsets = np.where(
                    mapped_pal | mapped_fwd, mid_hi, k - 1 - mid_hi
                )
            pos_primary = valid_starts + primary_offsets
            hit_primary = disc_mask[pos_primary]

            if np.any(even_pal_arr):
                pos_secondary = valid_starts + mid_lo
                hit_secondary = even_pal_arr & disc_mask[pos_secondary]
                hit_any = hit_primary | hit_secondary
            else:
                pos_secondary = None
                hit_secondary = None
                hit_any = hit_primary

            disc_canon = mapped_canon[hit_any]
            if disc_canon.size == 0:
                continue

            motif_discordant = np.bincount(disc_canon, minlength=n_canon)
            candidate_ids = np.flatnonzero(eligible_mask & (motif_discordant >= min_discordant_count))
            if candidate_ids.size == 0:
                continue

            candidate_set = set(candidate_ids.tolist())
            disc_fwd = mapped_fwd[hit_any]
            motif_fwd_disc = np.bincount(disc_canon[disc_fwd], minlength=n_canon)
            motif_rev_disc = motif_discordant - motif_fwd_disc

            motif_positions: dict[int, set[int]] = defaultdict(set)
            motif_majority_positions: dict[int, set[int]] = defaultdict(set)

            hit_indices = np.flatnonzero(hit_any)
            for w_idx in hit_indices:
                c_idx = int(mapped_canon[w_idx])
                if c_idx not in candidate_set:
                    continue
                if hit_primary[w_idx]:
                    p1 = int(pos_primary[w_idx])
                    motif_positions[c_idx].add(p1)
                    if maj_mask[p1]:
                        motif_majority_positions[c_idx].add(p1)
                if hit_secondary is not None and hit_secondary[w_idx]:
                    p2 = int(pos_secondary[w_idx])
                    motif_positions[c_idx].add(p2)
                    if maj_mask[p2]:
                        motif_majority_positions[c_idx].add(p2)

            for c_idx in candidate_ids:
                c_idx_int = int(c_idx)
                total = int(motif_totals[c_idx_int])
                discordant_count = int(motif_discordant[c_idx_int])
                kmer_key = labels[c_idx_int]
                is_pal = False
                if self.canonical:
                    f_part, r_part = kmer_key.split("/")
                    is_pal = f_part == r_part

                # Even-length palindromes test two central bases per genomic occurrence
                p_null = (
                    1.0 - (1.0 - global_error_rate) ** 2
                    if (is_pal and (k % 2) == 0)
                    else global_error_rate
                )
                #p_null = min(max(p_null, 1e-12), 1.0 - 1e-12)

                test = binomtest(
                    k=discordant_count,
                    n=total,
                    p=p_null,
                    alternative="greater",
                )
                error_rate = discordant_count / total
                fold_increase = error_rate / p_null
                lb_rate = wilson_lower_bound(discordant_count, total)
                ci_lower_fold = lb_rate / p_null

                f_disc = int(motif_fwd_disc[c_idx_int])
                r_disc = int(motif_rev_disc[c_idx_int])
                f_tot = int(fwd_totals[c_idx_int])
                r_tot = int(rev_totals[c_idx_int])

                orientation_state = "Undetermined"
                strand_biased = False
                strand_pval = 1.0
                if is_pal:
                    orientation_state = "Palindromic"
                elif (f_tot + r_tot) > 0 and discordant_count > 0:
                    p_fwd = f_tot / (f_tot + r_tot)
                    strand_pval = float(
                        binomtest(
                            k=f_disc, n=discordant_count, p=p_fwd, alternative="two-sided"
                        ).pvalue
                    )
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

                pos_set = motif_positions[c_idx_int]
                maj_set = motif_majority_positions[c_idx_int]
                maj_count = len(maj_set)

                candidate_tests.append(
                    {
                        "kmer_size": k,
                        "motif": kmer_key,
                        "is_palindrome": is_pal,
                        "occurrences": total,
                        "discordant_count": discordant_count,
                        "consensus_error_count": maj_count,
                        "subconsensus_bias_count": max(0, discordant_count - maj_count),
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
                        "p_value": float(test.pvalue),
                        "_positions": pos_set,
                        "_maj_positions": maj_set,
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

        logger.info(
            "   - Significant k-mers: %d -> Clustered into %d Motif Family(ies) explaining "
            "%d/%d (%.1f%%) of discordant sites.",
            len(significant_hits),
            len(families),
            self.explained_discordant_bases,
            self.discordant_bases,
            self.explained_discordant_pct,
        )
        return families

    def _compute_coord_stratification(self, coords: set[int]) -> dict[str, Any]:
        """Aggregate SNP vs indel counts and read-strand concordance over a set of coordinates."""
        if not coords or not self._primary_arrays:
            return {
                "snp_discordant": 0,
                "indel_discordant": 0,
                "read_fwd_depth": 0,
                "read_fwd_match": 0,
                "read_rev_depth": 0,
                "read_rev_match": 0,
                "read_strand_pvalue": 1.0,
                "read_strand_biased": False,
            }

        idx = np.fromiter(coords, dtype=np.intp)
        mm = self._primary_arrays["mismatch"][idx]
        ins_c = self._primary_arrays["ins"][idx]
        del_c = self._primary_arrays["del"][idx]
        fd = int(np.sum(self._primary_arrays["fwd_depth"][idx]))
        fm = int(np.sum(self._primary_arrays["fwd_match"][idx]))
        rd = int(np.sum(self._primary_arrays["rev_depth"][idx]))
        rm = int(np.sum(self._primary_arrays["rev_match"][idx]))

        indel_tot = ins_c + del_c
        snp_sites = int(np.sum((mm > 0) & (mm >= indel_tot)))
        indel_sites = int(np.sum(indel_tot > mm))

        if fd > 0 and rd > 0:
            table = [[fm, max(0, fd - fm)], [rm, max(0, rd - rm)]]
            _, pval = fisher_exact(table, alternative="two-sided")
            pval = float(pval)
            f_rate = fm / fd
            r_rate = rm / rd
            strand_biased = pval < 0.01 and abs(f_rate - r_rate) >= 0.15
        else:
            pval = 1.0
            strand_biased = False

        return {
            "snp_discordant": snp_sites,
            "indel_discordant": indel_sites,
            "read_fwd_depth": fd,
            "read_fwd_match": fm,
            "read_rev_depth": rd,
            "read_rev_match": rm,
            "read_strand_pvalue": round(pval, 6),
            "read_strand_biased": strand_biased,
        }

    def _cluster_into_families(self, significant_hits: list[dict[str, Any]]) -> list[dict[str, Any]]:
        sorted_hits = sorted(
            significant_hits,
            key=lambda x: (x["p_value"], -x["ci_lower_fold"], -x["discordant_count"]),
        )

        families: list[dict[str, Any]] = []
        for hit in sorted_hits:
            hit_pos = hit["_positions"]
            hit_maj_pos = hit["_maj_positions"]
            placed = False

            for fam in families:
                lead_pos = fam["_lead_positions"]
                min_size = min(len(hit_pos), len(lead_pos))
                coord_overlap = (len(hit_pos & lead_pos) / min_size) if min_size > 0 else 0.0

                if coord_overlap >= 0.25 or _share_core_or_substring(hit["motif"], fam["motif"]):
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
                        "raw_consensus_error_count": hit["consensus_error_count"],
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

        families.sort(
            key=lambda x: (len(x["_union_positions"]), x["ci_lower_fold"]),
            reverse=True,
        )

        # Partition BOTH discordant and consensus-error coordinates mutually exclusively
        claimed_coords: set[int] = set()
        claimed_maj_coords: set[int] = set()
        for fam in families:
            raw_union = fam.pop("_union_positions")
            maj_union = fam.pop("_union_maj_positions")
            fam.pop("_lead_positions", None)

            exclusive_coords = raw_union - claimed_coords
            claimed_coords.update(raw_union)

            exclusive_maj_coords = maj_union - claimed_maj_coords
            claimed_maj_coords.update(maj_union)

            raw_count = len(raw_union)
            assigned_count = len(exclusive_coords)
            raw_maj_count = len(maj_union)
            assigned_maj_count = len(exclusive_maj_coords)

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
            fam["consensus_error_count"] = assigned_maj_count
            fam["raw_consensus_error_count"] = raw_maj_count
            fam["subconsensus_bias_count"] = max(0, assigned_count - assigned_maj_count)

            fam.update(self._compute_coord_stratification(exclusive_coords))

            related = fam["member_motifs"][1:6]
            fam["related_motifs_str"] = ", ".join(related) if related else "-"

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
        families = [fam for fam in families if fam["assigned_unique_sites"] > 0]

        families.sort(
            key=lambda x: (
                x["assigned_unique_sites"],
                x["raw_family_sites"],
                x["ci_lower_fold"],
            ),
            reverse=True,
        )
        return families
