"""Biological replication architecture validation via GC skew and DnaA-box enrichment."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scipy.stats import poisson

logger = logging.getLogger("sbam")

_DNAA_FWD = "TTATCCACA"
_DNAA_REV = "TGTGGATAA"

_BASE_LUT = np.full(256, 4, dtype=np.uint8)
for _ch, _code in ((ord("A"), 0), (ord("C"), 1), (ord("G"), 2), (ord("T"), 3)):
    _BASE_LUT[_ch] = _code
    _BASE_LUT[_ch + 32] = _code


def _encode_kmer9(kmer: str) -> int:
    val = 0
    for ch in kmer:
        val = (val << 2) | int(_BASE_LUT[ord(ch)])
    return val


def _build_dnaa_lut() -> np.ndarray:
    """Precompute a 256 KB boolean lookup table for all 9-mers within 1 mismatch of DnaA boxes."""
    lut = np.zeros(1 << 18, dtype=bool)
    bases = ("A", "C", "G", "T")
    for motif in (_DNAA_FWD, _DNAA_REV):
        lut[_encode_kmer9(motif)] = True
        for pos in range(9):
            for base in bases:
                if base == motif[pos]:
                    continue
                variant = motif[:pos] + base + motif[pos + 1 :]
                lut[_encode_kmer9(variant)] = True
    return lut


_DNAA_LUT = _build_dnaa_lut()


def _find_dnaa_box_positions(
    seq: str, max_mismatches: int = 1, circular: bool = True
) -> list[int]:
    """Return 0-based start coordinates of 9-bp bacterial DnaA boxes across `seq`."""
    s_len = len(seq)
    if s_len < 9:
        return []

    if max_mismatches != 1:
        positions: list[int] = []
        scan_seq = (seq + seq[:8]) if circular else seq
        limit = s_len if circular else (s_len - 9 + 1)
        for i in range(limit):
            sub = scan_seq[i : i + 9]
            mm_f = sum(1 for a, b in zip(sub, _DNAA_FWD) if a != b)
            if mm_f <= max_mismatches:
                positions.append(i)
                continue
            mm_r = sum(1 for a, b in zip(sub, _DNAA_REV) if a != b)
            if mm_r <= max_mismatches:
                positions.append(i)
        return positions

    scan_seq = (seq + seq[:8]) if circular else seq
    raw = np.frombuffer(scan_seq.encode("ascii", errors="replace"), dtype=np.uint8)
    codes = _BASE_LUT[raw]
    invalid = (codes >= 4).astype(np.int32)
    inv_cumsum = np.zeros(len(codes) + 1, dtype=np.int32)
    np.cumsum(invalid, out=inv_cumsum[1:])
    valid = (inv_cumsum[9:] - inv_cumsum[:-9]) == 0

    c = np.where(codes < 4, codes, 0).astype(np.uint32)
    kmer_ids = (
        (c[:-8] << 16)
        | (c[1:-7] << 14)
        | (c[2:-6] << 12)
        | (c[3:-5] << 10)
        | (c[4:-4] << 8)
        | (c[5:-3] << 6)
        | (c[6:-2] << 4)
        | (c[7:-1] << 2)
        | c[8:]
    )
    hits = np.flatnonzero(valid & _DNAA_LUT[kmer_ids])
    return hits.tolist()


class BioVal:
    """Computes smoothed cumulative GC skew, replichore symmetry, and DnaA-box enrichment."""

    def __init__(
        self,
        fasta_path: str,
        min_contig_len: int = 500000,
        smooth_windows: int = 15,
        records: list[SeqRecord] | None = None,
        linear_mode: bool = False,
    ):
        self.fasta_path = fasta_path
        self.min_contig_len = min_contig_len
        self.smooth_windows = smooth_windows
        self.linear_mode = linear_mode
        self.records = records if records is not None else list(SeqIO.parse(self.fasta_path, "fasta"))
        self.results: dict[str, dict[str, Any]] = {}

    def calculate_symmetry(self, window_size: int = 1000) -> dict[str, dict[str, Any]]:
        logger.info(
            " > Calculating GC Skew, Replichore Symmetry, and DnaA-box enrichment "
            "(contigs >= %d kb)...",
            self.min_contig_len // 1000,
        )

        for record in self.records:
            seq = str(record.seq).upper()
            seq_len = len(seq)

            if seq_len < self.min_contig_len:
                continue

            gc_skews: list[float] = []
            weights: list[float] = []
            for i in range(0, seq_len, window_size):
                window = seq[i : i + window_size]
                w_len = len(window)
                g = window.count("G")
                c = window.count("C")
                skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
                gc_skews.append(skew)
                weights.append(w_len / window_size)

            skews_arr = np.array(gc_skews, dtype=float)
            weights_arr = np.array(weights, dtype=float)
            total_weight = float(np.sum(weights_arr))
            mean_skew = float(np.sum(skews_arr * weights_arr) / total_weight) if total_weight > 0 else 0.0
            deltas = (skews_arr - mean_skew) * weights_arr

            if self.linear_mode:
                raw_cum_skew = np.concatenate(([0.0], np.cumsum(deltas)))
                x_coords = [min(seq_len, idx * window_size) for idx in range(len(raw_cum_skew))]
            else:
                raw_cum_skew = np.concatenate(([0.0], np.cumsum(deltas[:-1])))
                x_coords = [idx * window_size for idx in range(len(raw_cum_skew))]

            k = min(self.smooth_windows, max(1, len(raw_cum_skew) // 10))
            if k > 1:
                if self.linear_mode:
                    padded = np.pad(raw_cum_skew, k, mode="edge")
                else:
                    padded = np.concatenate(
                        [raw_cum_skew[-k:], raw_cum_skew, raw_cum_skew[:k]]
                    )
                kernel = np.ones(k, dtype=float) / k
                smoothed_full = np.convolve(padded, kernel, mode="same")
                cumulative_skew = smoothed_full[k : k + len(raw_cum_skew)]
            else:
                cumulative_skew = raw_cum_skew

            oric_window = int(np.argmin(cumulative_skew))
            ter_window = int(np.argmax(cumulative_skew))

            oric_pos = int(x_coords[oric_window])
            ter_pos = int(x_coords[ter_window])

            flank = 25000
            all_dnaa_positions = _find_dnaa_box_positions(
                seq, max_mismatches=1, circular=not self.linear_mode
            )
            total_dnaa_boxes = len(all_dnaa_positions)

            if self.linear_mode:
                win_start = max(0, oric_pos - flank)
                win_end = min(seq_len, oric_pos + flank)
                window_bp = max(1, win_end - win_start)
                dnaa_boxes_near_ori = sum(
                    1 for pos in all_dnaa_positions if abs(pos - oric_pos) <= flank
                )
            else:
                window_bp = min(seq_len, 2 * flank)
                dnaa_boxes_near_ori = sum(
                    1
                    for pos in all_dnaa_positions
                    if min(abs(pos - oric_pos), seq_len - abs(pos - oric_pos)) <= flank
                )

            expected_dnaa_50kb = total_dnaa_boxes * (window_bp / seq_len) if seq_len > 0 else 0.0
            dnaa_fold_enrichment = (
                (dnaa_boxes_near_ori / expected_dnaa_50kb) if expected_dnaa_50kb > 0 else 0.0
            )
            dnaa_poisson_pval = (
                float(poisson.sf(dnaa_boxes_near_ori - 1, expected_dnaa_50kb))
                if expected_dnaa_50kb > 0 and dnaa_boxes_near_ori > 0
                else 1.0
            )

            num_windows = len(cumulative_skew)
            plot_theta = [(i / num_windows) * 360 for i in range(num_windows)]
            plot_r = cumulative_skew.tolist()
            plot_x_kb = [round(pos / 1000.0, 2) for pos in x_coords]

            if self.linear_mode:
                # In linear replicons, ideal bidirectional origin lies near chromosome center
                oric_dist_to_end = min(oric_pos, seq_len - oric_pos)
                symmetry_degrees = (oric_dist_to_end / (seq_len / 2.0)) * 180.0
                oric_offset_bp = oric_dist_to_end
                oric_offset_kb = round(oric_offset_bp / 1000.0, 1)
                oric_offset_deg = round(symmetry_degrees, 1)
                rotation_status = "LINEAR"
            else:
                distance = abs(oric_pos - ter_pos)
                if distance > (seq_len / 2):
                    distance = seq_len - distance
                symmetry_degrees = (distance / (seq_len / 2)) * 180

                oric_offset_bp = min(oric_pos, seq_len - oric_pos)
                oric_offset_kb = round(oric_offset_bp / 1000.0, 1)
                oric_offset_deg = round((oric_offset_bp / (seq_len / 2)) * 180, 1)

                if oric_offset_kb <= 25.0:
                    rotation_status = "ROTATED"
                elif oric_offset_kb <= 200.0 and oric_offset_deg <= 20.0:
                    rotation_status = "NEAR_ORIC"
                else:
                    rotation_status = "UNROTATED"

            if 165 <= symmetry_degrees <= 180:
                viability = "ACCEPTABLE"
            elif 150 <= symmetry_degrees < 165:
                viability = "ATYPICAL"
            else:
                viability = "WARNING"

            self.results[record.id] = {
                "length": seq_len,
                "oric_pos": oric_pos,
                "ter_pos": ter_pos,
                "dnaa_boxes_50kb": dnaa_boxes_near_ori,
                "expected_dnaa_50kb": round(expected_dnaa_50kb, 1),
                "dnaa_fold_enrichment": round(dnaa_fold_enrichment, 2),
                "dnaa_poisson_pval": round(dnaa_poisson_pval, 4),
                "oric_offset_bp": int(oric_offset_bp),
                "oric_offset_kb": oric_offset_kb,
                "oric_offset_deg": oric_offset_deg,
                "rotation_status": rotation_status,
                "symmetry": round(symmetry_degrees, 1),
                "viability": viability,
                "plot_theta": plot_theta,
                "plot_r": plot_r,
                "plot_x_kb": plot_x_kb,
            }

        return self.results
