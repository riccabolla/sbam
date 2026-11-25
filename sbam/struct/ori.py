import numpy as np
from Bio import SeqIO
from scipy.stats import poisson

_DNAA_FWD = "TTATCCACA"
_DNAA_REV = "TGTGGATAA"


def _find_dnaa_box_positions(seq: str, max_mismatches: int = 1) -> list[int]:
    """Returns 0-based start coordinates of 9-bp bacterial DnaA boxes (<=1 mismatch) across seq."""
    positions = []
    s_len = len(seq)
    for i in range(s_len - 9 + 1):
        sub = seq[i : i + 9]
        mm_f = sum(1 for a, b in zip(sub, _DNAA_FWD) if a != b)
        if mm_f <= max_mismatches:
            positions.append(i)
            continue
        mm_r = sum(1 for a, b in zip(sub, _DNAA_REV) if a != b)
        if mm_r <= max_mismatches:
            positions.append(i)
    return positions


class BioVal:
    def __init__(
        self,
        fasta_path,
        min_contig_len=500000,
        smooth_windows=15,
    ):
        self.fasta_path = fasta_path
        self.min_contig_len = min_contig_len
        self.smooth_windows = smooth_windows
        self.records = list(SeqIO.parse(self.fasta_path, "fasta"))
        self.results = {}

    def calculate_symmetry(self, window_size=1000):
        print(
            f" > Calculating GC Skew, Replichore Symmetry, and DnaA-box enrichment "
            f"(contigs >= {self.min_contig_len // 1000} kb)..."
        )

        for record in self.records:
            seq = str(record.seq).upper()
            seq_len = len(seq)

            if seq_len < self.min_contig_len:
                continue

            gc_skews = []
            for i in range(0, seq_len, window_size):
                window = seq[i : i + window_size]
                g = window.count("G")
                c = window.count("C")
                skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
                gc_skews.append(skew)

            skews_arr = np.array(gc_skews, dtype=float)
            centered_skews = skews_arr - np.mean(skews_arr)
            raw_cum_skew = np.cumsum(centered_skews)

            k = min(self.smooth_windows, max(1, len(raw_cum_skew) // 10))
            if k > 1:
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

            oric_pos = int(oric_window * window_size)
            ter_pos = int(ter_window * window_size)

            # Genome-wide DnaA-box baseline vs +/-25 kb (50 kb) oriC window enrichment
            flank = 25000
            window_bp = 2 * flank
            all_dnaa_positions = _find_dnaa_box_positions(seq, max_mismatches=1)
            total_dnaa_boxes = len(all_dnaa_positions)
            expected_dnaa_50kb = (
                total_dnaa_boxes * (window_bp / seq_len) if seq_len > 0 else 0.0
            )

            dnaa_boxes_near_ori = sum(
                1
                for pos in all_dnaa_positions
                if min(abs(pos - oric_pos), seq_len - abs(pos - oric_pos)) <= flank
            )
            dnaa_fold_enrichment = (
                (dnaa_boxes_near_ori / expected_dnaa_50kb)
                if expected_dnaa_50kb > 0
                else 0.0
            )
            dnaa_poisson_pval = (
                float(poisson.sf(dnaa_boxes_near_ori - 1, expected_dnaa_50kb))
                if expected_dnaa_50kb > 0 and dnaa_boxes_near_ori > 0
                else 1.0
            )

            num_windows = len(cumulative_skew)
            plot_theta = [(i / num_windows) * 360 for i in range(num_windows)]
            plot_r = cumulative_skew.tolist()

            distance = abs(oric_pos - ter_pos)
            if distance > (seq_len / 2):
                distance = seq_len - distance

            symmetry_degrees = (distance / (seq_len / 2)) * 180

            oric_offset_bp = min(oric_pos, seq_len - oric_pos)
            oric_offset_kb = round(oric_offset_bp / 1000.0, 1)
            oric_offset_deg = round((oric_offset_bp / (seq_len / 2)) * 180, 1)

            # Three-tier rotation classification
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
            }

        return self.results