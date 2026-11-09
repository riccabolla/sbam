import numpy as np
from Bio import SeqIO


class BioVal:
    def __init__(self, fasta_path, min_contig_len=500000):
        self.fasta_path = fasta_path
        self.min_contig_len = min_contig_len
        self.records = list(SeqIO.parse(self.fasta_path, "fasta"))
        self.results = {}

    def calculate_symmetry(self, window_size=1000):
        print(
            f" > Calculating GC Skew and Replichore Symmetry (contigs >= {self.min_contig_len // 1000} kb)..."
        )

        for record in self.records:
            seq = str(record.seq).upper()
            seq_len = len(seq)

            # Restrict bidirectional theta-replication geometry to chromosome-scale replicons
            if seq_len < self.min_contig_len:
                continue

            gc_skews = []
            for i in range(0, seq_len, window_size):
                window = seq[i : i + window_size]
                g = window.count("G")
                c = window.count("C")

                if (g + c) == 0:
                    skew = 0.0
                else:
                    skew = (g - c) / (g + c)
                gc_skews.append(skew)

            skews_arr = np.array(gc_skews, dtype=float)

            # Zero-center window skews so cumulative skew is rotation-invariant
            centered_skews = skews_arr - np.mean(skews_arr)
            cumulative_skew = np.cumsum(centered_skews)

            oric_window = int(np.argmin(cumulative_skew))
            ter_window = int(np.argmax(cumulative_skew))

            oric_pos = int(oric_window * window_size)
            ter_pos = int(ter_window * window_size)

            num_windows = len(cumulative_skew)
            plot_theta = [(i / num_windows) * 360 for i in range(num_windows)]
            plot_r = cumulative_skew.tolist()

            # Shortest circular distance between oriC and ter
            distance = abs(oric_pos - ter_pos)
            if distance > (seq_len / 2):
                distance = seq_len - distance

            symmetry_degrees = (distance / (seq_len / 2)) * 180

            # Shortest circular distance from oriC to coordinate 0 bp (rotation check)
            oric_offset_bp = min(oric_pos, seq_len - oric_pos)
            oric_offset_kb = round(oric_offset_bp / 1000.0, 1)
            oric_offset_deg = round((oric_offset_bp / (seq_len / 2)) * 180, 1)

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
                "oric_offset_bp": int(oric_offset_bp),
                "oric_offset_kb": oric_offset_kb,
                "oric_offset_deg": oric_offset_deg,
                "symmetry": round(symmetry_degrees, 1),
                "viability": viability,
                "plot_theta": plot_theta,
                "plot_r": plot_r,
            }

        return self.results