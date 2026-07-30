import numpy as np
from Bio import SeqIO


class BioVal:
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.records = list(SeqIO.parse(self.fasta_path, "fasta"))
        self.results = {}

    def calculate_symmetry(self, window_size=1000):
        print(" >Calculating GC Skew and Replichore Symmetry...")
        
        for record in self.records:
            seq = str(record.seq).upper()
            seq_len = len(seq)
            
            # GC skew is only biologically meaningful for chromosomes and large plasmids.
            # Right now it skips anything under 100kb.
            if seq_len < 100000:
                continue
                
            gc_skews = []
            
            # Calculate GC Skew for sliding windows
            for i in range(0, seq_len, window_size):
                window = seq[i:i+window_size]
                g = window.count('G')
                c = window.count('C')
                
                if (g + c) == 0:
                    skew = 0
                else:
                    skew = (g - c) / (g + c)
                gc_skews.append(skew)
                
            # Calculate cumulative GC skew
            cumulative_skew = np.cumsum(gc_skews)
            
            # Identify origin (global min) and terminus (global max)
            oric_window = np.argmin(cumulative_skew)
            ter_window = np.argmax(cumulative_skew)
            
            oric_pos = int(oric_window * window_size)
            ter_pos = int(ter_window * window_size)

            plot_theta = [(i / len(cumulative_skew)) * 360 for i in range(len(cumulative_skew))]
            plot_r = cumulative_skew.tolist()
            
            # Calculate the shortest distance along the circular genome
            distance = abs(oric_pos - ter_pos)
            if distance > (seq_len / 2):
                distance = seq_len - distance
                
            # Convert physical distance to degrees (perfect symmetry = 180 degrees)
            symmetry_degrees = (distance / (seq_len / 2)) * 180
            
            # Evaluate biological viability
            # If symmetry is between 165 and 180 degrees, it's highly viable.
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
                "symmetry": round(symmetry_degrees, 1),
                "viability": viability,
                "plot_theta": plot_theta,
                "plot_r": plot_r
            }
            
        return self.results