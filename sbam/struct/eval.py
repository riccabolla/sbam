import pysam
import os

class JunctionEvaluator:
    def __init__(self, bam_path, orig_lengths):
        self.bam_path = bam_path
        self.orig_lengths = orig_lengths
        self.results = {}

    def evaluate_junctions(self):
        print(" > Evaluating junction spanning scores and read depths")
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        
        for contig_id, orig_length in self.orig_lengths.items():
            target_contig = f"{contig_id}_cyclic"
            junction = orig_length
            
            # Dynamic overlap requirement based on contig size
            if orig_length < 20000:
                min_overlap = 500
            else:
                min_overlap = 1000
            
            spanning_reads = 0
            broken_reads = 0
            total_aligned_bases = 0  # New variable for depth
            
            try:
                # Calculate Average Depth
                for read in bam.fetch(target_contig):
                    if not read.is_unmapped and not read.is_secondary:
                        total_aligned_bases += read.query_alignment_length
                
                avg_depth = total_aligned_bases / orig_length

                # Evaluate Junction Spanning 
                for read in bam.fetch(target_contig, junction - 1, junction + 1):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    
                    left_extension = junction - read.reference_start
                    right_extension = read.reference_end - junction
                    
                    if left_extension >= min_overlap and right_extension >= min_overlap:
                        spanning_reads += 1
                    else:
                        broken_reads += 1
                        
                total_junction_reads = spanning_reads + broken_reads
                
                # handling low/zero coverage
                if avg_depth < 5.0:
                    score = 0.0
                    status = "NO_DATA"
                else:
                    score = spanning_reads / total_junction_reads if total_junction_reads > 0 else 0
                    
                    # Dynamic threshold
                    if orig_length < 50000:
                        status = "PASS" if score > 0.4 else "FAIL"
                    else:
                        status = "PASS" if score > 0.6 else "FAIL"
                        
                self.results[contig_id] = {
                    "length": orig_length,
                    "avg_depth": round(avg_depth, 1),
                    "junction_pos": junction,
                    "spanning_reads": spanning_reads,
                    "broken_reads": broken_reads,
                    "spanning_score": round(score, 3),
                    "status": status
                }
                
            except ValueError:
                print(f"[WARNING] Contig {target_contig} not found in BAM file.")
                
        bam.close()
        return self.results