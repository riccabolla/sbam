import pysam

class JunctionEvaluator:
    def __init__(self, bam_path, orig_lengths):
        self.bam_path = bam_path
        self.orig_lengths = orig_lengths
        self.results = {}

    def evaluate_junctions(self):
        print(" > Evaluating junction spanning scores and read depths")
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        
        # find chr depth and set as baseline for copy number calculations
        raw_stats = {}
        primary_contig = None
        max_len = 0
        
        for contig_id, orig_length in self.orig_lengths.items():
            if orig_length > max_len:
                max_len = orig_length
                primary_contig = contig_id
                
            target_contig = f"{contig_id}_cyclic"
            junction = orig_length
            
            # Dynamic overlap requirement based on contig size
            if orig_length < 20000:
                min_overlap = 500
            else:
                min_overlap = 1000
            
            spanning_reads = 0
            broken_reads = 0
            total_aligned_bases = 0 
            
            try:
                # Calculate Average Depth
                for read in bam.fetch(target_contig):
                    if not read.is_unmapped and not read.is_secondary:
                        total_aligned_bases += read.query_alignment_length
                
                avg_depth = total_aligned_bases / orig_length if orig_length > 0 else 0

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
                        
                raw_stats[contig_id] = {
                    "length": orig_length,
                    "avg_depth": avg_depth,
                    "spanning_reads": spanning_reads,
                    "broken_reads": broken_reads
                }
                
            except ValueError:
                print(f"[WARNING] Contig {target_contig} not found in BAM file.")
                
        bam.close()

        # Establish baseline depth from the chromosome
        baseline_depth = raw_stats.get(primary_contig, {}).get("avg_depth", 1.0)
        if baseline_depth < 1.0:
            baseline_depth = 1.0  # Prevent division by zero
            
        # Classify Plasmids and calculate copy number
        for contig_id, stats in raw_stats.items():
            orig_length = stats["length"]
            avg_depth = stats["avg_depth"]
            spanning_reads = stats["spanning_reads"]
            broken_reads = stats["broken_reads"]
            total_junction_reads = spanning_reads + broken_reads
            
            # Calculate Copy Number relative to the chromosome
            copy_number = avg_depth / baseline_depth
            
            # Basic Status Logic
            if avg_depth < 5.0:
                score = 0.0
                status = "NO_DATA"
            else:
                score = spanning_reads / total_junction_reads if total_junction_reads > 0 else 0
                threshold = 0.4 if orig_length < 50000 else 0.6
                status = "PASS" if score > threshold else "FAIL"
                
            # plasmid classification logic
            classification = "Unknown"
            if orig_length == max_len:
                classification = "Primary Chromosome"
            else:
                if status == "PASS":
                    if copy_number > 10.0:
                        classification = "High-Copy Plasmid (Rolling Circle)"
                    else:
                        classification = "Circular Plasmid"
                else:
                    if copy_number < 0.5:
                        classification = "Low-Coverage Debris"
                    elif spanning_reads == 0 and broken_reads == 0:
                        classification = "Linear Fragment"
                    elif broken_reads > 10:
                        if copy_number > 10.0:
                            classification = "Misassembled High-Copy Plasmid"
                        else:
                            classification = "Misassembled Junction"
                    else:
                        classification = "Incomplete / Linear"
                    
            self.results[contig_id] = {
                "length": orig_length,
                "avg_depth": round(avg_depth, 1),
                "copy_number": round(copy_number, 2),
                "junction_pos": orig_length,
                "spanning_reads": spanning_reads,
                "broken_reads": broken_reads,
                "spanning_score": round(score, 3),
                "status": status,
                "classification": classification
            }
            
        return self.results