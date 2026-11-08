import pysam


class JunctionEvaluator:
    def __init__(self, bam_path, orig_lengths):
        self.bam_path = bam_path
        self.orig_lengths = orig_lengths
        self.results = {}

    def evaluate_junctions(self):
        print(" > Evaluating junction spanning scores and read depths")

        raw_stats = {}
        primary_contig = None
        max_len = 0

        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            for contig_id, orig_length in self.orig_lengths.items():
                if orig_length > max_len:
                    max_len = orig_length
                    primary_contig = contig_id

                target_contig = f"{contig_id}_cyclic"
                junction = orig_length

                # Dynamic overlap requirement based on contig size
                min_overlap = 500 if orig_length < 20000 else 1000

                spanning_reads = 0
                broken_reads = 0
                total_aligned_bases = 0

                try:
                    # Calculate average depth across the original contig length.
                    # Because primary reads from the [0, buffer) region split between
                    # [0, buffer) and [orig_length, orig_length + buffer), summing
                    # primary alignments across target_contig and dividing by orig_length
                    # recovers the true unfolded depth.
                    for read in bam.fetch(target_contig):
                        if (
                            read.is_unmapped
                            or read.is_secondary
                            or read.is_supplementary
                        ):
                            continue
                        aligned_len = (
                            read.reference_length
                            if read.reference_length is not None
                            else read.query_alignment_length
                        )
                        if aligned_len:
                            total_aligned_bases += aligned_len

                    avg_depth = (
                        total_aligned_bases / orig_length if orig_length > 0 else 0.0
                    )

                    # Evaluate junction spanning at coordinate `orig_length`
                    for read in bam.fetch(target_contig, junction - 1, junction + 1):
                        if (
                            read.is_unmapped
                            or read.is_secondary
                            or read.is_supplementary
                        ):
                            continue
                        if read.reference_start is None or read.reference_end is None:
                            continue

                        left_extension = junction - read.reference_start
                        right_extension = read.reference_end - junction

                        if (
                            left_extension >= min_overlap
                            and right_extension >= min_overlap
                        ):
                            spanning_reads += 1
                        else:
                            broken_reads += 1

                    raw_stats[contig_id] = {
                        "length": orig_length,
                        "avg_depth": avg_depth,
                        "spanning_reads": spanning_reads,
                        "broken_reads": broken_reads,
                    }

                except ValueError:
                    print(f"[WARNING] Contig {target_contig} not found in BAM file.")

        # Establish baseline depth from the primary chromosome
        baseline_depth = raw_stats.get(primary_contig, {}).get("avg_depth", 1.0)
        if baseline_depth < 1.0:
            baseline_depth = 1.0  # Prevent division by zero

        # Classify contigs/plasmids and calculate copy number
        for contig_id, stats in raw_stats.items():
            orig_length = stats["length"]
            avg_depth = stats["avg_depth"]
            spanning_reads = stats["spanning_reads"]
            broken_reads = stats["broken_reads"]
            total_junction_reads = spanning_reads + broken_reads

            copy_number = avg_depth / baseline_depth

            if avg_depth < 5.0:
                score = 0.0
                status = "NO_DATA"
            else:
                score = (
                    spanning_reads / total_junction_reads
                    if total_junction_reads > 0
                    else 0.0
                )
                threshold = 0.4 if orig_length < 50000 else 0.6
                status = "PASS" if score > threshold else "FAIL"

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
                "classification": classification,
            }

        return self.results