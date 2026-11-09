import pysam


class JunctionEvaluator:
    def __init__(self, bam_path, orig_lengths, min_clip_len=50, boundary_tol=15):
        self.bam_path = bam_path
        self.orig_lengths = orig_lengths
        self.min_clip_len = min_clip_len
        self.boundary_tol = boundary_tol
        self.results = {}

    @staticmethod
    def _leading_clip_len(read):
        if not read.cigartuples:
            return 0
        op, length = read.cigartuples[0]
        return length if op in (4, 5) else 0

    @staticmethod
    def _trailing_clip_len(read):
        if not read.cigartuples:
            return 0
        op, length = read.cigartuples[-1]
        return length if op in (4, 5) else 0

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
                neutral_reads = 0
                total_aligned_bases = 0

                try:
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

                    # Fetch reads within a small window around the artificial junction
                    fetch_start = max(0, junction - self.boundary_tol)
                    fetch_end = junction + self.boundary_tol

                    for read in bam.fetch(target_contig, fetch_start, fetch_end):
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

                        # 1. Read solidly bridges the junction on both sides
                        if (
                            left_extension >= min_overlap
                            and right_extension >= min_overlap
                        ):
                            spanning_reads += 1
                            continue

                        # 2. Check if the read is genuinely broken/clipped at the junction
                        lead_clip = self._leading_clip_len(read)
                        trail_clip = self._trailing_clip_len(read)

                        clipped_near_left = (
                            left_extension < min_overlap
                            and lead_clip >= self.min_clip_len
                        )
                        clipped_near_right = (
                            right_extension < min_overlap
                            and trail_clip >= self.min_clip_len
                        )
                        stops_at_junction = (
                            abs(read.reference_end - junction) <= self.boundary_tol
                            or abs(read.reference_start - junction) <= self.boundary_tol
                        )

                        if clipped_near_left or clipped_near_right or stops_at_junction:
                            broken_reads += 1
                        else:
                            # Natural unclipped read end landing within min_overlap of junction
                            neutral_reads += 1

                    raw_stats[contig_id] = {
                        "length": orig_length,
                        "avg_depth": avg_depth,
                        "spanning_reads": spanning_reads,
                        "broken_reads": broken_reads,
                        "neutral_reads": neutral_reads,
                    }

                except ValueError:
                    print(f"[WARNING] Contig {target_contig} not found in BAM file.")

        baseline_depth = raw_stats.get(primary_contig, {}).get("avg_depth", 1.0)
        if baseline_depth < 1.0:
            baseline_depth = 1.0

        for contig_id, stats in raw_stats.items():
            orig_length = stats["length"]
            avg_depth = stats["avg_depth"]
            spanning_reads = stats["spanning_reads"]
            broken_reads = stats["broken_reads"]
            neutral_reads = stats["neutral_reads"]
            total_junction_reads = spanning_reads + broken_reads

            copy_number = avg_depth / baseline_depth

            if avg_depth < 5.0:
                score = 0.0
                status = "NO_DATA"
            elif spanning_reads == 0:
                score = 0.0
                status = "FAIL"
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
                "neutral_reads": neutral_reads,
                "spanning_score": round(score, 3),
                "status": status,
                "classification": classification,
            }

        return self.results