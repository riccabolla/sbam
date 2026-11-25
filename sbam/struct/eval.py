"""Structural continuity and circular junction / linear terminus evaluator."""

from __future__ import annotations

import logging
from typing import Any

import pysam

logger = logging.getLogger("sbam")


class JunctionEvaluator:
    """Evaluates read spanning across circular seams or terminal integrity on linear replicons."""

    def __init__(
        self,
        bam_path: str,
        orig_lengths: dict[str, int],
        min_clip_len: int = 50,
        boundary_tol: int = 15,
        clip_window: int = 500,
        chr_threshold: float = 0.60,
        plasmid_threshold: float = 0.40,
        linear_mode: bool = False,
        is_cyclic: bool | None = None,
    ):
        self.bam_path = bam_path
        self.orig_lengths = orig_lengths
        self.min_clip_len = min_clip_len
        self.boundary_tol = boundary_tol
        self.clip_window = clip_window
        self.chr_threshold = chr_threshold
        self.plasmid_threshold = plasmid_threshold
        self.linear_mode = linear_mode
        self.is_cyclic = is_cyclic
        self.results: dict[str, dict[str, Any]] = {}

    @staticmethod
    def _leading_clip_len(read: Any) -> int:
        if not read.cigartuples:
            return 0
        op, length = read.cigartuples[0]
        return length if op in (4, 5) else 0

    @staticmethod
    def _trailing_clip_len(read: Any) -> int:
        if not read.cigartuples:
            return 0
        op, length = read.cigartuples[-1]
        return length if op in (4, 5) else 0

    def _resolve_target_contig(self, bam: Any, contig_id: str) -> tuple[str, bool]:
        """Resolve the BAM reference name and whether it uses the `_cyclic` suffix."""
        if self.is_cyclic is True:
            return f"{contig_id}_cyclic", True
        if self.is_cyclic is False:
            return contig_id, False

        refs = getattr(bam, "references", None)
        if isinstance(refs, (list, tuple, set)):
            cyclic_name = f"{contig_id}_cyclic"
            if cyclic_name in refs:
                return cyclic_name, True
            if contig_id in refs:
                return contig_id, False

        if self.linear_mode:
            return contig_id, False
        return f"{contig_id}_cyclic", True

    def evaluate_junctions(self) -> dict[str, dict[str, Any]]:
        """Single-pass evaluation of depth, circular junction spanning, or linear termini."""
        logger.info(" > Evaluating junction spanning scores and read depths")

        raw_stats: dict[str, dict[str, Any]] = {}
        primary_contig = None
        max_len = 0

        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            for contig_id, orig_length in self.orig_lengths.items():
                if orig_length > max_len:
                    max_len = orig_length
                    primary_contig = contig_id

                target_contig, contig_is_cyclic = self._resolve_target_contig(bam, contig_id)
                junction = orig_length
                min_overlap = 1000 if orig_length >= 20000 else min(500, max(1, orig_length // 2))

                spanning_reads = 0
                broken_reads = 0
                neutral_reads = 0
                left_terminal_reads = 0
                right_terminal_reads = 0
                total_aligned_bases = 0

                try:
                    for read in bam.fetch(target_contig):
                        if (
                            read.is_unmapped
                            or read.is_secondary
                            or read.is_supplementary
                            or getattr(read, "is_qcfail", False)
                            or getattr(read, "is_duplicate", False)
                        ):
                            continue

                        aligned_len = (
                            read.reference_length
                            if read.reference_length is not None
                            else read.query_alignment_length
                        )
                        if aligned_len:
                            total_aligned_bases += aligned_len

                        r_start = read.reference_start
                        r_end = read.reference_end
                        if r_start is None or r_end is None:
                            continue

                        if self.linear_mode and not contig_is_cyclic:
                            lead_clip = self._leading_clip_len(read)
                            trail_clip = self._trailing_clip_len(read)

                            if r_start <= self.clip_window:
                                left_terminal_reads += 1
                                if lead_clip >= self.min_clip_len:
                                    broken_reads += 1
                                elif r_start <= self.boundary_tol:
                                    spanning_reads += 1
                                else:
                                    neutral_reads += 1

                            if r_end >= orig_length - self.clip_window:
                                right_terminal_reads += 1
                                if trail_clip >= self.min_clip_len:
                                    broken_reads += 1
                                elif r_end >= orig_length - self.boundary_tol:
                                    spanning_reads += 1
                                else:
                                    neutral_reads += 1
                            continue

                        # Circular junction evaluation across [L - clip_window, L + clip_window]
                        if r_end < junction - self.clip_window or r_start > junction + self.clip_window:
                            continue

                        left_extension = junction - r_start
                        right_extension = r_end - junction

                        if left_extension >= min_overlap and right_extension >= min_overlap:
                            spanning_reads += 1
                            continue

                        lead_clip = self._leading_clip_len(read)
                        trail_clip = self._trailing_clip_len(read)

                        clipped_near_left = (
                            abs(r_start - junction) <= self.clip_window
                            and left_extension < min_overlap
                            and lead_clip >= self.min_clip_len
                        )
                        clipped_near_right = (
                            abs(r_end - junction) <= self.clip_window
                            and right_extension < min_overlap
                            and trail_clip >= self.min_clip_len
                        )
                        stops_at_junction = (
                            abs(r_end - junction) <= self.boundary_tol
                            or abs(r_start - junction) <= self.boundary_tol
                        )

                        if clipped_near_left or clipped_near_right or stops_at_junction:
                            broken_reads += 1
                        elif r_end >= junction - self.boundary_tol and r_start <= junction + self.boundary_tol:
                            neutral_reads += 1

                    avg_depth = total_aligned_bases / orig_length if orig_length > 0 else 0.0
                    raw_stats[contig_id] = {
                        "length": orig_length,
                        "avg_depth": avg_depth,
                        "spanning_reads": spanning_reads,
                        "broken_reads": broken_reads,
                        "neutral_reads": neutral_reads,
                        "left_terminal_reads": left_terminal_reads,
                        "right_terminal_reads": right_terminal_reads,
                    }

                except ValueError:
                    logger.warning("Contig %s not found in BAM file.", target_contig)

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
            elif self.linear_mode:
                if total_junction_reads > 0:
                    score = spanning_reads / total_junction_reads
                else:
                    score = 1.0
                status = "LINEAR"
            elif spanning_reads == 0:
                score = 0.0
                status = "FAIL"
            else:
                score = spanning_reads / total_junction_reads if total_junction_reads > 0 else 0.0
                threshold = self.plasmid_threshold if orig_length < 50000 else self.chr_threshold
                status = "PASS" if score > threshold else "FAIL"

            if orig_length == max_len:
                classification = (
                    "Linear Primary Chromosome" if self.linear_mode else "Primary Chromosome"
                )
            elif orig_length >= 500000:
                classification = (
                    "Secondary Chromosome / Megaplasmid"
                    if status in ("PASS", "LINEAR")
                    else "Misassembled Secondary Replicon"
                )
            else:
                if status == "PASS":
                    if copy_number >= 10.0:
                        classification = "High-Copy Circular Plasmid"
                    elif copy_number < 0.7:
                        classification = "Sub-stoichiometric Circular Plasmid"
                    else:
                        classification = "Circular Plasmid"
                elif status == "LINEAR":
                    classification = "Linear Plasmid / Replicon"
                else:
                    if copy_number < 0.5:
                        classification = "Low-Coverage Debris"
                    elif spanning_reads == 0 and broken_reads == 0:
                        classification = "Linear Fragment"
                    elif broken_reads > 10:
                        if copy_number >= 10.0:
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
                "left_terminal_reads": stats["left_terminal_reads"],
                "right_terminal_reads": stats["right_terminal_reads"],
                "spanning_score": round(score, 3),
                "status": status,
                "classification": classification,
            }

        return self.results
