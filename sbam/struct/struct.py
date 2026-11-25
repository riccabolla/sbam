"""Cyclic buffer preparation and reference handling for SBAM."""

from __future__ import annotations

import logging
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger("sbam")


class GenomeStructure:
    """Prepares cyclic or linear assembly references for read alignment."""

    def __init__(self, fasta_path: str, outdir: str, records: list[SeqRecord] | None = None):
        self.fasta_path = fasta_path
        self.outdir = outdir
        raw_records = records if records is not None else list(SeqIO.parse(self.fasta_path, "fasta"))
        self.records: list[SeqRecord] = []
        for rec in raw_records:
            if len(rec.seq) == 0:
                logger.warning("Skipping 0-length contig '%s' in %s", rec.id, self.fasta_path)
                continue
            self.records.append(rec)

    def create_cyclic_buffer(self, buffer_size: int = 50000) -> tuple[str, dict[str, int]]:
        """
        Create a FASTA where the first `buffer_size` bases of each contig
        are appended to its 3' end to simulate circular continuity.
        """
        cyclic_records: list[SeqRecord] = []
        orig_lengths: dict[str, int] = {}

        for record in self.records:
            seq_len = len(record.seq)
            actual_buffer = min(buffer_size, seq_len)
            if actual_buffer == 0:
                continue

            orig_lengths[record.id] = seq_len
            buffer_seq = record.seq[:actual_buffer]
            cyclic_seq = record.seq + buffer_seq

            cyclic_record = SeqRecord(
                cyclic_seq,
                id=f"{record.id}_cyclic",
                description=f"original_length={seq_len} buffer={actual_buffer}",
            )
            cyclic_records.append(cyclic_record)

        out_fasta = os.path.join(self.outdir, "cyclic_assembly.fasta")
        SeqIO.write(cyclic_records, out_fasta, "fasta")
        return out_fasta, orig_lengths

    def prepare_reference(
        self, buffer_size: int = 50000, linear: bool = False
    ) -> tuple[str, dict[str, int], bool]:
        """
        Prepare the alignment reference FASTA.

        When `linear=True`, bypasses cyclic buffer creation and returns the original FASTA.
        Returns `(reference_fasta_path, orig_lengths, is_cyclic)`.
        """
        orig_lengths = {rec.id: len(rec.seq) for rec in self.records if len(rec.seq) > 0}
        if linear:
            return self.fasta_path, orig_lengths, False
        cyclic_fasta, orig_lengths = self.create_cyclic_buffer(buffer_size=buffer_size)
        return cyclic_fasta, orig_lengths, True
