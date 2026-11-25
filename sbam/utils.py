"""Shared sequence, math, FASTA loading, and serialization utilities for SBAM."""

from __future__ import annotations

import logging
import math
from dataclasses import asdict, dataclass, is_dataclass
from enum import Enum
from typing import Any

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger("sbam")

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def wilson_lower_bound(k: int, n: int, z: float = 1.96) -> float:
    """Compute the lower bound of the 95% Wilson score CI for a binomial proportion."""
    if n <= 0 or k <= 0:
        return 0.0
    p_hat = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = p_hat + z2 / (2.0 * n)
    rad = z * math.sqrt((p_hat * (1.0 - p_hat) + z2 / (4.0 * n)) / n)
    return max(0.0, (center - rad) / denom)


def compute_empirical_qv(error_count: int, evaluated_bases: int) -> tuple[float, bool]:
    """
    Compute an empirical read-concordance QV proxy: -10 * log10(error_count / evaluated_bases).
    If error_count == 0, returns the lower-bound QV using 0.5 / evaluated_bases and is_capped=True.
    """
    if evaluated_bases <= 0:
        return 0.0, False
    if error_count <= 0:
        return round(-10.0 * math.log10(0.5 / evaluated_bases), 1), True
    return round(-10.0 * math.log10(error_count / evaluated_bases), 1), False


def json_serializable(obj: Any) -> Any:
    """Convert NumPy scalars/arrays, dataclasses, and Enums into native Python JSON types."""
    if isinstance(obj, Enum):
        return obj.value
    if is_dataclass(obj) and not isinstance(obj, type):
        return asdict(obj)
    if hasattr(obj, "tolist"):
        return obj.tolist()
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


@dataclass
class Assembly:
    """In-memory representation of a parsed FASTA assembly shared across SBAM stages."""

    path: str
    records: list[SeqRecord]
    seqs: dict[str, str]
    lengths: dict[str, int]
    primary_id: str | None


def load_assembly(fasta_path: str) -> Assembly:
    """Parse a FASTA file once, filtering out 0-length contigs with a warning."""
    raw_records = list(SeqIO.parse(fasta_path, "fasta"))
    records: list[SeqRecord] = []
    seqs: dict[str, str] = {}
    lengths: dict[str, int] = {}

    for rec in raw_records:
        seq_len = len(rec.seq)
        if seq_len == 0:
            logger.warning("Skipping 0-length contig '%s' in %s", rec.id, fasta_path)
            continue
        records.append(rec)
        seqs[rec.id] = str(rec.seq).upper()
        lengths[rec.id] = seq_len

    primary_id = max(lengths, key=lambda k: lengths[k]) if lengths else None
    return Assembly(
        path=fasta_path,
        records=records,
        seqs=seqs,
        lengths=lengths,
        primary_id=primary_id,
    )
