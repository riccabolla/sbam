"""Typed data models and schema definitions for SBAM."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any

SCHEMA_VERSION = 2


class Verdict(str, Enum):
    ACCEPT = "ACCEPT"
    ACCEPT_WITH_NOTES = "ACCEPT (WITH NOTES)"
    ACCEPT_WITH_CAUTION = "ACCEPT WITH CAUTION"
    ACCEPT_WITH_CAUTION_LOW_DEPTH = "ACCEPT WITH CAUTION (LOW DEPTH <30x)"
    REVIEW = "REVIEW"
    REVIEW_SIZE_ANOMALY = "REVIEW (SIZE ANOMALY)"
    FAIL = "FAIL"
    FAIL_GROSS_BASE_DISCORDANCE = "FAIL (GROSS BASE DISCORDANCE)"
    FAIL_INCOMPLETE = "FAIL (INCOMPLETE)"
    INSUFFICIENT_COVERAGE = "INSUFFICIENT COVERAGE"


class CircularityStatus(str, Enum):
    PASS = "PASS"
    FAIL = "FAIL"
    LINEAR = "LINEAR"
    NO_DATA = "NO_DATA"


class ReplicationViability(str, Enum):
    ACCEPTABLE = "ACCEPTABLE"
    ATYPICAL = "ATYPICAL"
    WARNING = "WARNING"


class _DictMixin:
    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class JunctionMetrics(_DictMixin):
    length: int
    avg_depth: float
    copy_number: float
    junction_pos: int
    spanning_reads: int
    broken_reads: int
    neutral_reads: int
    spanning_score: float
    status: str
    classification: str
    left_terminal_reads: int = 0
    right_terminal_reads: int = 0


@dataclass
class ReplicationMetrics(_DictMixin):
    length: int
    oric_pos: int
    ter_pos: int
    dnaa_boxes_50kb: int
    expected_dnaa_50kb: float
    dnaa_fold_enrichment: float
    dnaa_poisson_pval: float
    oric_offset_bp: int
    oric_offset_kb: float
    oric_offset_deg: float
    rotation_status: str
    symmetry: float
    viability: str
    plot_theta: list[float] = field(default_factory=list)
    plot_r: list[float] = field(default_factory=list)
    plot_x_kb: list[float] = field(default_factory=list)


@dataclass
class MaskedRegion(_DictMixin):
    contig: str
    start: int
    end: int
    length: int
    mean_depth: float
    mean_concordance: float
    reason: str
    category: str = ""


@dataclass
class ConsensusErrorSite(_DictMixin):
    contig: str
    pos_1based: int
    ref_base: str
    depth: int
    concordance: float
    matches: int
    mismatches: int
    insertions: int
    deletions: int
    dominant_error_type: str
    fwd_depth: int
    fwd_concordance: float
    rev_depth: int
    rev_concordance: float
    strand_pvalue: float
    context: str


@dataclass
class ContigBaseMetrics(_DictMixin):
    length: int
    median_depth: float
    masked_bases: int
    masked_pct: float
    evaluated_bases: int
    discordant_bases: int
    consensus_error_bases: int
    subconsensus_bias_bases: int
    consensus_errors_per_mb: float
    empirical_qv: float
    empirical_qv_capped: bool
    low_concordance_blocks: list[dict[str, Any]]
    masked_regions: list[dict[str, Any]]
    total_mismatches: int = 0
    total_insertions: int = 0
    total_deletions: int = 0
    snp_rate: float = 0.0
    indel_rate: float = 0.0
    consensus_error_snp_sites: int = 0
    consensus_error_indel_sites: int = 0
    read_strand_biased_sites: int = 0


@dataclass
class MotifFamily(_DictMixin):
    kmer_size: int
    motif: str
    is_palindrome: bool
    member_motifs: list[str]
    family_size: int
    occurrences: int
    discordant_count: int
    consensus_error_count: int
    raw_consensus_error_count: int
    subconsensus_bias_count: int
    fwd_discordant: int
    rev_discordant: int
    fwd_occurrences: int
    rev_occurrences: int
    strand_biased: bool
    orientation_state: str
    error_rate: float
    fold_increase: float
    ci_lower_fold: float
    max_member_error_rate: float
    max_member_fold: float
    peak_member_motif: str
    p_value: float
    adj_p_value: float
    raw_family_sites: int
    assigned_unique_sites: int
    family_unique_sites: int
    pct_discordant_explained: float
    raw_pct_explained: float
    related_motifs_str: str
    artifact_class: str
    snp_discordant: int = 0
    indel_discordant: int = 0
    read_fwd_depth: int = 0
    read_fwd_match: int = 0
    read_rev_depth: int = 0
    read_rev_match: int = 0
    read_strand_pvalue: float = 1.0
    read_strand_biased: bool = False


@dataclass
class VerdictResult(_DictMixin):
    verdict: str
    color: str
    message: str
    info_notes: list[str]
    review_triggers: list[str]
    completeness_status: str
    completeness_msg: str
    completeness_pct: float | None
    contig_summaries: list[dict[str, Any]]
    primary_contig_id: str | None
    primary_length: int
    primary_depth: float
    total_length: int
    failed_contigs_count: int
    low_concord_blocks: list[dict[str, Any]]
    low_concord_block_count: int
    low_concord_total_bp: int
    evaluated_bases: int
    discordant_bases: int
    consensus_error_bases: int
    subconsensus_bias_bases: int
    baseline_discordance_rate: float
    consensus_error_rate: float
    subconsensus_bias_rate: float
    consensus_errors_per_mb: float
    subconsensus_bias_per_mb: float
    discordant_per_mb: float
    qv_display: str
    explained_discordant_bases: int
    explained_discordant_pct: float
    depth_sufficient_for_base_qc: bool
    primary_junc_score: str
    j_status: str
    primary_symmetry: str
    p_status: str
    primary_rotation: str
