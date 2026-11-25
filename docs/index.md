# sbam: Structural and Base-Level Analysis of Microbial Assemblies

`sbam` is a reference-free quality control and validation toolkit designed for long-read (Oxford Nanopore and PacBio) microbial genome assemblies.

Modern long-read assemblers (such as Flye, Raven, Canu, and consensus pipelines like Autocycler or Trycycler) routinely produce single-contig bacterial chromosomes and plasmids. However, standard assembly metrics ($N_{50}$, contig count, or database-driven gene completeness tools) cannot verify whether a circular contig is physically supported across its circularization seam, whether large-scale inversions have scrambled its replichores, or whether basecaller and polisher models have left consensus errors or systematic motif discordance.

`sbam` evaluates assemblies directly against the physical sequencing reads (either by aligning raw reads or ingesting a pre-aligned `--bam`) and intrinsic sequence geometry across four orthogonal dimensions:

1. **Physical Circularity, Linear Termini & Copy Number (`sbam.struct.eval`)**: In standard circular mode, appends a cyclic sequence buffer across contig boundaries to test continuous read spanning across coordinate $L$ while filtering out natural unclipped read ends—distinguishing intact circular chromosomes and plasmids (including sub-stoichiometric and high-copy replicons) from misassembled junctions and linear fragments. With `--linear`, bypasses cyclic buffer construction and evaluates $5'$ ($0$) and $3'$ ($L$) termini directly.
2. **Replication Structure, DnaA Enrichment & Origin Rotation (`sbam.struct.ori`)**: Computes length-weighted, zero-centered, circularly (or edge-padded linear) smoothed cumulative GC skew across chromosomal replicons ($\ge 500\text{ kb}$ by default, configurable via `--min-ori-length`) to verify bidirectional theta-replication symmetry ($\approx 180^\circ$), scans with a vectorized 2-bit lookup table for 9-bp DnaA-box enrichment around *oriC* against the genome-wide Poisson background (wrapping across the $0/L$ seam in circular mode), and classifies *oriC*-to-$0\text{ bp}$ rotation status (`ROTATED`, `NEAR_ORIC`, `UNROTATED`, or `LINEAR`).
3. **All-Contig Base-Level Concordance, Error Stratification & Motif Family Clustering (`sbam.base.base`)**: Performs multi-threaded, read-centric `--eqx` CIGAR sweeps (with cyclic buffer depth folding in circular mode and vectorized `M`-op fallback for external BAMs) across **all contigs** to quantify majority consensus errors ($< 50\%$ read support, reported as empirical read-concordance $\text{QV}$ and errors per $\text{Mb}$), stratify substitutions vs. indels and sequencing read-strand bias (two-sided Fisher's exact test), mask structural anomalies ($\ge 20\text{ bp}$), and cluster statistically enriched $k$-mers (including dual central-base testing on even-length palindromes) into mutually exclusive **Motif Families**.
4. **Verdict Engine, Completeness & Diagnostic Tracks (`sbam.verdict`, `sbam.report`)**: Evaluates total assembly length and primary chromosome proportion against expected species genome size (`-g / --genome_size`), integrates all structural and base-level signals into a calibrated rate-based verdict, generates an interactive HTML dashboard (`<prefix>_sbam_report.html`), versioned schema v2 JSON (`<prefix>_sbam_metrics.json`), and actionable IGV/polishing tracks (`<prefix>_masked_regions.bed`, `<prefix>_consensus_errors.tsv`, `<prefix>_consensus_errors.vcf`).

---

## Multi-Tier Rate-Based Verdict System

Rather than relying on brittle zero-tolerance raw counts that scale poorly with genome size, `sbam` integrates all structural, geometric, and base-level signals into a calibrated rate-based verdict (`sbam.verdict.evaluate_verdict`):

| Verdict | Criteria |
| :--- | :--- |
| **`ACCEPT`** | Primary contig passes circularity (or `--linear`), all secondary contigs pass circularity, replication symmetry is `ACCEPTABLE` ($165^\circ\text{--}180^\circ$), primary depth $\ge 30\times$, $0$ majority consensus errors ($< 50\%$ read support), $0$ low-concordance ($< 90\%$) excluded blocks, and structurally masked bases $\le 5.0\%$. |
| **`ACCEPT (WITH NOTES)`** | Structure and replication pass with depth $\ge 30\times$, majority consensus error rate is $\le 5.0\text{ / Mb}$ (empirical read-concordance $\text{QV} \ge 53.0$), and low-concordance ($< 90\%$) excluded blocks total $\le 1,000\text{ bp}$ ($\le 1$ block), with minor non-zero inspection notes highlighted in the summary banner. |
| **`ACCEPT WITH CAUTION`** | Structural continuity passes without review triggers, but the primary contig is smaller than `--min-ori-length` (default $500\text{ kb}$, skipping chromosomal replication check) **or** primary depth is between $5\times$ and $30\times$ (`ACCEPT WITH CAUTION (LOW DEPTH <30x)`). Low-depth assemblies still escalate to `REVIEW` or `REVIEW (SIZE ANOMALY)` if secondary contigs fail, masking exceeds $5.0\%$, or genome size is anomalous. |
| **`REVIEW`** | Triggered if any secondary contig fails circularity, primary or secondary contig structurally masked bases exceed $5.0\%$, primary or secondary majority consensus error rate exceeds $5.0\text{ / Mb}$ (on secondary contigs with $\ge 2$ sites), primary or secondary low-concordance ($< 90\%$) excluded blocks exceed $1,000\text{ bp}$ (or $\ge 2$ blocks), replication symmetry is `ATYPICAL` ($150^\circ\text{--}165^\circ$) or `WARNING` ($< 150^\circ$), or size completeness is anomalous (`REVIEW (SIZE ANOMALY)`). |
| **`FAIL`** | Primary contig fails circularity, genome completeness is $< 90\%$ (`FAIL (INCOMPLETE)`), or primary majority consensus error rate exceeds $50.0\text{ / Mb}$ (`FAIL (GROSS BASE DISCORDANCE)`, $\text{QV} < 43.0$). |
| **`INSUFFICIENT COVERAGE`** | Average read depth on the primary contig is $< 5.0\times$ (`NO_DATA`). |

---

## Pipeline Architecture & Module Reference

| Module | Class / Function | Primary Responsibility |
| :--- | :--- | :--- |
| `sbam.cli` | `main()`, `run_summarize()` | Parses CLI arguments (`-r`/`--bam`, `--prefix`, `--linear`, `--emit-tracks`, `-V/--version`), configures Python `logging` (`--log`), loads the FASTA once, orchestrates validation stages, exports JSON/HTML/BED/TSV/VCF outputs, cleans SBAM-created intermediates (`--clean`), and aggregates multi-sample cohorts (`sbam summarize`). |
| `sbam.utils` | `load_assembly()`, `revcomp()`, `compute_empirical_qv()` | Provides single-pass FASTA parsing (`Assembly`, filtering $0\text{ bp}$ contigs), sequence utilities, Wilson score lower bounds, empirical $\text{QV}$ calculations, and JSON serialization helpers. |
| `sbam.models` | `SCHEMA_VERSION`, `Verdict`, Dataclasses | Defines versioned JSON schema (`SCHEMA_VERSION = 2`), verdict/status Enums, and typed dataclasses for structural, replication, and base-fidelity metrics. |
| `sbam.struct.struct` | `GenomeStructure` | Filters empty contigs and constructs `cyclic_assembly.fasta` by appending $B\text{ bp}$ (`--buffer-size`) to each contig's $3'$ end, or bypasses buffer creation when `--linear` is set (`prepare_reference()`). |
| `sbam.struct.mapping` | `ReadAligner`, `AlignmentError` | Streams `minimap2 -a --eqx --secondary=no` directly into `samtools sort --write-index` to produce `alignment.sorted.bam` and `.bai` in a single disk pass, with isolated stderr logging and clean process cleanup. |
| `sbam.struct.eval` | `JunctionEvaluator` | Performs single-pass BAM evaluation of contig depth, relative copy number, and junction spanning scores at coordinate $L$ (or $5'/3'$ termini in `--linear` mode) while distinguishing true junction clips from natural unclipped read ends. |
| `sbam.struct.ori` | `BioVal` | Calculates length-weighted, zero-centered cumulative GC skew (circular or edge-padded linear), locates *oriC* and *ter*, computes replichore symmetry ($^\circ$), runs vectorized DnaA-box Poisson enrichment, and classifies rotation status. |
| `sbam.base.base` | `BaseAccuracy` | Profiles all contigs via flat-chunked CIGAR sweeps (`_profile_chunk`), excludes supplementary/QC-fail/duplicate alignments, folds cyclic buffer coverage, stratifies SNPs vs. indels and read-strand bias, masks $\ge 20\text{ bp}$ structural anomalies, and clusters enriched $k$-mers into mutually exclusive Motif Families. |
| `sbam.verdict` | `evaluate_verdict()` | Pure verdict engine computing assembly status, escalation triggers, completeness checks, and inspection notes independently of HTML rendering. |
| `sbam.report.tracks` | `export_masked_regions_bed()`, `export_consensus_errors_tsv()`, `export_consensus_errors_vcf()` | Exports machine-readable BED6, TSV, and VCFv4.2 diagnostic tracks for IGV inspection and automated polishing. |
| `sbam.report.report` | `DashboardBuilder` | Renders `sbam/report/templates/report.html.j2` via Jinja2 `PackageLoader` into `<prefix>_sbam_report.html` with interactive Plotly polar/linear maps and automatic inline SVG fallback (`--offline`). |

---

## Quick Start

### Single Isolate Evaluation

```bash
sbam -a assembly.fasta -r reads.fastq.gz -o sbam_out --prefix isolate_01 -t 8 -g 5.3m --clean --log
```

### Re-Running QC Thresholds on a Pre-Aligned BAM

```bash
sbam -a assembly.fasta --bam alignment.sorted.bam -o sbam_rerun --prefix isolate_01_strict --min-concordance 0.90 --offline
```

### Multi-Sample Cohort & Recurrent Motif Summary

After running `sbam` across multiple isolates in a batch or workflow (e.g., Nextflow or Snakemake), aggregate all sample JSON outputs into a cohort summary table and a cross-sample recurrent motif report:

```bash
sbam summarize "results/*/*_sbam_metrics.json" -o cohort_summary.tsv
```
