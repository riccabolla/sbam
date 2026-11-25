# sbam: Structural and Base-Level Analysis of Microbial Assemblies

`sbam` is a reference-free quality control and validation toolkit designed for long-read (Oxford Nanopore and PacBio) microbial genome assemblies.

Modern long-read assemblers (such as Flye, Raven, Canu, and consensus pipelines like Autocycler or Trycycler) routinely produce single-contig bacterial chromosomes and plasmids. However, standard assembly metrics ($N_{50}$, contig count, or database-driven gene completeness tools) cannot verify whether a circular contig is physically supported across its circularization seam, whether large-scale inversions have scrambled its replichores, or whether basecaller and polisher models have left consensus errors or systematic motif discordance.

`sbam` evaluates assemblies directly against the physical sequencing reads and intrinsic sequence geometry across four orthogonal dimensions:

1. **Physical Circularity & Replicon Copy Number (`sbam.struct.eval`)**: Appends a cyclic sequence buffer across contig boundaries to test continuous read spanning while filtering out natural unclipped read ends, distinguishing intact circular chromosomes and plasmids (including sub-stoichiometric and high-copy replicons) from misassembled junctions and linear fragments. Also supports naturally linear replicons via `--linear`.
2. **Replication Structure, DnaA Enrichment & Origin Rotation (`sbam.struct.ori`)**: Computes zero-centered, circularly smoothed cumulative GC skew across chromosomal replicons ($\ge 500\text{ kb}$) to verify bidirectional theta-replication symmetry ($\approx 180^\circ$), tests for 9-bp DnaA-box motif enrichment around *oriC* against the genome-wide Poisson background, and classifies *oriC*-to-$0\text{ bp}$ rotation status (`ROTATED`, `NEAR_ORIC`, `UNROTATED`).
3. **All-Contig Base-Level Concordance & Motif Family Clustering (`sbam.base.base`)**: Performs multi-threaded, indel-aware read pileups (with cyclic buffer depth folding) across **all contigs** to quantify majority consensus errors ($< 50\%$ read support, reported as an empirical read-concordance $\text{QV}$ and errors per $\text{Mb}$), masks structural anomalies ($\ge 20\text{ bp}$), and clusters statistically enriched $k$-mers into deduplicated **Motif Families** with 3-state reference-orientation classification.
4. **Genome Size Completeness (`sbam.report.report`)**: Evaluates total assembly length and primary chromosome proportion against expected species genome size when `-g / --genome_size` is provided.

---

## Multi-Tier Rate-Based Verdict System

Rather than relying on brittle zero-tolerance raw counts that scale poorly with genome size, `sbam` integrates all structural, geometric, and base-level signals into a calibrated rate-based verdict:

| Verdict | Criteria |
| :--- | :--- |
| **`ACCEPT`** | Primary contig passes circularity (or `--linear`), all secondary contigs pass circularity, replication symmetry is `ACCEPTABLE` ($165^\circ\text{--}180^\circ$), primary depth $\ge 30\times$, $0$ majority consensus errors ($< 50\%$ read support), $0$ low-concordance ($< 90\%$) excluded blocks, and structurally masked bases $\le 5.0\%$. |
| **`ACCEPT (WITH NOTES)`** | Structure and replication pass with depth $\ge 30\times$, majority consensus error rate is $\le 5.0\text{ / Mb}$ (empirical read-concordance $\text{QV} \ge 53.0$), and low-concordance ($< 90\%$) excluded blocks total $\le 1,000\text{ bp}$ ($\le 1$ block), but minor non-zero inspection notes exist (highlighted directly in the summary banner). |
| **`ACCEPT WITH CAUTION`** | Structural continuity passes, but the primary contig is $< 500\text{ kb}$ (no chromosomal replication check) **or** primary depth is between $5\times$ and $30\times$ (`LOW DEPTH <30x`, insufficient coverage for high-confidence base-level grading). |
| **`REVIEW`** | Triggered if any secondary contig fails circularity, structurally masked bases exceed $5.0\%$, majority consensus error rate exceeds $5.0\text{ / Mb}$ ($\text{QV} < 53.0$), low-concordance ($< 90\%$) excluded blocks exceed $1,000\text{ bp}$ (or $\ge 2$ blocks), a plasmid exceeds $5.0\text{ / Mb}$ consensus errors ($\ge 2$ sites), replication symmetry is `ATYPICAL` ($150^\circ\text{--}165^\circ$) or `WARNING` ($< 150^\circ$), or size completeness is anomalous (`REVIEW (SIZE ANOMALY)`). |
| **`FAIL`** | Primary contig fails circularity, genome completeness is $< 90\%$ (`FAIL (INCOMPLETE)`), or primary majority consensus error rate exceeds $50.0\text{ / Mb}$ (`FAIL (GROSS BASE DISCORDANCE)`, $\text{QV} < 43.0$). |
| **`INSUFFICIENT COVERAGE`** | Average read depth on the primary contig is $< 5.0\times$ (`NO_DATA`). |

---

## Pipeline Architecture & Module Reference

| Module | Class / Function | Primary Responsibility |
| :--- | :--- | :--- |
| `sbam.cli` | `main()`, `run_summarize()` | Parses CLI arguments, manages the unified logger (`--log`), orchestrates validation stages, exports JSON/HTML outputs, cleans intermediate files (`--clean`), and aggregates multi-sample cohorts and recurrent batch motifs (`sbam summarize`). |
| `sbam.struct.struct` | `GenomeStructure` | Reads the input FASTA and constructs the extended cyclic reference (`cyclic_assembly.fasta`) by appending the first $B\text{ bp}$ (`--buffer-size`) of each contig to its $3'$ end. |
| `sbam.struct.mapping` | `ReadAligner` | Streams `minimap2` alignments directly into `samtools sort` and `samtools index` to produce `alignment.sorted.bam` without intermediate SAM disk I/O. |
| `sbam.struct.eval` | `JunctionEvaluator` | Quantifies contig read depth, relative copy number, and junction spanning scores at coordinate $L$ while separating true clipped/broken alignments from natural unclipped read terminations. |
| `sbam.struct.ori` | `BioVal` | Calculates zero-centered, circularly smoothed GC skew, locates *oriC* and *ter*, computes replichore symmetry ($^\circ$), tests DnaA-box Poisson enrichment, and classifies *oriC*-to-$0\text{ bp}$ rotation status. |
| `sbam.base.base` | `BaseAccuracy` | Profiles all contigs via parallel `pysam` pileups, folds cyclic buffer coverage back into $[0, B)$, separates majority ($< 50\%$) consensus errors from sub-consensus ($50\%\text{--}85\%$) read bias, and clusters discordant $k$-mers into mutually exclusive Motif Families. |
| `sbam.report.report` | `DashboardBuilder` | Compiles all structural, biological, and base-level metrics into a standalone HTML dashboard (`*_sbam_report.html`) with Plotly polar plots and an automatic zero-dependency inline SVG fallback (`--offline`). |

---

## Quick Start

### Single Isolate Evaluation

```bash
sbam -a assembly.fasta -r reads.fastq.gz -o sbam_out -t 8 -g 5.3m --clean --log
```

### Multi-Sample Cohort & Recurrent Motif Summary

After running `sbam` across multiple isolates in a batch or workflow (e.g., Nextflow or Snakemake), aggregate all sample JSON outputs into a cohort summary table and a cross-sample recurrent motif report:

```bash
sbam summarize results/*/*_sbam_metrics.json -o cohort_summary.tsv
```