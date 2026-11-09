# sbam: Structural and Base-Level Analysis of Microbial Assemblies

`sbam` is a reference-free quality control and validation toolkit designed specifically for long-read (Oxford Nanopore and PacBio) microbial genome assemblies. 

Modern long-read assemblers (such as Flye, Raven, Canu, and consensus pipelines like Autocycler or Trycycler) routinely produce single-contig bacterial chromosomes and plasmids. However, standard assembly metrics ($N_{50}$, contig count, or database-driven gene completeness tools) cannot verify whether a circular contig is physically supported across its circularization seam, whether large-scale inversions have scrambled its replichores, or whether basecaller/polisher models have introduced systematic motif errors.

`sbam` evaluates assemblies directly against the physical sequencing reads and intrinsic sequence geometry across four orthogonal dimensions:

1. **Physical Circularity & Replicon Copy Number (`sbam.struct.eval`)**: Appends a cyclic sequence buffer across contig boundaries to test continuous read spanning while filtering out natural unclipped read ends, accurately distinguishing intact circular chromosomes and plasmids from misassembled junctions and linear fragments.
2. **Replication Structure & Origin Rotation (`sbam.struct.ori`)**: Computes zero-centered cumulative GC skew across chromosomal replicons ($\ge 500\text{ kb}$) to verify bidirectional theta-replication symmetry ($\approx 180^\circ$) and measures the physical offset of *oriC* relative to coordinate $0\text{ bp}$.
3. **Empirical Base-Level Fidelity (`sbam.base.base`)**: Performs multi-threaded, indel-aware read pileups (with cyclic buffer depth folding) to mask structural noise and identify statistically significant $k$-mer motif errors via Bonferroni-corrected binomial testing.
4. **Genome Size Completeness (`sbam.report.report`)**: Evaluates total assembly length and primary chromosome proportion against expected species genome size when `-g / --genome_size` is provided.

---

## Pipeline Architecture & Module Reference

`sbam` is organized into modular Python subpackages that execute sequentially during a standard run:

| Module | Class / Function | Primary Responsibility |
| :--- | :--- | :--- |
| `sbam.cli` | `main()`, `run_summarize()` | Parses CLI arguments, manages the unified logger (`--log`), orchestrates validation stages, exports JSON/HTML outputs, cleans intermediate files (`--clean`), and aggregates multi-sample cohorts (`sbam summarize`). |
| `sbam.struct.struct` | `GenomeStructure` | Reads the input FASTA and constructs the extended cyclic reference (`cyclic_assembly.fasta`) by appending the first $B\text{ bp}$ (`--buffer-size`) of each contig to its $3'$ end. |
| `sbam.struct.mapping` | `ReadAligner` | Streams `minimap2` alignments directly into `samtools sort` and `samtools index` to produce `alignment.sorted.bam` without intermediate SAM disk I/O. |
| `sbam.struct.eval` | `JunctionEvaluator` | Quantifies contig read depth, relative copy number, and junction spanning scores at coordinate $L$ while separating true clipped/broken alignments from natural unclipped read terminations. |
| `sbam.struct.ori` | `BioVal` | Calculates sliding-window GC skew, detrends global strand bias, locates *oriC* and *ter*, computes replichore symmetry ($^\circ$), and measures the *oriC*-to-$0\text{ bp}$ rotation offset. |
| `sbam.base.base` | `BaseAccuracy` | Executes parallel `pysam` pileups, folds cyclic buffer coverage back into $[0, B)$, masks coverage anomalies ($\ge 50\text{ bp}$), and flags discordant $4, 5, 6$-mer motifs. |
| `sbam.report.report` | `DashboardBuilder` | Compiles all structural, biological, and base-level metrics into a standalone, interactive Plotly HTML dashboard (`*_sbam_report.html`). |

---

## Quick Start

### Single Isolate Evaluation

```bash
sbam -a assembly.fasta -r reads.fastq.gz -o sbam_out -t 8 -g 5.3m --clean --log
```

### Multi-Sample Cohort Summary

After running `sbam` across multiple isolates in a sequencing batch or workflow (e.g., Nextflow or Snakemake), aggregate all sample JSON outputs into a single TSV or CSV table:

```bash
sbam summarize results/*/*_sbam_metrics.json -o cohort_summary.tsv
```

---

## Documentation Contents

* **[Installation](installation.md)**: Conda/Mamba setup, external dependencies, and test suite verification.
* **[Input, Usage & Outputs](input.md)**: Complete CLI flags, parameter tuning guide, HTML dashboard guide, JSON schema, and cohort summarization.
* **[Circularity Methodology](circularity.md)**: Cyclic buffer construction, unbiased junction scoring, and plasmid classification logic.
* **[Replication Structure Methodology](replication-structure.md)**: Cumulative GC skew mathematics, replichore symmetry, and *oriC* rotation verification.
* **[Base-Level Analysis Methodology](base-level-analysis.md)**: Buffer depth folding, indel-aware concordance, structural masking, and binomial motif testing.