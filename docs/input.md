# Input, Usage & Outputs

`sbam` takes as input a consensus genome assembly in `FASTA` format and the long reads in `FASTQ` (or `FASTA`, uncompressed or `.gz`) used to generate the assembly.

---

## Basic Usage

### Standard Single-Sample Run

```bash
sbam -a assembly.fasta -r reads.fastq.gz -o path/to/output_dir
```

### Production Run with Genome Size, Offline SVG Plots, Log File, and Cleanup

```bash
sbam -a assembly.fasta \
     -r reads.fastq.gz \
     -o path/to/output_dir \
     -t 8 \
     -g 5.3m \
     --read-type lr:hq \
     --offline \
     --clean \
     --log
```

### Cohort Summary Mode (`sbam summarize`)

When evaluating multiple isolates across a sequencing run, use `sbam summarize` to compile all generated `*_sbam_metrics.json` files into a cohort summary table (`.tsv` or `.csv`) and a cross-sample recurrent motif table:

```bash
sbam summarize path/to/results/*/*_sbam_metrics.json -o cohort_summary.tsv
```

---

## Command-Line Parameters

### Required Arguments (`sbam`)

| Flag | Argument | Description |
| :--- | :--- | :--- |
| `-a`, `--assembly` | `FILE` | Path to the polished consensus assembly (`FASTA`). |
| `-r`, `--reads` | `FILE` | Path to the long reads (`FASTQ` / `FASTA`, plain or `.gz`). |
| `-o`, `--outdir` | `DIR` | Output directory for reports, metrics, and alignments. |

### Optional Arguments (`sbam`)

| Flag | Default | Description |
| :--- | :--- | :--- |
| `-t`, `--threads` | `4` | Number of CPU threads allocated to `minimap2`, `samtools sort`, and `pysam` base pileup workers. |
| `-g`, `--genome_size`, `--genome-size` | `None` | Expected species genome size. Accepts suffixes such as `4.8m`, `4.8mbp`, `4800kb`, or raw integer base pairs `4800000`. |
| `--buffer-size` | `50000` | Length of the cyclic sequence buffer ($\text{bp}$) appended to the $3'$ end of each contig for circular read mapping ($\ge 1000\text{ bp}$). |
| `--read-type` | `map-ont` | `minimap2` alignment preset. Choices: `map-ont` (standard ONT), `lr:hq` (ONT R10.4.1 / Q20+ SUP), `map-pb` (PacBio CLR), `map-hifi` (PacBio HiFi). |
| `--min-concordance` | `0.85` | Minimum empirical read-to-assembly base concordance ($0.0\text{--}1.0$) before a position is counted as discordant ($< 85\%$ read support). Sites with $< 50\%$ read support are tracked separately as majority consensus errors. |
| `--depth-variance` | `2.5` | Fold-change threshold relative to median contig depth ($\text{median} \times 2.5$ and $\text{median} / 2.5$) used to mask collapsed repeats or coverage dips. |
| `--k-list` | `4,5,6` | Comma-separated $k$-mer sizes ($2\text{--}10$) to test for systematic motif discordance. |
| `--min-ori-length` | `500000` | Minimum contig length ($\text{bp}$) required to run chromosomal GC-skew replication symmetry and DnaA-box analysis. |
| `--linear` | `False` | Treat replicons as naturally linear (e.g., *Streptomyces*, *Borrelia*, or linear phages/plasmids) rather than requiring circular seam spanning. |
| `--offline` | `False` | Omit external Plotly CDN scripts and render all polar GC-skew architecture maps using `sbam`'s built-in zero-dependency inline SVG engine (ideal for air-gapped HPC nodes). Note: even without `--offline`, the HTML report automatically falls back to the inline SVG renderer if the browser is offline. |
| `--clean`, `--remove-intermediates` | `False` | Deletes intermediate files (`cyclic_assembly.fasta`, `alignment.sorted.bam`, and `alignment.sorted.bam.bai`) after generating the HTML and JSON reports. |
| `--log`, `--log-file` | `None` | Writes a complete unified execution log (capturing all `sbam` console output, `minimap2` stderr, and `samtools` stderr). If passed without a filename (`--log`), saves to `<outdir>/sbam.log`. |

---

## How to Choose Key Parameters

To get the most accurate and efficient results from `sbam`, it helps to understand how the optional parameters interact with your sequencing data.

### `-g, --genome_size`

`sbam` evaluates the assembly against the expected biological size of the organism. It calculates a total **Size Completeness** score, allowing it to warn you if the assembler dropped massive regions of DNA ($< 90\%$ of expected size), if accessory plasmids, uncollapsed repeats, or contamination inflated the assembly ($> 110\%$), or if the primary chromosome is fragmented ($< 85\%$ of expected genome size despite total length being within $90\%\text{--}110\%$).

It does not require an exact base-pair value, as it tolerates a $\pm 10\%$ size window before triggering a warning. You can provide a known species reference size (e.g., `-g 5.3m` for *Klebsiella pneumoniae*) or estimate it directly from your raw reads using tools like [lrge](https://github.com/mbhall88/lrge).

### `--buffer-size`

Bacterial chromosomes and plasmids are physical circles, but a `FASTA` file represents each sequence as a linear string of bases. The coordinate where a circular sequence is artificially opened is called the "seam" or junction. If a sequencing read spans this artificial breakpoint, a standard aligner mapping against a linearized reference will split the read into supplementary alignments or soft-clip the overhanging bases.

To solve this, `sbam` copies the first $B$ bases of each contig and appends them to the $3'$ end of the sequence, creating a continuous runway across the seam. Reads spanning the boundary align continuously across coordinate $L$ (the original contig length) rather than clipping.

**Rule of thumb:** Ensure `--buffer-size` is larger than your longest sequencing reads. For standard Nanopore libraries with read $N_{50} \approx 10\text{--}25\text{ kb}$, the default `50000` ($50\text{ kb}$) provides ample margin. If you are analyzing ultra-long Nanopore libraries with reads exceeding $100\text{ kb}$, increase the buffer accordingly (e.g., `--buffer-size 150000`).

### `-t, --threads`

`sbam` parallelizes two compute-intensive stages: `minimap2` read alignment and `pysam` base-level pileup profiling across all contigs. We recommend running `sbam` with **4 to 8 threads**, after which disk I/O and BAM decompression yield diminishing returns.

![Figure 1. Threads usage](figure/thread_scalability.png)

*Figure 1. CPU thread scalability across alignment and pileup stages.*

---

## Output Files

Running `sbam` produces the following files inside `--outdir`:

1. **`<assembly>_sbam_report.html`**:
   A self-contained HTML dashboard displaying:
   * **Executive Summary & Inspection Notes**: Overall rate-based assembly verdict (`ACCEPT`, `ACCEPT (WITH NOTES)`, `ACCEPT WITH CAUTION`, `REVIEW`, `FAIL`, or `INSUFFICIENT COVERAGE`), prominent inspection notes for minor findings within acceptance thresholds, and an **Assembly Overview** table covering circularity and base-level consensus QC (`<50% Err | Masked`) across **all** contigs (chromosomes and plasmids).
   * **Confidence Profile Cards**: At-a-glance status for Seam & Contig Continuity, Replication Structure, Primary Base Concordance (displaying empirical read-concordance $\text{QV}$ and consensus errors per $\text{Mb}$), and Size Completeness.
   * **Circularity Table**: Per-contig length, average depth, relative copy number (e.g., $1.0\text{c}$, $0.6\text{c}$, $4.4\text{c}$), spanning vs. clipped read counts, and junction scores.
   * **Replication Structure & Polar Plots**: Replichore symmetry ($^\circ$), predicted *oriC* and *ter* coordinates, DnaA-box enrichment ($\pm 25\text{ kb}$ observed vs. expected counts, fold enrichment, and Poisson $p$-value), *oriC* offset from $0\text{ bp}$ (`ROTATED`, `NEAR_ORIC`, `UNROTATED`), and interactive Plotly (or inline SVG) polar GC-skew maps.
   * **Base-Level Analysis**: Breakdown of structurally masked bases ($\%$), **Majority Consensus Disagreement ($< 50\%$ read support)** vs. **Minority Sub-Consensus Read Disagreement ($50\%\text{--}85\%$ read support)**, largest excluded structural regions ($\ge 50\text{ bp}$), and the **Clustered Motif Families** table with mutually exclusive site assignment and 3-state orientation classification.

2. **`<assembly>_sbam_metrics.json`**:
   A structured machine-readable JSON file containing the verdict, completeness percentage, per-contig junction metrics, replication symmetry and DnaA enrichment statistics, per-contig base-level metrics (`contig_base_metrics`), primary chromosome empirical $\text{QV}$, masked intervals, and clustered motif families.

3. **`cohort_summary.tsv` & `cohort_summary_recurrent_motifs.tsv`** *(generated via `sbam summarize`)*:
   * `cohort_summary.tsv`: Aggregates 27 structural, geometric, and base-level metrics per isolate across a study.
   * `cohort_summary_recurrent_motifs.tsv`: Aggregates motif families across all samples in the batch, calculating cohort prevalence ($\%$), total unique discordant sites, total consensus errors, and classifying families as `Recurrent Batch/Basecaller Signature` ($\ge 2$ samples) vs. `Sample-Specific Motif`.

4. **`sbam.log`** *(optional, created when `--log` is passed)*:
   A complete chronological log combining `sbam` progress messages, `minimap2` alignment statistics, and `samtools` output.

5. **Intermediate Files** *(automatically deleted when `--clean` is passed)*:
   * `cyclic_assembly.fasta`: Contigs extended with the cyclic $3'$ buffer.
   * `alignment.sorted.bam` & `alignment.sorted.bam.bai`: Sorted and indexed read alignments against the cyclic assembly.