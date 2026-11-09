# Input, Usage & Outputs

`sbam` takes as input a consensus genome assembly in `FASTA` format and the long reads in `FASTQ` (or `FASTA`, uncompressed or `.gz`) used to generate the assembly.

---

## Basic Usage

### Standard Single-Sample Run

```bash
sbam -a assembly.fasta -r reads.fastq.gz -o path/to/output_dir
```

### Production Run with Genome Size, Log File, and Cleanup

```bash
sbam -a assembly.fasta \
     -r reads.fastq.gz \
     -o path/to/output_dir \
     -t 8 \
     -g 5.3m \
     --read-type lr:hq \
     --clean \
     --log
```

### Cohort Summary Mode (`sbam summarize`)

When evaluating dozens of isolates across a study, use `sbam summarize` to compile all generated `*_sbam_metrics.json` files into a single tab-separated (`.tsv`) or comma-separated (`.csv`) summary table:

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
| `--min-concordance` | `0.75` | Minimum empirical read-to-assembly base concordance ($0.0 - 1.0$) before a position is flagged as discordant. |
| `--depth-variance` | `2.5` | Fold-change threshold relative to median contig depth ($\text{median} \times 2.5$ and $\text{median} / 2.5$) used to mask collapsed repeats or coverage dropouts. |
| `--min-ori-length` | `500000` | Minimum contig length ($\text{bp}$) required to run bidirectional GC-skew replication symmetry analysis. |
| `--clean`, `--remove-intermediates` | `False` | Deletes intermediate files (`cyclic_assembly.fasta`, `alignment.sorted.bam`, and `alignment.sorted.bam.bai`) after generating the HTML and JSON reports. |
| `--log`, `--log-file` | `None` | Writes a complete unified execution log (capturing all `sbam` console output, `minimap2` stderr, and `samtools` stderr). If passed without a filename (`--log`), saves to `<outdir>/sbam.log`. |

---

## How to Choose Key Parameters

To get the most accurate and efficient results from `sbam`, it helps to understand how the optional parameters interact with your sequencing data.

### `-g, --genome_size`

`sbam` evaluates the assembly against the expected biological size of the organism. It calculates a total **Genome Completeness** score, allowing it to warn you if the assembler dropped massive regions of DNA ($< 90\%$ of expected size), if uncollapsed haplotypes or contamination inflated the assembly ($> 110\%$), or if the primary chromosome is severely fragmented ($< 85\%$ of expected genome size despite total length being normal).

It does not require an exact base-pair value, as it tolerates a $\pm 10\%$ biological accessory-genome variance before triggering a warning. You can provide a known species reference size (e.g., `-g 5.3m` for *Klebsiella pneumoniae*) or estimate it directly from your raw reads using tools like [lrge](https://github.com/mbhall88/lrge).

### `--buffer-size`

Bacterial chromosomes and plasmids are physical circles, but a `FASTA` file represents each sequence as a linear string of bases. The coordinate where a circular sequence is artificially opened is called the "seam" or junction. If a sequencing read spans this artificial breakpoint, a standard aligner mapping against a linearized reference will split the read into two supplementary alignments or soft-clip the overhanging bases.

To solve this, `sbam` copies the first $N$ bases of each contig and appends them to the $3'$ end of the sequence, creating a continuous runway across the seam. Reads spanning the boundary align continuously across coordinate $L$ (the original contig length) rather than clipping.

**Rule of thumb:** Ensure `--buffer-size` is larger than your longest sequencing reads. For standard Nanopore libraries with read $N_{50} \approx 10\text{--}25\text{ kb}$, the default `50000` ($50\text{ kb}$) provides ample margin. If you are analyzing ultra-long Nanopore libraries with reads exceeding $100\text{ kb}$, increase the buffer accordingly (e.g., `--buffer-size 150000`).

### `-t, --threads`

`sbam` parallelizes two compute-intensive stages: `minimap2` read alignment and `pysam` base-level pileup profiling. We recommend running `sbam` with **4 to 8 threads**, after which disk I/O and BAM decompression yield diminishing returns.

![Figure 1. Threads usage](figure/thread_scalability.png)

*Figure 1. CPU thread scalability across alignment and pileup stages.*

---

## Output Files

Running `sbam` produces the following files inside `--outdir`:

1. **`<assembly>_sbam_report.html`**:
   A self-contained, interactive HTML dashboard displaying:
   * **Executive Summary**: Overall assembly verdict (`ACCEPT`, `REVIEW`, `ACCEPT WITH CAUTION`, or `FAIL`) and a per-contig overview table with inferred replicon types and circularity evaluations.
   * **Confidence Profile Cards**: At-a-glance status for Structural Integrity, Replication Structure, Base-Level Fidelity, and Genome Completeness.
   * **Circularity Table**: Per-contig length, average depth, relative copy number (e.g., $1.0\text{c}$, $2.8\text{c}$), spanning vs. clipped read counts, and junction scores.
   * **Replication Structure & Polar Plots**: Replichore symmetry ($^\circ$), predicted *oriC* and *ter* coordinates, *oriC* offset from $0\text{ bp}$, and interactive Plotly polar GC-skew maps.
   * **Base-Level Analysis**: Total masked bases ($\%$), a table of the largest excluded structural regions ($\ge 50\text{ bp}$, including mean depth and mean concordance), and any statistically significant discordant $k$-mer motifs.

2. **`<assembly>_sbam_metrics.json`**:
   A structured machine-readable JSON file containing the verdict, completeness percentage, per-contig junction metrics, replication symmetry coordinates, masked intervals, and systematic error motifs. Used as input for `sbam summarize` or custom workflow parsers.

3. **`cohort_summary.tsv` / `.csv`** *(generated via `sbam summarize`)*:
   A multi-sample summary table containing 19 columns per isolate:
   `sample`, `verdict`, `num_contigs`, `total_length_bp`, `completeness_pct`, `primary_contig_id`, `primary_length_bp`, `primary_depth_x`, `primary_junction_score`, `primary_circularity_status`, `primary_symmetry_deg`, `primary_oric_offset_kb`, `primary_replication_status`, `circular_plasmids_count`, `failed_contigs_count`, `masked_bases_bp`, `masked_pct`, `systematic_error_motifs_count`, and `top_error_motifs`.

4. **`sbam.log`** *(optional, created when `--log` is passed)*:
   A complete chronological log combining `sbam` progress messages, `minimap2` alignment statistics, and `samtools` output.

5. **Intermediate Files** *(automatically deleted when `--clean` is passed)*:
   * `cyclic_assembly.fasta`: Contigs extended with the cyclic $3'$ buffer.
   * `alignment.sorted.bam` & `alignment.sorted.bam.bai`: Sorted and indexed read alignments against the cyclic assembly.