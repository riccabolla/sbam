# sbam

[![Pytest](https://github.com/riccabolla/SBAM/actions/workflows/pytest.yml/badge.svg?branch=main)](https://github.com/riccabolla/SBAM/actions/workflows/pytest.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


SBAM is a quality control pipeline for long-read microbial genome assemblies.

It maps the original sequencing reads back to the assembly (or reads a pre-aligned BAM) to empirically prove structural circularity, validate biological replication architecture, and detect systematic sequence motif errors.

You can find all the details in the *[wiki](https://sbam.readthedocs.io/en/latest/)*.

## Installation

### Prerequisites

Ensure you have the following in your system's $PATH:

* Minimap2
* Samtools

### Conda

We recommend using a conda environment to prevent dependency conflicts:

```bash
conda create -n sbam-env -c bioconda sbam -y

conda activate sbam-env

sbam -h
```

### Source
```bash
# Clone the repository
git clone https://github.com/riccabolla/sbam.git
cd sbam

# Create environment
conda env create -f environment.yml
conda activate sbam-env

# Install SBAM
pip install -e .

# Verify Installation
python -m sbam -h
```

## Quick Start
```bash
sbam \
    -a /path/to/assembly.fasta \
    -r /path/to/polished_reads.fastq.gz \
    -o /path/to/output_dir \
    -t 8
```

You can also re-run QC thresholds without re-aligning using a pre-aligned BAM, or aggregate cohorts:

```bash
# Re-run using an existing sorted, indexed BAM and custom sample prefix
sbam -a /path/to/assembly.fasta --bam /path/to/alignment.sorted.bam -o /path/to/output_dir --prefix isolate_01

# Aggregate cohort JSON reports into TSV + recurrent motif table
sbam summarize "results/*/*_sbam_metrics.json" -o cohort_summary.tsv
```

## Parameters

```text
-a, --assembly
  [Required] Path to the consensus assembly (FASTA).

-r, --reads
  Path to the long reads (FASTQ/FASTA, optionally gzipped). Mutually exclusive with --bam.

--bam
  Path to a pre-aligned coordinate-sorted BAM file (auto-indexed if .bai is missing). Mutually exclusive with -r/--reads.

-o, --outdir
  [Required] Directory to save outputs, diagnostic tracks, and the HTML report.

--prefix, --sample-name
  Sample prefix used for output filenames and report titles (defaults to assembly filename).

-t, --threads
  Number of CPU threads to use (Default: 4).

-g, --genome_size, --genome-size
  Expected genome size (e.g., 4.8m, 5000kb, 4800000; Default: None).

--buffer-size BUFFER_SIZE
  Size of the cyclic buffer edge (bp) for mapping (Default: 50000).

--read-type {map-ont,lr:hq,map-pb,map-hifi}
  Sequencing technology preset for minimap2 (Default: map-ont).

--min-concordance
  Minimum empirical read-to-assembly base concordance before flagging a site as discordant (Default: 0.85, must be >0.50 and <=1.0).

--depth-variance
  Fold-change threshold relative to median contig depth before structurally masking a region (Default: 2.5, must be >1.0).

--k-list
  Comma-separated k-mer sizes to test for systematic motif discordance (Default: 4,5,6).

--min-ori-length
  Minimum contig length (bp) to run chromosomal GC-skew replication analysis (Default: 500000).

--linear
  Treat replicons as naturally linear rather than circular (Default: False).

--offline
  Render all GC-skew architecture maps using SBAM's inline SVG engine without external Plotly CDN (Default: False).

--emit-tracks / --no-emit-tracks
  Write IGV/polishing diagnostic tracks (<prefix>_masked_regions.bed, <prefix>_consensus_errors.tsv, <prefix>_consensus_errors.vcf) (Default: enabled).

--clean, --remove-intermediates
  Delete intermediate alignment files created by SBAM (cyclic FASTA, BAM, and BAI) upon completion (Default: False).

--log, --log-file
  Save a unified SBAM execution log file (defaults to <outdir>/sbam.log when passed without a path).

-V, --version
  Print SBAM version and exit.
```
