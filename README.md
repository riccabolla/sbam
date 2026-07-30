# sbam

[![Pytest](https://github.com/riccabolla/SBAM/actions/workflows/pytest.yml/badge.svg?branch=main)](https://github.com/riccabolla/SBAM/actions/workflows/pytest.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


SBAM is quality control pipeline for long-read microbial genome assemblies. 

It maps the original sequencing reads back to the assembly to empirically prove structural circularity, validate biological replication architecture, and detect systematic sequence motif errors. 

You can find all the details in the *[wiki]()*.

## Installation

### Prerequisites

Ensure you have the following in your system's $PATH:

* Minimap2

* Samtools

### Conda

We recommend using a conda environment to prevent dependency conflicts:

```bash
conda create -n sbam-env -c bioconda sbam -y

conda activate sbam

sbam -h
```

### Source
```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/sbam.git
cd sbam

# Create environment
conda create -n sbam-env python=3.10 minimap2 samtools
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
## Parameters

```bash
-a, --assembly

[Required] Path to the consensus assembly (FASTA).

-r, --reads

[Required] Path to the long reads (FASTQ/FASTQ.gz).

-o, --outdir

[Required] Directory to save the BAM files and HTML report.

-t, --threads

Number of CPU threads to use. (Default: 4)
```