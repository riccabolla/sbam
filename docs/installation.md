# Installation

`sbam` requires Python $\ge 3.9$ and relies on two standard bioinformatics command-line utilities (`minimap2` and `samtools`) for long-read alignment and BAM indexing.

---

## 1. Install via Conda / Mamba (Recommended)

Installing with Conda or Mamba automatically resolves both the Python libraries and the external C binaries (`minimap2` and `samtools`):

```bash
git clone https://github.com/riccabolla/sbam.git
cd sbam

# Create and activate the environment
conda env create -f environment.yml
conda activate sbam

# Install sbam in editable or standard mode
pip install -e .
```

---

## 2. Install via `pip`

If you already have `minimap2` and `samtools` available in your system `PATH`, you can install `sbam` directly via `pip`:

```bash
git clone https://github.com/riccabolla/sbam.git
cd sbam
pip install .
```

### Python Dependencies

* `biopython`
* `numpy`
* `scipy`
* `pysam`
* `jinja2`

---

## 3. Verify Installation

Confirm that the `sbam` CLI and external aligner dependencies are accessible:

```bash
sbam --help
minimap2 --version
samtools --version
```

To run the unit test suite:

```bash
pytest
```