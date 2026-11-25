# Installation

`sbam` requires Python $\ge 3.9$ and relies on two standard bioinformatics command-line utilities (`minimap2` and `samtools`) for long-read alignment and BAM indexing. Package metadata and dependencies are defined in `pyproject.toml` (PEP 621).

---

## 1. Install via Conda / Mamba (Recommended)

Installing with Conda or Mamba automatically resolves both the Python libraries (`numpy`, `pysam`, `scipy`, `biopython`, `jinja2`) and the external C binaries (`minimap2` and `samtools`):

```bash
git clone https://github.com/riccabolla/sbam.git
cd sbam

# Create and activate the environment
conda env create -f environment.yml
conda activate sbam-env

# Install sbam in editable or standard mode (add [dev] for pytest, ruff, mypy)
pip install -e ".[dev]"
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

* `biopython>=1.81`
* `numpy>=1.24.0`
* `scipy>=1.11.0`
* `pysam>=0.22.0`
* `jinja2>=3.1.2`

---

## 3. Verify Installation

Confirm that the `sbam` CLI and external aligner dependencies are accessible:

```bash
sbam --version
sbam --help
minimap2 --version
samtools --version
```

To run the unit and BAM integration test suite:

```bash
pytest tests/ -v
```
