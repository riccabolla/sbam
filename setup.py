from setuptools import setup, find_packages

setup(
    name="sbam",
    version="0.1.0",
    description="Structural and Base-Level analysis of microbial assemblies",
    packages=find_packages(),
    install_requires=[
        "pysam",
        "numpy",
        "scipy",
        "biopython",
        "jinja2"
    ],
    entry_points={
        "console_scripts": [
            "sbam=sbam.__main__:main",
        ],
    },
)