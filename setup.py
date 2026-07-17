from setuptools import setup, find_packages

setup(
    name="sbam",
    version="1.0.0",
    description="Structural and Base-Level Assessment of Microbes",
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