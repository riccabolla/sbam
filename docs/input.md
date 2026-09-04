## Input

``sbam`` will take as input the consensu assembly in ``fasta`` format with the ``fastq`` used to generate the file. 

### Basic usage

```python
sbam -a assembly.fasta -r read.fastq.gz -o path/output/dir
```

### Parameters

#### Input/Output

``-a assembly.fasta, --assembly assembly.fasta
``

Input assembly file

``
-r read.fastq.gz, --read read.fastq.gz
``

Input ``fastq`` file

``
-o outputdir, --outdir outputdir
``

Path to the output directory

#### Optional

``
-t cpus, --threads cpus
``

Number of cpus provided

``
-g genome_size, --genome_size genome_size
``

Expected genome size based on the specie

``
--buffer-size buffer_size
``

Size of the cyclic buffer size for mapping

``
--read-type read_type
``

Specify sequencing used, required by minimap (default is map-ont)


#### How to choose

To get the most accurate and efficient results from SBAM, it helps to understand how the optional parameters affect the pipeline.

``-g, --genome_size``

SBAM evaluates the assembly against expected biological size. It calculates a total Genome Completeness score, allowing it to warn you if the assembler dropped massive regions of DNA or if the main chromosome is severely fragmented compared to the expected size.

It does not expect an exact value, since it tolerates a 10% of difference before calling a warn.

You can provide it by yourself if you already know the common genome size for the analyzed specie, or you can use tools like [lrge](https://github.com/mbhall88/lrge) to calculate genome size based on ``fastq``.

``---buffer-size``

Bacterial chromosomes are physical circles, but a ``FASTA`` file represents a sequence as a linear string of bases. The point where a circular sequence is artificially linearized is often called the "seam" or breakpoint. If a sequencing read spans this artificial breakpoint, a standard aligner mapping against the linearized reference may represent the read as two separate alignment segments or as a primary alignment with a soft-clipped portion.

To address this, SBAM copies the first N bases of the assembly and appends them to the end of the sequence, creating a continuous "runway" for reads that cross the seam. This extended sequence acts as a cyclic buffer, allowing reads spanning the artificial breakpoint to align continuously rather than being treated as split or partially clipped alignments

It is important that the buffer size is larger than your longest sequencing reads. Indeed. if your reads average 20kb, the default buffer of 50000 (50kb) provides plenty of room. However, if you are using ultra-long Oxford Nanopore sequencing (where reads can easily reach 100kb+), you must increase this (e.g., --buffer-size 150000). If the buffer is too small, long reads will clip, and SBAM will falsely flag the chromosome as broken!

``-t, --threads``

SBAM uses threads during two steps: the minimap2 alignment and the base-level motif profiling.

When possible, we recommend to run it with 4 to 8 threads, after that value there is a progressive diminishing returns.
