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