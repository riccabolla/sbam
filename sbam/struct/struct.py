import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class GenomeArchitecture:
    def __init__(self, fasta_path, outdir):
        self.fasta_path = fasta_path
        self.outdir = outdir
        # Load all contigs in memory 
        self.records = list(SeqIO.parse(self.fasta_path, "fasta"))
        
    def create_cyclic_buffer(self, buffer_size=50000):
        """
        Creates a new FASTA where the first 'buffer_size' bases 
        are appended to the end of each contig.
        """
        cyclic_records = []
        
        for record in self.records:
            seq_len = len(record.seq)
            
            # Dynamic buffer: If a plasmid is smaller than the requested buffer,
            # sbam just buffer its entire length to avoid index out-of-bounds.
            actual_buffer = min(buffer_size, seq_len)
            
            if actual_buffer == 0:
                continue
            
            # Extract the front buffer and append to the back
            buffer_seq = record.seq[:actual_buffer]
            cyclic_seq = record.seq + buffer_seq
            
            # Create a new sequence record with metadata tracking the original length
            cyclic_record = SeqRecord(
                cyclic_seq,
                id=f"{record.id}_cyclic",
                description=f"original_length={seq_len} buffer={actual_buffer}"
            )
            cyclic_records.append(cyclic_record)
            
        # Write the expanded sequences to a new file for minimap2
        out_fasta = os.path.join(self.outdir, "cyclic_assembly.fasta")
        SeqIO.write(cyclic_records, out_fasta, "fasta")
        
        return out_fasta, {rec.id: len(rec.seq) for rec in self.records}