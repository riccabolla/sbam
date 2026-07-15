import subprocess
import os
import sys

class ReadAligner:
    def __init__(self, reads_path, cyclic_fasta_path, outdir, threads=4):
        self.reads = reads_path
        self.reference = cyclic_fasta_path
        self.outdir = outdir
        self.threads = threads
        self.bam_path = os.path.join(self.outdir, "alignment.sorted.bam")

    def run_mapping(self, read_type="map-ont"):
        """
        Runs minimap2 and pipes output to samtools to create a sorted, indexed BAM.
        read_type: 'map-ont' for Nanopore, 'map-pb' for PacBio.
        """
        print(f" > [minimap2] Aligning reads to cyclic reference...")
        
        # minimap2 command
        minimap_cmd = [
            "minimap2", "-a", "-x", read_type, "-t", str(self.threads),
            self.reference, self.reads
        ]
        
        # samtools sort 
        sort_cmd = [
            "samtools", "sort", "-@", str(self.threads), "-o", self.bam_path, "-"
        ]
        
        try:
            # Execute minimap2 | samtools sort
            p_minimap = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p_sort = subprocess.Popen(sort_cmd, stdin=p_minimap.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Allow p_minimap to receive a SIGPIPE if p_sort exits
            p_minimap.stdout.close()
            
            # Wait for samtools sort to finish and capture errors
            _, sort_err = p_sort.communicate()
            
            if p_sort.returncode != 0:
                print(f"[ERROR] samtools sort failed: {sort_err.decode('utf-8')}", file=sys.stderr)
                sys.exit(1)
                
            print(f" > [samtools] Alignment sorted and saved to: {self.bam_path}")
            
            # Create BAM index
            index_cmd = ["samtools", "index", self.bam_path]
            subprocess.run(index_cmd, check=True)
            print(f" > [samtools] BAM indexed successfully.")
            
            return self.bam_path
            
        except FileNotFoundError as e:
            print(f"[ERROR] Dependency missing. Ensure minimap2 and samtools are in your PATH. Details: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"[ERROR] Alignment pipeline failed: {e}")
            sys.exit(1)