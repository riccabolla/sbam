import os
import subprocess
import sys


class ReadAligner:
    def __init__(self, reads_path, cyclic_fasta_path, outdir, threads=4):
        self.reads = reads_path
        self.reference = cyclic_fasta_path
        self.outdir = outdir
        self.threads = threads
        self.bam_path = os.path.join(self.outdir, "alignment.sorted.bam")
        self.minimap_log = os.path.join(self.outdir, "minimap2.log")

    def run_mapping(self, read_type="map-ont"):
        """
        Runs minimap2 and pipes output to samtools to create a sorted, indexed BAM.
        read_type: 'map-ont' for Nanopore, 'map-pb' for PacBio.
        """
        print(" > [minimap2] Aligning reads to cyclic reference...")

        minimap_cmd = [
            "minimap2",
            "-a",
            "-x",
            read_type,
            "-t",
            str(self.threads),
            self.reference,
            self.reads,
        ]

        sort_cmd = [
            "samtools",
            "sort",
            "-@",
            str(self.threads),
            "-o",
            self.bam_path,
            "-",
        ]

        try:
            with open(self.minimap_log, "w") as mm_err_log:
                p_minimap = subprocess.Popen(
                    minimap_cmd,
                    stdout=subprocess.PIPE,
                    stderr=mm_err_log,
                )
                p_sort = subprocess.Popen(
                    sort_cmd,
                    stdin=p_minimap.stdout,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )

                # Allow p_minimap to receive a SIGPIPE if p_sort exits early
                if p_minimap.stdout is not None:
                    p_minimap.stdout.close()

                _, sort_err = p_sort.communicate()
                minimap_rc = p_minimap.wait()

            if minimap_rc != 0:
                with open(self.minimap_log, "r") as mm_err_log:
                    err_tail = "".join(mm_err_log.readlines()[-20:])
                print(
                    f"[ERROR] minimap2 failed (exit code {minimap_rc}):\n{err_tail}",
                    file=sys.stderr,
                )
                sys.exit(1)

            if p_sort.returncode != 0:
                print(
                    f"[ERROR] samtools sort failed: {sort_err.decode('utf-8', errors='replace')}",
                    file=sys.stderr,
                )
                sys.exit(1)

            print(f" > [samtools] Alignment sorted and saved to: {self.bam_path}")

            index_cmd = ["samtools", "index", self.bam_path]
            subprocess.run(index_cmd, check=True)
            print(" > [samtools] BAM indexed successfully.")

            return self.bam_path

        except FileNotFoundError as e:
            print(
                f"[ERROR] Dependency missing. Ensure minimap2 and samtools are in your PATH. Details: {e}",
                file=sys.stderr,
            )
            sys.exit(1)
        except Exception as e:
            print(f"[ERROR] Alignment pipeline failed: {e}", file=sys.stderr)
            sys.exit(1)