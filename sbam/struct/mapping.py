import os
import subprocess
import sys
import tempfile


class ReadAligner:
    def __init__(self, reads_path, cyclic_fasta_path, outdir, threads=4, log_path=None):
        self.reads = reads_path
        self.reference = cyclic_fasta_path
        self.outdir = outdir
        self.threads = threads
        self.bam_path = os.path.join(self.outdir, "alignment.sorted.bam")
        self.log_path = log_path

    def run_mapping(self, read_type="map-ont"):
        """
        Runs minimap2 and pipes output to samtools to create a sorted, indexed BAM.
        If self.log_path is set, appends minimap2 and samtools stderr to the unified log file;
        otherwise uses a temporary file that is automatically removed on completion.
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

        temp_err_file = None
        try:
            if self.log_path:
                err_stream = open(self.log_path, "a", encoding="utf-8")
                err_stream.write("\n--- [minimap2 / samtools alignment log] ---\n")
                err_stream.flush()
            else:
                temp_err_file = tempfile.NamedTemporaryFile(
                    mode="w+", suffix="_align.log", delete=False, encoding="utf-8"
                )
                err_stream = temp_err_file

            try:
                p_minimap = subprocess.Popen(
                    minimap_cmd,
                    stdout=subprocess.PIPE,
                    stderr=err_stream,
                )
                p_sort = subprocess.Popen(
                    sort_cmd,
                    stdin=p_minimap.stdout,
                    stdout=subprocess.PIPE,
                    stderr=err_stream,
                )

                if p_minimap.stdout is not None:
                    p_minimap.stdout.close()

                p_sort.communicate()
                minimap_rc = p_minimap.wait()
            finally:
                err_stream.flush()
                if self.log_path:
                    err_stream.close()

            if minimap_rc != 0 or p_sort.returncode != 0:
                err_source = self.log_path if self.log_path else temp_err_file.name
                with open(err_source, "r", encoding="utf-8", errors="replace") as ef:
                    err_tail = "".join(ef.readlines()[-25:])
                print(
                    f"[ERROR] Alignment failed (minimap2={minimap_rc}, samtools={p_sort.returncode}):\n{err_tail}",
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
        finally:
            if temp_err_file is not None:
                try:
                    temp_err_file.close()
                    if os.path.exists(temp_err_file.name):
                        os.remove(temp_err_file.name)
                except OSError:
                    pass