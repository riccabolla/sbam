"""Long-read alignment wrapper using minimap2 and samtools."""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile

logger = logging.getLogger("sbam")


class AlignmentError(RuntimeError):
    """Raised when read alignment or BAM indexing fails."""


def _read_tail(path: str, lines: int = 25) -> str:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return "".join(fh.readlines()[-lines:])
    except OSError:
        return ""


def _log_tool_output(label: str, path: str) -> None:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            content = fh.read().strip()
        if content:
            logger.debug("--- [%s stderr] ---\n%s", label, content)
    except OSError:
        pass


class ReadAligner:
    """Aligns long reads to a cyclic or linear reference and writes a sorted, indexed BAM."""

    def __init__(
        self,
        reads_path: str,
        reference_path: str,
        outdir: str,
        threads: int = 4,
        log_path: str | None = None,
        is_cyclic: bool = True,
    ):
        self.reads = reads_path
        self.reference = reference_path
        self.outdir = outdir
        self.threads = threads
        self.bam_path = os.path.join(self.outdir, "alignment.sorted.bam")
        self.log_path = log_path
        self.is_cyclic = is_cyclic

    def run_mapping(self, read_type: str = "map-ont") -> str:
        """
        Run minimap2 (`--eqx --secondary=no`) piped to `samtools sort --write-index`.
        Raises `AlignmentError` on missing dependencies or non-zero exit status.
        """
        ref_kind = "cyclic" if self.is_cyclic else "linear"
        logger.info(" > [minimap2] Aligning reads to %s reference...", ref_kind)

        minimap_cmd = [
            "minimap2",
            "-a",
            "--eqx",
            "--secondary=no",
            "-x",
            read_type,
            "-t",
            str(self.threads),
            self.reference,
            self.reads,
        ]

        bai_path = f"{self.bam_path}.bai"
        sort_cmd = [
            "samtools",
            "sort",
            "--write-index",
            "-@",
            str(self.threads),
            "-o",
            f"{self.bam_path}##idx##{bai_path}",
            "-",
        ]

        mm_err_path = None
        sort_err_path = None
        p_minimap = None
        p_sort = None

        try:
            with tempfile.NamedTemporaryFile(
                mode="w+", suffix="_minimap2.log", delete=False, encoding="utf-8"
            ) as mm_err_file, tempfile.NamedTemporaryFile(
                mode="w+", suffix="_samtools_sort.log", delete=False, encoding="utf-8"
            ) as sort_err_file:
                mm_err_path = mm_err_file.name
                sort_err_path = sort_err_file.name

                try:
                    p_minimap = subprocess.Popen(
                        minimap_cmd,
                        stdout=subprocess.PIPE,
                        stderr=mm_err_file,
                    )
                except FileNotFoundError as exc:
                    raise AlignmentError(
                        f"Dependency missing (minimap2). Ensure minimap2 is in PATH. Details: {exc}"
                    ) from exc

                try:
                    p_sort = subprocess.Popen(
                        sort_cmd,
                        stdin=p_minimap.stdout,
                        stdout=subprocess.PIPE,
                        stderr=sort_err_file,
                    )
                except FileNotFoundError as exc:
                    p_minimap.terminate()
                    p_minimap.wait()
                    raise AlignmentError(
                        f"Dependency missing (samtools). Ensure samtools is in PATH. Details: {exc}"
                    ) from exc

                if p_minimap.stdout is not None:
                    p_minimap.stdout.close()

                p_sort.communicate()
                minimap_rc = p_minimap.wait()

                mm_err_file.flush()
                sort_err_file.flush()
                _log_tool_output("minimap2", mm_err_path)
                _log_tool_output("samtools sort", sort_err_path)

                if minimap_rc != 0 or p_sort.returncode != 0:
                    mm_tail = _read_tail(mm_err_path)
                    sort_tail = _read_tail(sort_err_path)
                    raise AlignmentError(
                        f"Alignment failed (minimap2={minimap_rc}, samtools={p_sort.returncode}).\n"
                        f"{mm_tail}\n{sort_tail}".strip()
                    )

            logger.info(" > [samtools] Alignment sorted and saved to: %s", self.bam_path)

            if not os.path.exists(bai_path):
                index_cmd = ["samtools", "index", "-@", str(self.threads), self.bam_path]
                idx_proc = subprocess.run(
                    index_cmd,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if idx_proc.stderr:
                    logger.debug("--- [samtools index stderr] ---\n%s", idx_proc.stderr.strip())
                if idx_proc.returncode != 0:
                    raise AlignmentError(
                        f"samtools index failed ({idx_proc.returncode}): {idx_proc.stderr.strip()}"
                    )
            logger.info(" > [samtools] BAM indexed successfully.")

            return self.bam_path

        except AlignmentError:
            raise
        except Exception as exc:
            if p_minimap is not None and p_minimap.poll() is None:
                p_minimap.terminate()
                p_minimap.wait()
            if p_sort is not None and p_sort.poll() is None:
                p_sort.terminate()
                p_sort.wait()
            raise AlignmentError(f"Alignment pipeline failed: {exc}") from exc
        finally:
            for tmp_path in (mm_err_path, sort_err_path):
                if tmp_path is not None:
                    try:
                        if os.path.exists(tmp_path):
                            os.remove(tmp_path)
                    except OSError:
                        pass
