# Circularity, Linear Termini & Replicon Evaluation

Prokaryotic chromosomes and most plasmids are covalently closed circular DNA molecules, though notable lineages (e.g., *Streptomyces*, *Borrelia*) and linear phages/plasmids possess naturally linear replicons. When an assembler outputs a circular contig in `FASTA` format, it must cut the circle at an arbitrary coordinate ($0\text{ bp}$ to $L\text{ bp}$).

However, assemblers frequently misreport linear fragments, unresolved repeat boundaries, or chimeric joins as "circular," or leave short terminal overlaps or gaps at the linearization seam. `sbam` empirically tests physical circularity by mapping long reads across a cyclically buffered reference—or, when `--linear` is passed, bypasses buffering and evaluates $5'$ and $3'$ linear termini directly.

---

## 1. Reference Preparation (`sbam.struct.struct`)

`GenomeStructure` filters out any $0\text{-bp}$ FASTA records (logging a warning so empty records never create length mismatches downstream) and prepares the alignment reference via `prepare_reference(buffer_size, linear)`:

* **Circular Mode (default)**: `GenomeStructure.create_cyclic_buffer()` extracts a prefix of length $B = \min(\text{buffer\_size}, L)$ from the $5'$ end of each contig $c$ (of original length $L$) and appends it to the $3'$ end:

  $$S_{\text{cyclic}} = S[0 : L] + S[0 : B]$$

  The resulting sequence $S_{\text{cyclic}}$ has total length $L + B$, and the original $5'/3'$ circularization seam sits internally at coordinate $x = L$:
  * If the contig is a true biological circle, long reads straddling the seam align smoothly across coordinate $L$ into the $[L, L+B)$ buffer runway.
  * If the contig is linear, chimeric, or contains a seam indel/repeat artifact, reads approaching coordinate $L$ terminate abruptly or exhibit soft/hard clipping near $x = L$.

* **Linear Mode (`--linear`)**: Cyclic buffer creation is bypassed completely and reads are aligned directly against the original FASTA (`is_cyclic=False`).

---

## 2. Single-Pass Continuity & Junction Scoring (`sbam.struct.eval`)

`JunctionEvaluator.evaluate_junctions()` auto-detects BAM contig naming (`<id>_cyclic` vs. `<id>`) from the BAM header and scans each contig in a **single BAM pass**, skipping unmapped, secondary, supplementary (`0x800`), QC-fail, and duplicate alignments (`is_unmapped`, `is_secondary`, `is_supplementary`, `is_qcfail`, `is_duplicate`).

### Circular Junction Mode: Required Spanning Overlap ($M$)

To count as a **spanning read** ($N_{\text{span}}$) across the circular seam $x = L$, a read must align continuously across $x = L$ with a minimum flanking anchor of $M\text{ bp}$ on both the left and right sides of the junction:

$$\text{left\_extension} = L - \text{reference\_start} \ge M, \qquad \text{right\_extension} = \text{reference\_end} - L \ge M$$

where the dynamic overlap threshold $M$ scales cleanly even for very small contigs ($< 500\text{ bp}$):

$$M = \begin{cases} \min(500,\; \max(1, \lfloor L / 2 \rfloor)) & \text{if } L < 20,000\text{ bp} \\ 1,000\text{ bp} & \text{if } L \ge 20,000\text{ bp} \end{cases}$$

### Filtering Natural Unclipped Read Ends & Distant Adapter Clips

A critical flaw in naive junction evaluators is treating *all* reads that cross $x = L$ with $< M\text{ bp}$ of overhang as "broken." Because physical DNA molecules shear at random positions during extraction, a read of aligned length $L_{\text{read}}$ crossing a perfectly intact junction will naturally terminate within $M\text{ bp}$ of the junction with probability:

$$P(\text{natural end within } M) \approx \frac{2M}{L_{\text{read}}}$$

For typical Nanopore read lengths ($L_{\text{read}} \approx 5.5\text{ kb}$) and $M = 1\text{ kb}$, roughly $36\%$ of valid reads crossing an intact circle end naturally within $1\text{ kb}$ of the junction. Conversely, a broken junction with a local indel or repeat divergence may cause reads to soft-clip $16\text{--}500\text{ bp}$ away from $x = L$, whereas untrimmed reads may carry terminal adapter clips $501\text{--}999\text{ bp}$ away.

To handle both cases accurately, `sbam` inspects reads overlapping the junction window $[L - W_{\text{clip}},\, L + W_{\text{clip}}]$ (default `clip_window = 500 bp`) and enforces clip proximity:

1. **Clipped / Broken Reads ($N_{\text{broken}}$)**:
   * The read has $< M\text{ bp}$ left extension, starts within `clip_window` of the seam ($|\text{reference\_start} - L| \le 500\text{ bp}$), and exhibits leading soft/hard clipping (`CIGAR` op `4` or `5`) $\ge 50\text{ bp}$; **OR**
   * The read has $< M\text{ bp}$ right extension, ends within `clip_window` of the seam ($|\text{reference\_end} - L| \le 500\text{ bp}$), and exhibits trailing soft/hard clipping $\ge 50\text{ bp}$; **OR**
   * The read alignment starts or terminates abruptly within `boundary_tol = 15 bp` of the seam ($|\text{reference\_start} - L| \le 15\text{ bp}$ or $|\text{reference\_end} - L| \le 15\text{ bp}$).
2. **Neutral Natural Read Ends ($N_{\text{neutral}}$)**:
   * Reads overlapping the seam boundary ($[L - 15,\, L + 15]$) with $< M\text{ bp}$ overhang on one side, no local junction clipping, and no abrupt boundary termination (including reads with distant terminal adapter clips $> 500\text{ bp}$ away from $L$). These represent normal physical molecule ends and are excluded from the denominator.

### Linear Replicon Mode (`--linear`)

When `--linear` is active and contigs are unbuffered (`not contig_is_cyclic`), `JunctionEvaluator` evaluates physical continuity at the $5'$ ($0\text{ bp}$) and $3'$ ($L\text{ bp}$) chromosome/plasmid termini instead of an artificial seam:

* Reads reaching within `clip_window = 500 bp` of coordinate $0$ or $L$ are tallied (`left_terminal_reads`, `right_terminal_reads`).
* Terminal alignments with leading/trailing clips $\ge 50\text{ bp}$ at the contig ends are counted as `broken_reads`, while unclipped reads reaching within $15\text{ bp}$ of $0$ or $L$ count toward clean terminal support (`spanning_reads`).
* Contigs with depth $\ge 5.0\times$ receive status `LINEAR` and a terminal continuity score $N_{\text{span}} / (N_{\text{span}} + N_{\text{broken}})$ (or `1.0` if no clipped terminal reads are observed).

### Continuity Score & Status Thresholds

The continuity score $J$ is computed as:

$$J = \frac{N_{\text{span}}}{N_{\text{span}} + N_{\text{broken}}}$$

Contigs are assigned a status based on coverage, mode, and contig size:

* **`NO_DATA`**: Average contig depth $< 5.0\times$ (triggers `INSUFFICIENT COVERAGE` if on the primary chromosome).
* **`LINEAR`**: Assigned when `--linear` mode is enabled and average depth $\ge 5.0\times$.
* **`PASS`** *(circular mode)*:
  * $J > 0.60$ for contigs $\ge 50\text{ kb}$
  * $J > 0.40$ for small plasmids $< 50\text{ kb}$
* **`FAIL`** *(circular mode)*: $N_{\text{span}} = 0$ or $J$ falls below the size-adjusted threshold. Any secondary contig with `FAIL` status automatically downgrades an otherwise passing assembly to `REVIEW`.

---

## 3. Depth, Copy Number & Replicon Classification

Because primary reads originating from the first $B\text{ bp}$ of a circular contig map with equal probability to $[0, B)$ and $[L, L+B)$, summing the aligned reference span (`read.reference_length`) of all primary alignments across `target_contig` and dividing by the true biological length $L$ recovers the exact unfolded average depth:

$$\text{Depth}_c = \frac{\sum \text{aligned\_bases}_c}{L_c}$$

Setting the longest contig's depth as the chromosomal baseline ($\text{Depth}_{\text{chr}}$), `sbam` calculates the relative copy number for every contig:

$$\text{Copy Number}_c = \frac{\text{Depth}_c}{\text{Depth}_{\text{chr}}}$$

Using length, status, and copy number, each contig is classified automatically:

| Condition | Inferred Classification | Dashboard Evaluation |
| :--- | :--- | :--- |
| Longest contig (circular mode) | `Primary Chromosome` | `Reliable Replicon Graph` (if symmetry acceptable) |
| Longest contig (`--linear` mode) | `Linear Primary Chromosome` | `Reliable Replicon Graph` (if symmetry acceptable) |
| Secondary contig $\ge 500\text{ kb}$ (`PASS` / `LINEAR`) | `Secondary Chromosome / Megaplasmid` | `Structurally Intact` (or replication evaluation) |
| Secondary contig $\ge 500\text{ kb}$ (`FAIL`) | `Misassembled Secondary Replicon` | `Contradicted Junction` |
| Secondary $< 500\text{ kb}$, `PASS`, $\text{Copy Number} \ge 10.0\text{c}$ | `High-Copy Circular Plasmid` | `Structurally Intact` |
| Secondary $< 500\text{ kb}$, `PASS`, $0.7\text{c} \le \text{Copy Number} < 10.0\text{c}$ | `Circular Plasmid` | `Structurally Intact` |
| Secondary $< 500\text{ kb}$, `PASS`, $\text{Copy Number} < 0.7\text{c}$ | `Sub-stoichiometric Circular Plasmid` | `Structurally Intact (<0.7c: partial plasmid loss or extraction bias)` |
| Secondary $< 500\text{ kb}$, `LINEAR` | `Linear Plasmid / Replicon` | `Structurally Intact` |
| Secondary $< 500\text{ kb}$, `FAIL`, $\text{Copy Number} < 0.5\text{c}$ | `Low-Coverage Debris` | `Likely Low-Coverage Artifact` |
| Secondary $< 500\text{ kb}$, `FAIL`, $N_{\text{span}} = 0, N_{\text{broken}} = 0$ | `Linear Fragment` | `Unbridged / Linear Fragment` |
| Secondary $< 500\text{ kb}$, `FAIL`, $N_{\text{broken}} > 10$ ($\ge 10.0\text{c}$) | `Misassembled High-Copy Plasmid` | `Contradicted Junction` |
| Secondary $< 500\text{ kb}$, `FAIL`, $N_{\text{broken}} > 10$ ($< 10.0\text{c}$) | `Misassembled Junction` | `Contradicted Junction` |
