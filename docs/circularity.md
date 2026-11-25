# Circularity & Replicon Evaluation

Prokaryotic chromosomes and most plasmids are covalently closed circular DNA molecules. When an assembler outputs a circular contig in `FASTA` format, it must cut the circle at an arbitrary coordinate ($0\text{ bp}$ to $L\text{ bp}$).

However, assemblers frequently misreport linear fragments, unresolved repeat boundaries, or chimeric joins as "circular," or leave short terminal overlaps or gaps at the linearization seam. `sbam` empirically tests physical circularity by mapping long reads across a cyclically buffered reference.

---

## 1. Cyclic Buffer Construction (`sbam.struct.struct`)

For every contig $c$ of original length $L$, `GenomeStructure.create_cyclic_buffer()` extracts a prefix of length $B = \min(\text{buffer\_size}, L)$ from the $5'$ end of the contig and appends it to the $3'$ end:

$$S_{\text{cyclic}} = S[0 : L] + S[0 : B]$$

The resulting sequence $S_{\text{cyclic}}$ has total length $L + B$, and the original $5'/3'$ circularization seam is now located internally at coordinate $x = L$.

* If the contig is a true biological circle, long reads sequenced from molecules straddling the seam will align smoothly across coordinate $L$ into the $[L, L+B)$ buffer runway.
* If the contig is linear, chimeric, or contains an insertion/deletion artifact at the seam, reads approaching coordinate $L$ will terminate abruptly or exhibit soft/hard clipping at $x \approx L$.

---

## 2. Unbiased Junction Spanning Score (`sbam.struct.eval`)

`JunctionEvaluator.evaluate_junctions()` inspects all primary alignments (`not is_unmapped`, `not is_secondary`, `not is_supplementary`) within a boundary window around the junction coordinate $x = L$.

### Required Spanning Overlap ($M$)

To count as a **spanning read** ($N_{\text{span}}$), a read must align continuously across $x = L$ with a minimum flanking anchor of $M\text{ bp}$ on both the left and right sides of the junction:

$$\text{left\_extension} = L - \text{reference\_start} \ge M$$

$$\text{right\_extension} = \text{reference\_end} - L \ge M$$

where the dynamic overlap threshold $M$ scales with contig length:

$$M = \begin{cases} 500\text{ bp} & \text{if } L < 20,000\text{ bp} \\ 1,000\text{ bp} & \text{if } L \ge 20,000\text{ bp} \end{cases}$$

### Filtering Natural Unclipped Read Ends

A critical flaw in naive junction evaluators is treating *all* reads that cross $x = L$ with $< M\text{ bp}$ of overhang as "broken." Because physical DNA molecules shear at random positions during extraction, a read of aligned length $L_{\text{read}}$ crossing a perfectly intact junction will naturally terminate within $M\text{ bp}$ of the junction with probability:

$$P(\text{natural end within } M) \approx \frac{2M}{L_{\text{read}}}$$

For typical Nanopore read lengths ($L_{\text{read}} \approx 5.5\text{ kb}$) and $M = 1\text{ kb}$, roughly $36\%$ of valid reads crossing an intact circle end naturally within $1\text{ kb}$ of the junction, which would artificially cap the junction score at $\approx 0.64$.

To eliminate this bias, `sbam` inspects the CIGAR string (`cigartuples`) and terminal coordinates of every non-spanning read at the junction and classifies it into one of two categories:

1. **Clipped / Broken Reads ($N_{\text{broken}}$)**:
   * The read has $< M\text{ bp}$ extension on the left or right side **and** exhibits soft-clipping (`CIGAR` op `4`) or hard-clipping (`CIGAR` op `5`) of $\ge 50\text{ bp}$ on that end; **OR**
   * The read alignment starts or terminates abruptly at the junction boundary ($|\text{reference\_start} - L| \le 15\text{ bp}$ or $|\text{reference\_end} - L| \le 15\text{ bp}$).
2. **Neutral Natural Read Ends ($N_{\text{neutral}}$)**:
   * The read crosses $x = L$ with $< M\text{ bp}$ extension on one side, has **no terminal clipping** ($< 50\text{ bp}$), and does not pile up at the exact seam coordinate. These represent normal physical molecule ends and are excluded from the denominator.

### Junction Score & Status Thresholds

The junction spanning score $J$ is computed as:

$$J = \frac{N_{\text{span}}}{N_{\text{span}} + N_{\text{broken}}}$$

Contigs are assigned a circularity status based on coverage, mode, and contig size:

* **`NO_DATA`**: Average contig depth $< 5.0\times$ (triggers the `INSUFFICIENT COVERAGE` verdict if on the primary chromosome).
* **`LINEAR`**: Assigned when `--linear` mode is enabled for organisms with naturally linear chromosomes or plasmids (e.g., *Streptomyces*, *Borrelia*).
* **`PASS`**:
  * $J > 0.60$ for contigs $\ge 50\text{ kb}$
  * $J > 0.40$ for small plasmids $< 50\text{ kb}$
* **`FAIL`**: $N_{\text{span}} = 0$ or $J$ falls below the size-adjusted threshold. Any secondary contig with `FAIL` status automatically downgrades an otherwise passing assembly to `REVIEW`.

---

## 3. Depth, Copy Number & Replicon Classification

Because primary reads originating from the first $B\text{ bp}$ of a circular contig map with equal probability to $[0, B)$ and $[L, L+B)$, summing the aligned reference span (`read.reference_length`) of all primary alignments across `target_contig` and dividing by the true biological length $L$ recovers the exact unfolded average depth:

$$\text{Depth}_c = \frac{\sum \text{aligned\_bases}_c}{L_c}$$

Setting the longest contig's depth as the chromosomal baseline ($\text{Depth}_{\text{chr}}$), `sbam` calculates the relative copy number for every contig:

$$\text{Copy Number}_c = \frac{\text{Depth}_c}{\text{Depth}_{\text{chr}}}$$

Using length, circularity status, and copy number, each contig is classified automatically:

| Condition | Inferred Classification | Dashboard Evaluation |
| :--- | :--- | :--- |
| Longest contig (standard mode) | `Primary Chromosome` | `Reliable Replicon Graph` (if symmetry acceptable) |
| Longest contig (`--linear` mode) | `Linear Primary Chromosome` | `Reliable Replicon Graph` |
| Secondary contig $\ge 500\text{ kb}$ (`PASS` / `LINEAR`) | `Secondary Chromosome / Megaplasmid` | `Structurally Intact` (or replication evaluation) |
| Secondary contig $\ge 500\text{ kb}$ (`FAIL`) | `Misassembled Secondary Replicon` | `Contradicted Junction` |
| Secondary $< 500\text{ kb}$, `PASS`, $\text{Copy Number} \ge 10.0\text{c}$ | `High-Copy Circular Plasmid` | `Structurally Intact` |
| Secondary $< 500\text{ kb}$, `PASS`, $0.7\text{c} \le \text{Copy Number} < 10.0\text{c}$ | `Circular Plasmid` | `Structurally Intact` |
| Secondary $< 500\text{ kb}$, `PASS`, $\text{Copy Number} < 0.7\text{c}$ | `Sub-stoichiometric Circular Plasmid` | `Structurally Intact (<0.7c: partial plasmid loss or extraction bias)` |
| Secondary $< 500\text{ kb}$, `FAIL`, $\text{Copy Number} < 0.5\text{c}$ | `Low-Coverage Debris` | `Likely Low-Coverage Artifact` |
| Secondary $< 500\text{ kb}$, `FAIL`, $N_{\text{span}} = 0, N_{\text{broken}} = 0$ | `Linear Fragment` | `Unbridged / Linear Fragment` |
| Secondary $< 500\text{ kb}$, `FAIL`, $N_{\text{broken}} > 10$ ($\ge 10.0\text{c}$) | `Misassembled High-Copy Plasmid` | `Contradicted Junction` |
| Secondary $< 500\text{ kb}$, `FAIL`, $N_{\text{broken}} > 10$ ($< 10.0\text{c}$) | `Misassembled Junction` | `Contradicted Junction` |