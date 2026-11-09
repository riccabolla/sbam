# Base-Level Fidelity & Motif Analysis

Even when a long-read assembly is structurally circular and geometrically symmetric, basecaller neural networks and consensus polishers can introduce systematic base-level errors—most commonly homopolymer length contractions/expansions or substitution/indel dropouts at methylated sequence motifs (such as Dcm `CCWGG` or Dam `GATC` sites).

Reference-free consensus quality estimators that rely on `FASTQ` base quality scores fail to catch these errors when the basecaller is *confidently wrong* across reads sharing the same modification context. `sbam` bypasses `FASTQ` Q-scores by measuring empirical read-to-assembly pileup concordance across the primary chromosome and testing $k$-mer motifs for statistically significant error enrichment.

---

## 1. Cyclic Buffer Depth & Match Folding (`sbam.base.base`)

When reads are aligned to the cyclically buffered reference $S_{\text{cyclic}}$ of length $L + B$, reads originating from the first $B\text{ bp}$ of the chromosome ($[0, B)$) can map equally well to $[0, B)$ or to the identical cyclic tail $[L, L+B)$. Consequently, primary read depth across $[0, B)$ is split roughly in half between the two regions.

To prevent the first $B\text{ bp}$ of the chromosome from suffering an artificial $50\%$ coverage drop, `BaseAccuracy.build_fidelity_profile()` piles up the entire buffered length $L + B$ across `-t / --threads` parallel `multiprocessing` workers and folds the cyclic tail back into the start of the chromosome:

$$\text{Depth}[0 : B] \leftarrow \text{Depth}_{\text{raw}}[0 : B] + \text{Depth}_{\text{raw}}[L : L + B]$$

$$\text{Matches}[0 : B] \leftarrow \text{Matches}_{\text{raw}}[0 : B] + \text{Matches}_{\text{raw}}[L : L + B]$$

---

## 2. Indel-Aware Empirical Concordance

In Oxford Nanopore assemblies, systematic basecalling and polishing errors predominantly manifest as $1\text{--}2\text{ bp}$ insertions or deletions rather than pure base substitutions.

In `_multi_thread()`, `sbam` extracts aligned query sequences at each chromosomal coordinate $x \in [0, L)$ using `pileupcolumn.get_query_sequences(add_indels=True)`. With `add_indels=True`:

* A read with an exact base match and no adjacent indel yields a single character (e.g., `'A'`), matching `ref_base`.
* A read carrying an insertion immediately following position $x$ yields an indel string (e.g., `'A+1A'` or `'A+2TT'`), which fails equality (`qb.upper() == ref_base`).
* A read carrying a deletion at or immediately after position $x$ yields `'*'` or `'A-1G'`, which also fails equality.

Thus, the empirical concordance at coordinate $x$ captures substitutions, insertions, and deletions simultaneously:

$$\text{Concordance}[x] = \frac{\text{Matches}[x]}{\text{Depth}[x]}$$

---

## 3. Structural Noise Masking & Interval Merging

Before scanning for sequence-specific $k$-mer errors, `sbam` masks chromosomal regions affected by large-scale structural anomalies (such as collapsed multi-copy rRNA operons, mobile IS elements, or prophage coverage dropouts) so they do not confound basecaller motif statistics.

Given the median non-zero chromosomal depth $D_{\text{med}}$, the depth variance limit $V$ (default `--depth-variance 2.5`), and minimum concordance threshold $C_{\min}$ (default `--min-concordance 0.75`), a base $x$ is masked if:

$$\text{Depth}[x] > D_{\text{med}} \times V \quad \lor \quad \text{Depth}[x] < \frac{D_{\text{med}}}{V} \quad \lor \quad \text{Concordance}[x] < C_{\min}$$

### Reporting Excluded Structural Regions

To summarize structural anomalies cleanly in the HTML report:

1. Isolated single-base/short discordance blips ($< 20\text{ bp}$) are excluded from structural interval merging so they remain available for $k$-mer motif testing without bridging healthy intervening sequence.
2. Contiguous masked blocks ($\ge 20\text{ bp}$) separated by $< 500\text{ bp}$ are merged into unified structural intervals.
3. Merged intervals with total length $\ge 50\text{ bp}$ are sorted by size (up to the top 50) and reported alongside their **Mean Depth ($\times$)**, **Mean Concordance ($\%$)**, and primary classification:
   * `High Depth (Collapsed Repeat / Multi-copy Operon)`
   * `Low Depth (Gap / Poor Coverage)`
   * `Low Concordance (Severe Basecalling Discordance)`

If the total fraction of masked bases across the primary chromosome exceeds $5.0\%$, the dashboard highlights an **Excluded Structural Regions Warning**, indicating potential fragmentation or unresolved repeats.

---

## 4. Statistical $k$-mer Motif Error Discovery

Across all unmasked regions of the primary chromosome, `BaseAccuracy.analyze_motifs()` scans sliding windows of length $k \in \{4, 5, 6\}$ (`k_list`). Any $k$-mer window overlapping a masked coordinate or containing an ambiguous `'N'` base is skipped in $O(1)$ time using a prefix-sum array (`np.cumsum(mask_arr)`).

For each $k$-mer size:

1. **Count Occurrences & Center-Base Discordance**:
   For every valid $k$-mer at index $i$, `sbam` increments its total occurrence count $n_m$ and checks whether the center nucleotide of the motif ($x_{\text{center}} = i + \lfloor k / 2 \rfloor$) has empirical concordance below $C_{\min}$ ($\text{Concordance}[x_{\text{center}}] < 0.75$), incrementing its discordant count $d_m$.
2. **Global Baseline Error Rate ($p_0$)**:
   The background discordance rate across all evaluated $k$-mer windows of size $k$ is:

   $$p_0 = \frac{\sum_{m} d_m}{\sum_{m} n_m}$$

3. **Bonferroni-Corrected One-Sided Binomial Test**:
   For every motif $m$ with at least $n_m \ge 500$ occurrences on the chromosome, `sbam` tests the null hypothesis that the motif's error rate $p_m = d_m / n_m$ is no greater than the global baseline error rate $p_0$:

   $$P(X \ge d_m) = \sum_{j=d_m}^{n_m} \binom{n_m}{j} p_0^j (1 - p_0)^{n_m - j}$$

   To control the family-wise error rate across all $M_k$ eligible motifs of length $k$ (up to $4^6 = 4,096$ hypotheses), a Bonferroni multiple-testing threshold is applied:

   $$\alpha_{\text{adj}} = \frac{0.01}{M_k}$$

4. **Effect-Size Filter**:
   A motif is flagged as a **Systematic Motif Error** if and only if:

   $$p\text{-value} < \alpha_{\text{adj}} \quad \text{and} \quad \text{Fold Increase} = \frac{p_m}{p_0} > 2.0\times$$

Flagged motifs are ranked by fold increase and displayed in both the terminal summary, the HTML dashboard, and the `*_sbam_metrics.json` output.