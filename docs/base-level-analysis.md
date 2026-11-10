# Base-Level Concordance, Empirical QV & Motif Family Analysis

Even when a long-read assembly is structurally circular and geometrically symmetric, basecaller models and consensus polishers can leave consensus-level errors (such as homopolymer indels or uncollapsed repeat variants) or exhibit sub-consensus read disagreement at specific sequence contexts (such as methylated motifs or G/C-rich homopolymers).

`sbam` bypasses `FASTQ` quality scores by profiling empirical read-to-assembly pileup concordance across **all contigs** in the assembly and performing a statistically rigorous motif screen on the primary chromosome.

---

## 1. All-Contig Profiling & Cyclic Buffer Depth Folding (`sbam.base.base`)

When reads are aligned to the cyclically buffered reference $S_{\text{cyclic}}$ of length $L + B$, reads originating from the first $B\text{ bp}$ of a circular contig ($[0, B)$) map with equal probability to $[0, B)$ or to the identical cyclic tail $[L, L+B)$.

To prevent an artificial $50\%$ coverage drop across the first $B\text{ bp}$, `BaseAccuracy` piles up the full buffered length $L + B$ across `-t / --threads` parallel `multiprocessing` workers (with `max_depth=100000` so high-copy plasmids are never truncated at `pysam`'s default $8,000\times$ limit) and folds the cyclic tail back onto $[0, B)$:

$$\text{Depth}[0 : B] \leftarrow \text{Depth}_{\text{raw}}[0 : B] + \text{Depth}_{\text{raw}}[L : L + B]$$

$$\text{Matches}[0 : B] \leftarrow \text{Matches}_{\text{raw}}[0 : B] + \text{Matches}_{\text{raw}}[L : L + B]$$

Every contig in the FASTA (primary chromosome and all secondary plasmids/replicons) is profiled, and its majority consensus error count ($< 50\%$ read support), error rate per $\text{Mb}$, and structurally masked percentage are reported in the **Assembly Overview** table.

---

## 2. Indel-Aware Concordance Without Deletion Double-Counting

In Oxford Nanopore assemblies, systematic basecalling and polishing errors predominantly manifest as $1\text{--}2\text{ bp}$ insertions or deletions rather than pure substitutions.

In `_multi_thread()`, `sbam` extracts aligned query sequences at each coordinate $x \in [0, L)$ using `pileupcolumn.get_query_sequences(add_indels=True)`. In `pysam`, `add_indels=True` attaches insertion tokens (e.g., `'A+1T'`) and deletion lookahead tokens (e.g., `'A-1G'`) to the base *preceding* the indel, while also placing a deletion placeholder `'*'` at the deleted coordinate itself.

To prevent a single deletion from being double-counted at both the preceding base and the deleted base:

* `sbam` strips any deletion lookahead suffix (`qb.split("-")[0].upper() == ref_base`) at the preceding coordinate.
* Thus, **deletions** fail equality once at the deleted coordinate (`'*'`), **insertions** fail equality once at the preceding coordinate (`'A+1T'`), and **substitutions** fail equality at the substituted coordinate:

$$\text{Concordance}[x] = \frac{\text{Matches}[x]}{\text{Depth}[x]}$$

---

## 3. Separating Structural Masking from Point & Motif Discordance

A critical requirement for motif discovery is keeping isolated point and motif discordances unmasked while excluding large structural artifacts (collapsed repeats, prophage dropouts, or unpolished structural blocks) so they do not distort background statistics.

Given the median non-zero contig depth $D_{\text{med}}$, the depth variance limit $V$ (default `--depth-variance 2.5`), and minimum concordance threshold $C_{\min}$ (default `--min-concordance 0.85`), `sbam` constructs `structural_mask` from two sources:

1. **Depth Outliers**: $\text{Depth}[x] > D_{\text{med}} \times V$ or $\text{Depth}[x] < D_{\text{med}} / V$.
2. **Contiguous Low-Concordance Blocks ($\ge 20\text{ bp}$)**: Contiguous runs of $\ge 20\text{ bp}$ where $\text{Concordance}[x] < C_{\min}$.

**Isolated discordant sites and short motif runs ($< 20\text{ bp}$) are intentionally NOT structurally masked**, ensuring they remain in the unmasked evaluated set ($E_{\text{eval}}$) for majority/minority error accounting and $k$-mer motif testing.

### Classification of Excluded Structural Regions ($\ge 50\text{ bp}$)

Contiguous structural blocks ($\ge 20\text{ bp}$) separated by $< 500\text{ bp}$ are merged, and regions $\ge 50\text{ bp}$ are classified by their mean depth and concordance:

* **`High Depth (...x of median; Collapsed Repeat / Multi-copy Operon)`**: Mean depth $> D_{\text{med}} \times V$.
* **`Low Depth (<10x; Coverage Gap / Dropout)`**: Mean depth $< 10.0\times$.
* **`Relative Depth Dip (...x of median; Structural Heterogeneity / Prophage / Mappability)`**: Mean depth is between $10.0\times$ and $D_{\text{med}} / V$ (distinguishing local population heterogeneity at e.g. $84\times$ depth from true $0\times$ gaps).
* **`Contiguous Low Concordance (Unpolished Block / Structural Variant)`**: Normal depth, but contiguous low read agreement. Any excluded block with mean concordance $< 90.0\%$ is highlighted in red; if total length exceeds $1,000\text{ bp}$ (or $\ge 2$ blocks), the assembly is flagged for `REVIEW`.

---

## 4. Majority Consensus Errors ($< 50\%$), Empirical $\text{QV}$ & Depth Guard

Across all unmasked evaluated bases ($E_{\text{eval}}$), discordant sites ($\text{Concordance}[x] < C_{\min}$, default $85\%$) are partitioned into two biologically distinct categories:

1. **Majority Consensus Disagreement ($\text{Concordance}[x] < 50\%$)**:
   Fewer than half of the aligned reads support the base written in the FASTA assembly ($E_{\text{maj}}$). These represent candidate polishing errors, homopolymer frameshifts, or uncollapsed repeat variants. `sbam` reports both the rate per megabase ($\text{Errors / Mb}$) and an **Empirical Read-Concordance Quality Value ($\text{QV}$)**:

   $$\text{QV}_{\text{emp}} = \begin{cases} -10 \log_{10}\left( \dfrac{E_{\text{maj}}}{E_{\text{eval}}} \right) & \text{if } E_{\text{maj}} > 0 \\[8pt] -10 \log_{10}\left( \dfrac{0.5}{E_{\text{eval}}} \right) \;\text{(lower bound)} & \text{if } E_{\text{maj}} = 0 \end{cases}$$

   *Note:* Because $\text{QV}_{\text{emp}}$ is reference-free and derived from read-to-assembly agreement, it is an optimistic upper bound when reads and the assembly share a systematic consensus error.

2. **Minority Sub-Consensus Read Disagreement ($50\% \le \text{Concordance}[x] < 85\%$)**:
   The assembly base is supported by the majority of reads, but a significant minority ($15\%\text{--}50\%$) disagree. In high-depth ONT assemblies, this reflects strand-specific basecaller jitter or DNA modification (methylation) signatures rather than wrong assembly bases, and therefore does not penalize the assembly verdict.

### Low-Depth Guard ($< 30\times$)

At sequencing depths below $30\times$, binomial sampling variance increases the chance of random $< 50\%$ fluctuations. When primary contig depth is between $5\times$ and $30\times$, `sbam` caps the verdict at `ACCEPT WITH CAUTION (LOW DEPTH <30x)` and notes that coverage is insufficient for high-confidence base-level grading.

---

## 5. Single-Counted $k$-mer Testing, 3-State Orientation & Motif Family Clustering

Across all unmasked regions of the primary chromosome, `BaseAccuracy.analyze_motifs()` scans sliding windows of length $k \in \{4, 5, 6\}$ (`--k-list`).

### Single-Counted Canonical & Palindromic Windows

Each genomic window $[i, i+k)$ is counted **once** ($n_m \mathrel{+}= 1$):

* **Non-palindromic canonical pairs (`FWD/REV`)**: If the reference window matches the lexicographically smaller forward orientation (`FWD`), concordance is tested at the center offset $i + \lfloor k / 2 \rfloor$. If the reference window matches the reverse-complement orientation (`REV`), concordance is tested at the exact base-paired homologous coordinate $i + k - 1 - \lfloor k / 2 \rfloor$.
* **Palindromic $k$-mers (`raw_kmer == revcomp(raw_kmer)`, e.g., `TCGA/TCGA` or `GTCGAC/GTCGAC`)**: Each genomic occurrence is counted once at $i + \lfloor k / 2 \rfloor$ and marked as `Palindromic`, preventing artificial doubling of $n_m$ and inflation of binomial $p$-values.

### Consistent Baseline, Bonferroni & Wilson $95\%$ Lower-CI Filter

Every candidate motif $m$ across all $k$ is tested against the single chromosome-wide baseline discordance rate:

$$p_0 = \frac{D_{\text{total}}}{E_{\text{eval}}}$$

To eliminate low-count noise (such as $4$ hits in $100$ occurrences), a motif must satisfy four criteria to be called significant:

1. Minimum genomic occurrences: $n_m \ge 100$.
2. Minimum discordant sites: $d_m \ge 10$.
3. One-sided binomial $p$-value passing a single global Bonferroni threshold across all $M_{\text{total}}$ tested motifs across all $k$:

   $$P(X \ge d_m \mid n_m, p_0) < \frac{0.01}{M_{\text{total}}}$$

4. **Wilson $95\%$ Confidence Interval Lower Bound on Fold Enrichment $\ge 2.0\times$**:
   Using the lower bound $p_{\text{LB}}$ of the $95\%$ Wilson score interval for $d_m / n_m$, the conservative fold enrichment must satisfy:

   $$\text{CI}_{\text{lower fold}} = \frac{p_{\text{LB}}}{p_0} \ge 2.0\times$$

### 3-State Reference-Orientation Classification

For each non-palindromic motif, `sbam` tests whether discordant hits ($d_+$, $d_-$) are distributed proportionally to the motif's forward and reverse-complement reference occurrences ($n_+$, $n_-$) using a two-sided binomial test ($p_{\text{orient}}$) with null probability $p_+ = n_+ / (n_+ + n_-)$:

* **`Orientation-Asymmetric`**: $p_{\text{orient}} < 0.01$ **and** $\max(d_+, d_-) / d_m \ge 80\%$ (strong orientation skew on the reference strand, e.g., $0/54$ or $36/0$).
* **`Orientation-Balanced`**: $d_m \ge 20$, $p_{\text{orient}} \ge 0.05$, **and** $\max(d_+, d_-) / d_m \le 65\%$ (well-powered and balanced across both reference orientations).
* **`Undetermined (Low Count / Intermediate)`**: $d_m < 20$ or intermediate split ($0.01 \le p_{\text{orient}} < 0.05$ or $65\% < \text{share} < 80\%$, e.g., $9/1$ or $8/4$), avoiding false claims of symmetry when statistical power is low.
* **`Palindromic`**: Motif is identical to its reverse complement.
* **`Consensus Error (<50% Support)`**: Flagged in red whenever the family contains one or more majority consensus error sites ($< 50\%$ read agreement).

### Mutually Exclusive Motif Family Clustering

Because a single biological context (e.g., `GTCGAC`) triggers overlapping $4$-, $5$-, and $6$-mers (`TCGA`, `GTCGA`, `CGTCGA`, `GTCGAC`), `sbam` greedily clusters motifs that share $\ge 25\%$ of their discordant genomic coordinates or share a nested substring/core into a single **Motif Family**.

To ensure the table rows add up cleanly without double-counting overlapping boundaries between families, each discordant coordinate is partitioned **mutually exclusively** to its highest-ranked family (`Assigned Sites (% Explained)`). The sum of `Assigned Sites` across all rows equals the exact deduplicated headline union of motif-explained discordant sites, while raw pre-deduplication family counts are shown in brackets whenever families overlap.