# Base-Level Concordance, Error Stratification & Motif Family Analysis

Even when a long-read assembly is structurally circular and geometrically symmetric, basecaller models and consensus polishers can leave consensus-level errors (such as homopolymer indels or uncollapsed repeat variants) or exhibit sub-consensus read disagreement at specific sequence contexts (such as methylated motifs or G/C-rich homopolymers).

`sbam` bypasses `FASTQ` quality scores by profiling empirical read-to-assembly concordance across **all contigs** in the assembly, stratifying substitutions vs. indels and sequencing read-strand bias, and performing a vectorized motif screen on the primary chromosome.

---

## 1. All-Contig CIGAR Sweeps & Cyclic Buffer Depth Folding (`sbam.base.base`)

In circular mode, reads originating from the first $B\text{ bp}$ of a circular contig ($[0, B)$) map with equal probability to $[0, B)$ or to the identical cyclic tail $[L, L+B)$. To prevent an artificial $50\%$ coverage drop across $[0, B)$, `BaseAccuracy` profiles the full buffered length $L + B$ (`buffer_len = 0` when `--linear` is active or when a pre-aligned `--bam` lacks `_cyclic` suffixes) and folds the cyclic tail back onto $[0, B)$ across all per-base arrays:

$$\text{Depth}[0 : B] \leftarrow \text{Depth}_{\text{raw}}[0 : B] + \text{Depth}_{\text{raw}}[L : L + B]$$

$$\text{Matches}[0 : B] \leftarrow \text{Matches}_{\text{raw}}[0 : B] + \text{Matches}_{\text{raw}}[L : L + B]$$

### High-Throughput Read-Centric `--eqx` CIGAR Sweep

Rather than running slow column-by-column `pysam.pileup()` loops, `ReadAligner` invokes `minimap2` with `--eqx --secondary=no`, and `_profile_chunk()` walks each primary read's `cigartuples` over flat $\sim 250\text{ kb}$ chunks using $O(1)$-per-run NumPy difference arrays (`np.cumsum`), bypassing `multiprocessing.Pool` when `--threads 1`:

* **Excluded Alignments**: Unmapped, secondary, supplementary (`0x800`), QC-fail (`0x200`), and duplicate (`0x400`) alignments are strictly excluded so split/chimeric supplementary alignments never double-count cyclic seam bases.
* **`=` (CIGAR op `7`) & `X` (CIGAR op `8`)**: Directly accumulate per-base depth, match, and mismatch counts in $O(1)$ per CIGAR run.
* **`M` (CIGAR op `0` Fallback)**: For pre-aligned external `--bam` files generated without `--eqx`, `_profile_chunk()` automatically falls back to vectorized `np.frombuffer` sequence equality against the reference chunk.

Every contig in the FASTA (primary chromosome and all secondary plasmids/replicons) is profiled, and its majority consensus error count ($< 50\%$ read support), error rate per $\text{Mb}$, and structurally masked percentage are reported in the **Assembly Overview** table.

---

## 2. Indel-Aware Concordance, Error Stratification & Read-Strand Bias

In Oxford Nanopore assemblies, systematic basecalling and polishing errors predominantly manifest as $1\text{--}2\text{ bp}$ insertions or deletions rather than pure substitutions.

### Indel Accounting Without Deletion Double-Counting

* **Deletions (`D`, CIGAR op `2`)**: Increment `depth` and `deletions` once at each deleted reference coordinate without adding a match.
* **Insertions (`I`, CIGAR op `1`)**: Increment `insertions` at the immediately preceding reference anchor (`rpos - 1`) and—if that anchor was counted as a base match—subtract $1$ from `matches` once at that anchor coordinate.
* **Substitutions (`X` or mismatching `M`)**: Increment `depth` and `mismatches` at the substituted reference coordinate.

$$\text{Concordance}[x] = \frac{\text{Matches}[x]}{\text{Depth}[x]}$$

### SNP vs. Indel & Sequencing Read-Strand Stratification

During the same CIGAR sweep, `_profile_chunk()` tracks forward (`not read.is_reverse`) vs. reverse (`read.is_reverse`) read depth and matches separately (`fwd_depth`, `fwd_match`, `rev_depth`, `rev_match`):

* **Genome-Wide Error Spectrum**: Reports unmasked read-base `snp_rate` and `indel_rate`, plus counts of `< 50%` consensus error sites dominated by SNPs vs. indels (`DEL` / `INS` / `MIXED`).
* **Read-Strand Bias Fisher Test**: At every `< 50%` consensus error site and for every clustered motif family, `sbam` runs a two-sided Fisher's exact test (`scipy.stats.fisher_exact`) comparing forward vs. reverse sequencing read concordance (`p < 0.01` with $\ge 15\%$ strand concordance delta flagged as read-strand biased), separating methylation-induced ONT strand-specific basecaller bias from genuine assembly errors.

---

## 3. Separating Structural Masking from Point & Motif Discordance

A critical requirement for motif discovery is keeping isolated point and motif discordances unmasked while excluding large structural artifacts (collapsed repeats, prophage dropouts, or unpolished structural blocks) so they do not distort background statistics.

Given the median non-zero contig depth $D_{\text{med}}$, the depth variance limit $V$ (default `--depth-variance 2.5`, required $> 1.0$), and minimum concordance threshold $C_{\min}$ (default `--min-concordance 0.85`, required $> 0.50$), `sbam` constructs `structural_mask` from two sources:

1. **Depth Outliers**: $\text{Depth}[x] > D_{\text{med}} \times V$ or $\text{Depth}[x] < D_{\text{med}} / V$.
2. **Contiguous Low-Concordance Blocks ($\ge 20\text{ bp}$)**: Contiguous runs of $\ge 20\text{ bp}$ where $\text{Concordance}[x] < C_{\min}$.

**Isolated discordant sites and short motif runs ($< 20\text{ bp}$) are intentionally NOT structurally masked**, ensuring they remain in the unmasked evaluated set ($E_{\text{eval}}$) for majority/minority error accounting and $k$-mer motif testing.

### Category-Aware Merging & Reporting of Excluded Structural Regions ($\ge 20\text{ bp}$)

Contiguous structural blocks ($\ge 20\text{ bp}$) are assigned a category (`HIGH_DEPTH`, `LOW_DEPTH`, or `LOW_CONCORDANCE`). Only adjacent blocks of the **same category** separated by $< 500\text{ bp}$ are merged, and `mean_depth` / `mean_concordance` are computed strictly over the masked bases inside `[start, end)` so intervening normal sequence never dilutes $< 90\%$ unpolished blocks:

* **`High Depth (...x of median; Collapsed Repeat / Multi-copy Operon)`**: Mean depth $> D_{\text{med}} \times V$.
* **`Low Depth (<10x; Coverage Gap / Dropout)`**: Mean depth $< 10.0\times$.
* **`Relative Depth Dip (...x of median; Structural Heterogeneity / Prophage / Mappability)`**: Mean depth is between $10.0\times$ and $D_{\text{med}} / V$.
* **`Contiguous Low Concordance (Unpolished Block / Structural Variant)`**: Normal depth, but contiguous low read agreement.

All `< 90.0%` concordance blocks (`low_concordance_blocks`) are retained without truncation in the verdict engine, JSON metrics, and `sbam summarize`—alongside the top 50 longest remaining regions in the HTML report and **all** $\ge 20\text{ bp}$ intervals across all contigs in `<prefix>_masked_regions.bed`. If total `< 90%` excluded blocks exceed $1,000\text{ bp}$ (or $\ge 2$ blocks) or total masked bases exceed $5.0\%$ on either the primary chromosome or any secondary contig, the assembly is flagged for `REVIEW`.

---

## 4. Majority Consensus Errors ($< 50\%$), Empirical $\text{QV}$ & Depth Guard

Across all unmasked evaluated bases ($E_{\text{eval}}$), discordant sites ($\text{Concordance}[x] < C_{\min}$, default $85\%$) are partitioned into two biologically distinct categories:

1. **Majority Consensus Disagreement ($\text{Concordance}[x] < 50\%$)**:
   Fewer than half of the aligned reads support the base written in the FASTA assembly ($E_{\text{maj}}$). These represent candidate polishing errors, homopolymer frameshifts, or uncollapsed repeat variants. Every `< 50%` site across all contigs is exported to `<prefix>_consensus_errors.tsv` and `<prefix>_consensus_errors.vcf`. `sbam` reports both the rate per megabase ($\text{Errors / Mb}$) and an **Empirical Read-Concordance Quality Value ($\text{QV}$)**:

   $$\text{QV}_{\text{emp}} = \begin{cases} -10 \log_{10}\left( \dfrac{E_{\text{maj}}}{E_{\text{eval}}} \right) & \text{if } E_{\text{maj}} > 0 \\[8pt] -10 \log_{10}\left( \dfrac{0.5}{E_{\text{eval}}} \right) \;\text{(lower bound)} & \text{if } E_{\text{maj}} = 0 \end{cases}$$

   *Note:* Because $\text{QV}_{\text{emp}}$ is reference-free and derived from read-to-assembly agreement, it is an optimistic upper bound when reads and the assembly share a systematic consensus error.

2. **Minority Sub-Consensus Read Disagreement ($50\% \le \text{Concordance}[x] < 85\%$)**:
   The assembly base is supported by the majority of reads, but a significant minority ($15\%\text{--}50\%$) disagree. In high-depth ONT assemblies, this reflects strand-specific basecaller jitter or DNA modification (methylation) signatures rather than wrong assembly bases, and therefore does not penalize the assembly verdict.

### Low-Depth Guard ($< 30\times$)

At sequencing depths below $30\times$, binomial sampling variance increases the chance of random $< 50\%$ fluctuations. When primary contig depth is between $5\times$ and $30\times$, `sbam` caps an otherwise passing verdict at `ACCEPT WITH CAUTION (LOW DEPTH <30x)` and notes that coverage is insufficient for high-confidence base-level grading (while still escalating to `REVIEW` or `REVIEW (SIZE ANOMALY)` if secondary contigs fail, masking exceeds $5\%$, or genome size is anomalous).

---

## 5. Vectorized $k$-mer Testing, Dual-Center Palindromes & Motif Family Clustering

Across all unmasked regions of the primary chromosome, `BaseAccuracy.analyze_motifs()` encodes the sequence into 2-bit rolling integer hashes, counts genome-wide motif occurrences in C speed via precomputed $4^k$ lookup tables (`k <= 7`, sparse `np.unique` for `k >= 8`) and `np.bincount()`, and evaluates sparse discordant coordinates across $k \in \{4, 5, 6\}$ (`--k-list`).

### Single-Counted Canonical & Dual-Center Palindromic Windows

Each valid unmasked genomic window $[i, i+k)$ is counted **once** ($n_m \mathrel{+}= 1$):

* **Non-palindromic canonical pairs (`FWD/REV`)**: If the reference window matches the lexicographically smaller forward orientation (`FWD`, or in non-canonical mode `canonical=False`), concordance is tested at $i + \lfloor k / 2 \rfloor$. If it matches the reverse-complement orientation (`REV`), concordance is tested at the exact base-paired homologous coordinate $i + k - 1 - \lfloor k / 2 \rfloor$.
* **Palindromic $k$-mers (`raw_kmer == revcomp(raw_kmer)`, e.g., `GATC/GATC`, `TCGA/TCGA`, `GTCGAC/GTCGAC`)**: Each genomic occurrence is counted **once** ($n_m \mathrel{+}= 1$) and marked `Palindromic`. For **even-length palindromes** ($k \bmod 2 = 0$), both symmetric central bases ($i + k/2 - 1$ and $i + k/2$, such as both `A` and `T` in `GATC`) are tested per occurrence (so 6mA modifications at position 1 of `GATC` are detected just like position 2), using the two-base null discordance probability $p_{\text{null}} = 1 - (1 - p_0)^2$ in the binomial test.

### Consistent Baseline, Bonferroni & Wilson $95\%$ Lower-CI Filter

Every candidate motif $m$ across all $k$ is tested against the chromosome-wide baseline discordance rate $p_0 = D_{\text{total}} / E_{\text{eval}}$ (or $1 - (1 - p_0)^2$ for even-length palindromes). A motif must satisfy four criteria to be called significant:

1. Minimum genomic occurrences: $n_m \ge 100$.
2. Minimum discordant sites: $d_m \ge 10$.
3. One-sided binomial $p$-value passing a global Bonferroni threshold across all $M_{\text{total}}$ tested motifs across all $k$:

   $$P(X \ge d_m \mid n_m, p_{\text{null}}) < \frac{0.01}{M_{\text{total}}}$$

4. **Wilson $95\%$ Confidence Interval Lower Bound on Fold Enrichment $\ge 2.0\times$**:
   Using the lower bound $p_{\text{LB}}$ of the $95\%$ Wilson score interval for $d_m / n_m$, the conservative fold enrichment must satisfy:

   $$\text{CI}_{\text{lower fold}} = \frac{p_{\text{LB}}}{p_0} \ge 2.0\times$$

### 3-State Reference-Orientation Classification

For each non-palindromic motif, `sbam` tests whether discordant hits ($d_+$, $d_-$) are distributed proportionally to the motif's forward and reverse-complement reference occurrences ($n_+$, $n_-$) using a two-sided binomial test ($p_{\text{orient}}$) with null probability $p_+ = n_+ / (n_+ + n_-)$:

* **`Orientation-Asymmetric`**: $p_{\text{orient}} < 0.01$ **and** $\max(d_+, d_-) / d_m \ge 80\%$ (strong orientation skew on the reference strand).
* **`Orientation-Balanced`**: $d_m \ge 20$, $p_{\text{orient}} \ge 0.05$, **and** $\max(d_+, d_-) / d_m \le 65\%$ (well-powered and balanced across both reference orientations).
* **`Undetermined (Low Count / Intermediate)`**: $d_m < 20$ or intermediate split ($0.01 \le p_{\text{orient}} < 0.05$ or $65\% < \text{share} < 80\%$).
* **`Palindromic`**: Motif is identical to its reverse complement.
* **`Consensus Error (<50% Support)`**: Flagged in red whenever the family contains one or more exclusively assigned majority consensus error sites ($< 50\%$ read agreement).

### Mutually Exclusive Motif Family Clustering

Because a single biological context (e.g., `GTCGAC`) triggers overlapping $4$-, $5$-, and $6$-mers (`TCGA`, `GTCGA`, `CGTCGA`, `GTCGAC`), `sbam` greedily clusters motifs that share $\ge 25\%$ of their discordant genomic coordinates or share a nested substring or $\ge 5\text{ bp}$ shifted/trimmed core (`len >= 6`, avoiding over-merging of distinct 5-mers with $0\%$ overlap) into a single **Motif Family**.

To guarantee exact reconciliation between the Motif Families table, headline metrics, HTML dashboard, and `sbam summarize`, **both** discordant coordinates (`assigned_unique_sites`) and $< 50\%$ consensus error coordinates (`consensus_error_count`) are partitioned **mutually exclusively** to the highest-ranked family—preventing shared consensus errors from inflating `consensus_error_bases` or corrupting $\text{QV}$.
