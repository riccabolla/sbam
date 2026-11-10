# Replication Structure, DnaA Corroboration & Rotation

Most bacterial chromosomes replicate bidirectionally from a single origin of replication (*oriC*) toward a diametrically opposite terminus (*ter*) via theta ($\theta$) replication. Because the leading and lagging strands experience different single-stranded exposure times during replication, cytosine deamination ($C \rightarrow T$) preferentially enriches the leading strand in Guanine ($G$) over Cytosine ($C$).

When a bacterial chromosome is correctly assembled, this mutational bias creates a smooth, V-shaped cumulative GC-skew profile whose global minimum marks *oriC* and whose global maximum marks *ter*, separated by $\approx 180^\circ$ around the circular chromosome. Large chimeric joins, misresolved rRNA/IS-element inversions, or uncollapsed tandem duplications distort this geometry.

---

## 1. Chromosomal Size Filter ($\ge 500\text{ kb}$)

By default, `BioVal.calculate_symmetry()` (`sbam.struct.ori`) only evaluates contigs with length $L \ge 500,000\text{ bp}$ (configurable via `--min-ori-length`).

Plasmids and megaplasmids ($< 500\text{ kb}$) frequently replicate via unidirectional theta or rolling-circle mechanisms and carry mosaic horizontal gene transfer islands (e.g., composite transposons and AMR cassettes) that do not follow bidirectional chromosomal replichore symmetry. Excluding plasmids $< 500\text{ kb}$ prevents intact plasmids from triggering false replication warnings.

---

## 2. Zero-Centered & Circularly Smoothed Cumulative GC Skew

Each chromosomal contig is partitioned into non-overlapping sliding windows of size $W = 1,000\text{ bp}$. For each window $i$, the local GC skew $s_i$ is calculated as:

$$s_i = \begin{cases} \dfrac{G_i - C_i}{G_i + C_i} & \text{if } G_i + C_i > 0 \\ 0 & \text{if } G_i + C_i = 0 \end{cases}$$

### Detrending Global Strand Bias & Circular Smoothing

In raw cumulative GC skew ($\sum s_i$), any slight genome-wide imbalance between total $G$ and total $C$ on the forward FASTA strand introduces a linear ramp across the contig. That linear drift causes the global minimum and maximum to shift depending on where the circular chromosome was arbitrarily linearized.

To make *oriC* and *ter* detection mathematically invariant to circular rotation, `sbam` zero-centers the window skews by subtracting the mean window skew $\bar{s}$ prior to cumulative summation:

$$\tilde{s}_i = s_i - \frac{1}{N}\sum_{j=1}^{N} s_j$$

$$R_k = \sum_{i=1}^{k} \tilde{s}_i \quad \text{for } k \in \{1, \dots, N\}$$

To prevent small mobile elements or local genomic islands from jittering the extrema, $R_k$ is smoothed across the circular boundary using a circular moving-average kernel of width $K = \min(15, \lfloor N / 10 \rfloor)$ windows ($15\text{ kb}$):

$$C_k = \frac{1}{K} \sum_{j = k - \lfloor K/2 \rfloor}^{k + \lfloor K/2 \rfloor} R_{(j \bmod N)}$$

The predicted origin and terminus coordinates are then identified from the global extrema of the smoothed curve $C_k$:

$$\text{oriC}_{\text{pos}} = W \times \arg\min_{k} (C_k), \qquad \text{ter}_{\text{pos}} = W \times \arg\max_{k} (C_k)$$

---

## 3. Replichore Symmetry & Plausibility Thresholds

The shortest physical distance $d$ between *oriC* and *ter* along the circular chromosome of length $L$ is:

$$d = \min\left( |\text{oriC}_{\text{pos}} - \text{ter}_{\text{pos}}|,\; L - |\text{oriC}_{\text{pos}} - \text{ter}_{\text{pos}}| \right)$$

Converting this physical distance into angular degrees (where $L / 2$ corresponds to a perfect $180^\circ$ opposition) yields the **Replichore Symmetry** score:

$$\text{Symmetry} = \left( \frac{d}{L / 2} \right) \times 180^\circ$$

`sbam` classifies the biological plausibility of the chromosome into three tiers:

| Symmetry Range | Viability Status | Badge Color | Biological Interpretation |
| :--- | :--- | :--- | :--- |
| $165.0^\circ \le \text{Symmetry} \le 180.0^\circ$ | `ACCEPTABLE` | Green | Balanced replichores consistent with an intact, bidirectionally replicating chromosome. |
| $150.0^\circ \le \text{Symmetry} < 165.0^\circ$ | `ATYPICAL` | Amber / Orange | Moderate replichore asymmetry (triggers `REVIEW`); may reflect recent large phage integration, natural rearrangements, or a minor inversion. |
| $\text{Symmetry} < 150.0^\circ$ | `WARNING` | Red | Severe geometric distortion (triggers `REVIEW`); strongly indicative of a major misassembly, chimeric contig join, or large inversion across replichores. |

---

## 4. Genome-Wide DnaA-Box Baseline & Poisson Enrichment

To provide an independent sequence-motif corroboration of the GC-skew minimum, `sbam` scans the entire chromosome for occurrences of the 9-bp bacterial DnaA box consensus (`TTATCCACA` on the forward strand or `TGTGGATAA` on the reverse strand, allowing $\le 1$ mismatch).

Let $N_{\text{total}}$ be the total number of DnaA boxes across the chromosome of length $L$. Under a uniform genomic background, the expected number of DnaA boxes inside the $50\text{ kb}$ window ($\pm 25\text{ kb}$) centered on $\text{oriC}_{\text{pos}}$ is:

$$\lambda_{\text{exp}} = N_{\text{total}} \times \frac{50,000}{L}$$

Given the observed count $k_{\text{obs}}$ within $\pm 25\text{ kb}$ of $\text{oriC}_{\text{pos}}$, `sbam` reports the fold enrichment ($k_{\text{obs}} / \lambda_{\text{exp}}$) and the one-sided Poisson survival $p$-value:

$$P(X \ge k_{\text{obs}}) = 1 - \sum_{j=0}^{k_{\text{obs}} - 1} \frac{\lambda_{\text{exp}}^j e^{-\lambda_{\text{exp}}}}{j!}$$

---

## 5. *oriC* Offset to Coordinate $0\text{ bp}$ (3-Tier Rotation Check)

Standard long-read post-processing pipelines (such as `Autocycler` or `dnaapler`) rotate circularized bacterial chromosomes so that `dnaA` (located at or near *oriC*) begins at coordinate $0\text{ bp}$. Because public genes (`dnaA`, `repA`) and the physical GC-skew minimum can be separated by $10\text{--}150\text{ kb}$ depending on local gene arrangement around the origin, `sbam` computes the shortest circular distance from $\text{oriC}_{\text{pos}}$ to coordinate $0\text{ bp}$:

$$\text{Offset}_{\text{bp}} = \min\left( \text{oriC}_{\text{pos}},\; L - \text{oriC}_{\text{pos}} \right), \qquad \text{Offset}_{\text{deg}} = \left( \frac{\text{Offset}_{\text{bp}}}{L / 2} \right) \times 180^\circ$$

and assigns an informational 3-tier **Rotation Status**:

* **`ROTATED`** ($\text{Offset} \le 25.0\text{ kb}$): The GC-skew minimum sits directly at the $0\text{ bp}$ boundary.
* **`NEAR_ORIC`** ($25.0\text{ kb} < \text{Offset} \le 200.0\text{ kb}$ and $\text{Offset}_{\text{deg}} \le 20.0^\circ$): *oriC* is located near $0\text{ bp}$ (typical when rotated to a neighboring `dnaA` start codon).
* **`UNROTATED`** ($\text{Offset} > 200.0\text{ kb}$ or $\text{Offset}_{\text{deg}} > 20.0^\circ$): The chromosome was linearized far from *oriC*. `sbam` displays an informational notice recommending reorientation with `dnaapler` followed by re-verification.