# Replication Structure, DnaA Corroboration & Rotation

Most bacterial chromosomes replicate bidirectionally from a single origin of replication (*oriC*) toward a diametrically opposite terminus (*ter*) via theta ($\theta$) replication. Because the leading and lagging strands experience different single-stranded exposure times during replication, cytosine deamination ($C \rightarrow T$) preferentially enriches the leading strand in Guanine ($G$) over Cytosine ($C$).

When a bacterial chromosome is correctly assembled, this mutational bias creates a smooth, V-shaped cumulative GC-skew profile whose global minimum marks *oriC* and whose global maximum marks *ter*, separated by $\approx 180^\circ$ around a circular chromosome (or centered on a linear chromosome with *oriC* near $L / 2$). Large chimeric joins, misresolved rRNA/IS-element inversions, or uncollapsed tandem duplications distort this geometry.

---

## 1. Chromosomal Size Filter ($\ge 500\text{ kb}$ Default)

By default, `BioVal.calculate_symmetry()` (`sbam.struct.ori`) evaluates contigs with length $L \ge 500,000\text{ bp}$ (configurable via `--min-ori-length`; the chosen threshold is propagated dynamically to the verdict engine and HTML report).

Plasmids and megaplasmids ($< 500\text{ kb}$) frequently replicate via unidirectional theta or rolling-circle mechanisms and carry mosaic horizontal gene transfer islands (e.g., composite transposons and AMR cassettes) that do not follow bidirectional chromosomal replichore symmetry. Excluding plasmids below `--min-ori-length` prevents intact plasmids from triggering false replication warnings.

---

## 2. Length-Weighted, Zero-Centered & Smoothed Cumulative GC Skew

Each chromosomal contig is partitioned into non-overlapping windows of nominal size $W = 1,000\text{ bp}$ (with the trailing remainder window of length $L_N \le W$ assigned fractional weight $w_N = L_N / W$). For each window $i$, the local GC skew $s_i$ and weight $w_i$ are calculated as:

$$s_i = \begin{cases} \dfrac{G_i - C_i}{G_i + C_i} & \text{if } G_i + C_i > 0 \\ 0 & \text{if } G_i + C_i = 0 \end{cases}, \qquad w_i = \frac{L_i}{W}$$

### Detrending Global Strand Bias & Exact Boundary Coordinates

In raw cumulative GC skew ($\sum s_i$), any slight genome-wide imbalance between total $G$ and total $C$ on the forward FASTA strand introduces a linear ramp across the contig, and giving a tiny trailing remainder window equal weight distorts the end of the profile.

To make *oriC* and *ter* detection rotation-invariant and accurately aligned to $0\text{ bp}$, `sbam` subtracts the length-weighted mean skew $\bar{s}_w$ prior to cumulative summation and evaluates $R_k$ at exact window-boundary coordinates $x_k = k \cdot W$ starting with $R_0 = 0$:

$$\bar{s}_w = \frac{\sum_{j=0}^{N-1} w_j s_j}{\sum_{j=0}^{N-1} w_j}, \qquad \Delta_i = w_i (s_i - \bar{s}_w)$$

$$R_0 = 0, \qquad R_k = \sum_{i=0}^{k-1} \Delta_i$$

* **Circular Mode (default)**: $R$ has length $N$ ($k \in \{0, \dots, N-1\}$) and is smoothed with a circular moving-average kernel of width $K = \min(15, \lfloor N / 10 \rfloor)$ windows ($15\text{ kb}$), wrapping around the $0/L$ seam.
* **Linear Mode (`--linear`)**: $R$ includes the terminal $3'$ endpoint ($x_N = L$) and is smoothed with edge padding (`np.pad(..., mode="edge")`) rather than circular wrap-around.

The predicted origin and terminus coordinates are identified from the global extrema of the smoothed curve $C$:

$$\text{oriC}_{\text{pos}} = x_{\arg\min_k(C_k)}, \qquad \text{ter}_{\text{pos}} = x_{\arg\max_k(C_k)}$$

---

## 3. Replichore Symmetry & Plausibility Thresholds

* **Circular Mode (default)**: The shortest circular distance $d$ between *oriC* and *ter* is:

  $$d = \min\left( |\text{oriC}_{\text{pos}} - \text{ter}_{\text{pos}}|,\; L - |\text{oriC}_{\text{pos}} - \text{ter}_{\text{pos}}| \right)$$

* **Linear Mode (`--linear`)**: Linear chromosomes with bidirectional replication from a central origin place *oriC* near the midpoint $L / 2$ so both linear replichores have equal length. Symmetry is measured from *oriC* to the nearest chromosome end:

  $$d = \min\left(\text{oriC}_{\text{pos}},\; L - \text{oriC}_{\text{pos}}\right)$$

Converting $d$ into angular degrees (where $L / 2$ corresponds to $180^\circ$) yields the **Replichore Symmetry** score:

$$\text{Symmetry} = \left( \frac{d}{L / 2} \right) \times 180^\circ$$

`sbam` classifies the biological plausibility of the chromosome into three tiers:

| Symmetry Range | Viability Status | Badge Color | Biological Interpretation |
| :--- | :--- | :--- | :--- |
| $165.0^\circ \le \text{Symmetry} \le 180.0^\circ$ | `ACCEPTABLE` | Green | Balanced replichores consistent with an intact, bidirectionally replicating chromosome. |
| $150.0^\circ \le \text{Symmetry} < 165.0^\circ$ | `ATYPICAL` | Amber / Orange | Moderate replichore asymmetry (triggers `REVIEW`); may reflect recent large phage integration, natural rearrangements, or a minor inversion. |
| $\text{Symmetry} < 150.0^\circ$ | `WARNING` | Red | Severe geometric distortion (triggers `REVIEW`); strongly indicative of a major misassembly, chimeric contig join, or large inversion across replichores. |

---

## 4. Vectorized DnaA-Box Baseline & Poisson Enrichment

To provide an independent sequence-motif corroboration of the GC-skew minimum, `sbam` scans the entire chromosome using a precomputed 256 KB 2-bit lookup table (`_DNAA_LUT`) for occurrences of the 9-bp bacterial DnaA box consensus (`TTATCCACA` on the forward strand or `TGTGGATAA` on the reverse strand, allowing $\le 1$ mismatch). In circular mode, the scan appends `seq[:8]` so DnaA boxes spanning the $0/L$ seam are captured accurately.

Let $N_{\text{total}}$ be the total number of DnaA boxes across the chromosome of length $L$, and let $W_{\text{ori}}$ be the size of the $\pm 25\text{ kb}$ *oriC* window ($50,000\text{ bp}$ in circular mode, or $\min(L, \text{oriC}_{\text{pos}} + 25\text{ kb}) - \max(0, \text{oriC}_{\text{pos}} - 25\text{ kb})$ in linear mode). Under a uniform genomic background, the expected count inside the *oriC* window is:

$$\lambda_{\text{exp}} = N_{\text{total}} \times \frac{W_{\text{ori}}}{L}$$

Given the observed count $k_{\text{obs}}$ within $\pm 25\text{ kb}$ of $\text{oriC}_{\text{pos}}$, `sbam` reports the fold enrichment ($k_{\text{obs}} / \lambda_{\text{exp}}$) and the one-sided Poisson survival $p$-value:

$$P(X \ge k_{\text{obs}}) = 1 - \sum_{j=0}^{k_{\text{obs}} - 1} \frac{\lambda_{\text{exp}}^j e^{-\lambda_{\text{exp}}}}{j!}$$

---

## 5. *oriC* Offset & Rotation Check

Standard long-read post-processing pipelines (such as `Autocycler` or `dnaapler`) rotate circularized bacterial chromosomes so that `dnaA` (located at or near *oriC*) begins at coordinate $0\text{ bp}$. In circular mode, `sbam` computes the shortest circular distance from $\text{oriC}_{\text{pos}}$ to coordinate $0\text{ bp}$:

$$\text{Offset}_{\text{bp}} = \min\left( \text{oriC}_{\text{pos}},\; L - \text{oriC}_{\text{pos}} \right), \qquad \text{Offset}_{\text{deg}} = \left( \frac{\text{Offset}_{\text{bp}}}{L / 2} \right) \times 180^\circ$$

and assigns an informational **Rotation Status**:

* **`ROTATED`** ($\text{Offset} \le 25.0\text{ kb}$): The GC-skew minimum sits directly at the $0\text{ bp}$ boundary.
* **`NEAR_ORIC`** ($25.0\text{ kb} < \text{Offset} \le 200.0\text{ kb}$ and $\text{Offset}_{\text{deg}} \le 20.0^\circ$): *oriC* is located near $0\text{ bp}$ (typical when rotated to a neighboring `dnaA` start codon).
* **`UNROTATED`** ($\text{Offset} > 200.0\text{ kb}$ or $\text{Offset}_{\text{deg}} > 20.0^\circ$): The chromosome was linearized far from *oriC*. `sbam` displays an informational notice recommending reorientation with `dnaapler` followed by re-verification.
* **`LINEAR`** (`--linear` mode): Rotation to coordinate $0\text{ bp}$ is not applicable to naturally linear chromosomes; `oriC` offset records distance to the nearest contig end and `dnaapler` rotation notices are suppressed.
