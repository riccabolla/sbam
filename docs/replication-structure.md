# Replication Structure & Geometry

Most bacterial chromosomes replicate bidirectionally from a single origin of replication (*oriC*) toward a diametrically opposite terminus (*ter*) via theta ($\theta$) replication. Because the leading and lagging strands experience different single-stranded exposure times during replication, cytosine deamination ($C \rightarrow T$) preferentially enriches the leading strand in Guanine ($G$) over Cytosine ($C$).

When a bacterial chromosome is correctly assembled, this mutational bias creates a smooth, V-shaped cumulative GC-skew profile whose global minimum marks *oriC* and whose global maximum marks *ter*, separated by $\approx 180^\circ$ around the circular chromosome. Large chimeric joins, misresolved rRNA/IS-element inversions, or uncollapsed tandem duplications distort this geometry.

---

## 1.Chromosomal Size Filter ($\ge 500\text{ kb}$)

By default, `BioVal.calculate_symmetry()` (`sbam.struct.ori`) only evaluates contigs with length $L \ge 500,000\text{ bp}$ (configurable via `--min-ori-length`). 

Plasmids and megaplasmids ($< 500\text{ kb}$) frequently replicate via unidirectional theta or rolling-circle mechanisms and carry mosaic horizontal gene transfer islands (e.g., composite transposons and AMR cassettes) that do not follow bidirectional chromosomal replichore symmetry. Excluding plasmids $< 500\text{ kb}$ prevents intact plasmids from triggering false replication warnings.

---

## 2. Zero-Centered Cumulative GC Skew

Each chromosomal contig is partitioned into non-overlapping sliding windows of size $W = 1,000\text{ bp}$. For each window $i$, the local GC skew $s_i$ is calculated as:

$$s_i = \begin{cases} \dfrac{G_i - C_i}{G_i + C_i} & \text{if } G_i + C_i > 0 \\ 0 & \text{if } G_i + C_i = 0 \end{cases}$$

### Detrending Global Strand Bias for Rotation Invariance

In raw cumulative GC skew ($\sum s_i$), any slight genome-wide imbalance between total $G$ and total $C$ on the forward FASTA strand introduces a linear ramp across the contig. That linear drift causes the global minimum and maximum to shift depending on where the circular chromosome was arbitrarily linearized.

To make *oriC* and *ter* detection mathematically invariant to circular rotation, `sbam` zero-centers the window skews by subtracting the mean window skew $\bar{s}$ prior to cumulative summation:

$$\tilde{s}_i = s_i - \frac{1}{N}\sum_{j=1}^{N} s_j$$

$$C_k = \sum_{i=1}^{k} \tilde{s}_i \quad \text{for } k \in \{1, \dots, N\}$$

Because $\sum_{i=1}^{N} \tilde{s}_i = 0$, the cumulative curve $C_k$ closes smoothly across the $0\text{ bp} / L\text{ bp}$ circular boundary. The predicted origin and terminus coordinates are then identified from the global extrema:

$$\text{oriC}_{\text{pos}} = W \times \arg\min_{k} (C_k), \qquad \text{ter}_{\text{pos}} = W \times \arg\max_{k} (C_k)$$

---

## 3. Replichore Symmetry & Plausibility Thresholds

The shortest physical distance $d$ between *oriC* and *ter* along the circular chromosome of length $L$ is:

$$d = \min\left( |\text{oriC}_{\text{pos}} - \text{ter}_{\text{pos}}|,\; L - |\text{oriC}_{\text{pos}} - \text{ter}_{\text{pos}}| \right)$$

Converting this physical distance into angular degrees (where $L / 2$ corresponds to a perfect $180^\circ$ opposition) yields the **Replichore Symmetry** score:

$$\text{Symmetry} = \left( \frac{d}{L / 2} \right) \times 180^\circ$$

`sbam` classifies the biological plausibility of the chromosome into three tiers:

| Symmetry Range | Viability Status | Biological Interpretation |
| :--- | :--- | :--- |
| $165.0^\circ \le \text{Symmetry} \le 180.0^\circ$ | `ACCEPTABLE` | Balanced replichores consistent with an intact, bidirectionally replicating chromosome. |
| $150.0^\circ \le \text{Symmetry} < 165.0^\circ$ | `ATYPICAL` | Moderate replichore asymmetry; may reflect recent large phage integration, natural genomic rearrangements, or a minor inversion. |
| $\text{Symmetry} < 150.0^\circ$ | `WARNING` | Severe geometric distortion; strongly indicative of a major misassembly, chimeric contig join, or large inversion across replichores. |

---

## 4. *oriC* Offset to Coordinate $0\text{ bp}$ (Rotation Check)

Standard long-read post-processing pipelines (such as `Autocycler` or `dnaapler`) rotate circularized bacterial chromosomes so that `dnaA` (located at or immediately adjacent to *oriC*) begins at coordinate $0\text{ bp}$.

To verify whether an assembly has been properly rotated, `sbam` calculates the shortest circular distance from $\text{oriC}_{\text{pos}}$ to coordinate $0\text{ bp}$:

$$\text{Offset}_{\text{bp}} = \min\left( \text{oriC}_{\text{pos}},\; L - \text{oriC}_{\text{pos}} \right)$$

$$\text{Offset}_{\text{deg}} = \left( \frac{\text{Offset}_{\text{bp}}}{L / 2} \right) \times 180^\circ$$

An offset near $0\text{--}150\text{ kb}$ ($\le 10^\circ$) confirms that the chromosome is standardized near *oriC*, whereas a large offset (e.g., $> 1,000\text{ kb}$) indicates the chromosome was left at an arbitrary assembly graph breakpoint.