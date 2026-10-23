# Validation

What has been checked, how, and with which result. Every number below can be
regenerated: the real-data numbers with `Rscript examples/validate_example_A.R`
(and `Rscript examples/run_example_A.R`), the synthetic ones with the test
suite (`Rscript run_tests.R`).

Data: **dataset A** (`examples/dataset_A`) — one filter measured with all four
instruments: FTIR PerkinElmer Spotlight (33 particles), FTIR Bruker LUMOS (42),
WITec Raman (43), Agilent 8700 LDIR (40). Almost all particles are PET, so
dataset A tests geometry well and material agreement only weakly. Plus
synthetic data with known answers, and the Raman file of a field sample
(Faanui SW1L R4 B, 14,562 particles) for the Raman image placement.

## 1. Alignment

Every instrument is aligned into the Raman frame by a similarity transform
(rotation, uniform scale, translation, scale within 0.8-1.25): landmarks
(large particles and fibres), else RANSAC / global registration on the
quality-passing plastics, then ICP. (Without Raman the frame
is FTIR (PerkinElmer), else FTIR (Bruker); the numbers below are all for the
Raman frame.)

**Dataset A**

| Pair | rotation | scale | tier used |
|---|---|---|---|
| FTIR (PerkinElmer) → Raman | −178.69° | 0.9965 | landmarks (18 pairs; 0.33 expected by chance) |
| FTIR (Bruker) → Raman | +1.62° | 1.0021 | landmarks (17; 0.28) |
| LDIR → Raman | −178.45° | 1.0384 | landmarks (16; 0.32) |

Without a transform no FTIR particle matches its Raman partner (0 pairs): the
instruments' frames differ by ~180° and millimetres of offset. ICP, on top of
the landmark transform, leaves the matched pairs unchanged on this sparse
filter but tightens them (PerkinElmer median distance 12.9 → 8.7 µm, Bruker
19.0 → 15.8 µm, measured at the former 100 µm gate).

The LDIR scale (+3.8%) reflects the LDIR overlay's scale calibration from the
Excel areas; a < 5% size discrepancy is accepted for LDIR.

**Stability** (`tools/validation.R`, one parameter changed at a time, every
pair re-aligned): landmark tier on/off (forcing RANSAC vs global
registration), landmark size 50–200 µm, Raman size filter 10–40 µm, TPS
on/off, RANSAC inlier distance 100–300 µm with Tier 2 forced. **Every variant
keeps 100% of the baseline pairs.** Different coarse starts (landmark vs
RANSAC; coarse rotations differing by up to 0.6°) converge to the same final
transform after ICP. (Tier 2 with the default geometric anchors, i.e. no
material anchor filter.)

**Synthetic** (`tests/testthat/test-pipeline-synthetic.R`): 120 particles in
the Raman frame; FTIR (PerkinElmer) rotated 178° and scaled 0.997, written in
millimetres, with 4 particles 20 mm outside the field; FTIR (Bruker) rotated
3° with image-style (downward) y; LDIR rotated −90° and scaled 1.03. Recovered:
rotations within 0.2°, scales within 0.2% (LDIR 0.5%); every matched pair is a
true correspondence (100 of 100 per instrument).

**Every instrument combination** (`test-pipeline-synthetic.R`): the same
synthetic filter run with each of the 11 combinations of two to four
instruments. In every run each instrument is paired with the reference and
every true correspondence is found. Without Raman, LDIR particles outside the
reference instrument's field can pair by chance at the 150 µm LDIR gate: one
chance pair in each of the three runs where LDIR is aligned to an FTIR, as
expected at this density. Raman covers the whole synthetic field, so with
Raman there are no such particles.

## 2. Matching and chance matches

Only particles that pass their instrument's quality cutoff are matched (FTIR
AAU ≥ 0.70, Raman HQI ≥ 70, LDIR ≥ 0.85, the cutoffs of the agreement
scoring) and whose Feret max is at least 20 µm, the common detection limit;
a pair whose Feret max differ by more than 5× is not made. Matching is
one-to-one (Hungarian) within a distance gate: 50 µm, and for LDIR against
the reference 3σ of the LDIR positional error fitted to that run's alignment
(within 10–50 µm; 20 µm without a fit), wider for large particles (0.15 ×
major axis). The fit uses the nearest-neighbour distances (a true-pair peak
with a 2-D normal error plus chance neighbours at the reference density):
σ = 7.1 µm on dataset A (gate 21 µm), 7.3 µm on field sample A with the
image-blob join, and 3.0 µm with the LDIR software's analysed image (gate
10 µm). A pair of two aligned instruments (neither the
reference) carries both alignment errors and takes the wider gate, 50 µm. The cost is the centre distance in µm plus unitless size
and shape terms (0.3·|ln area ratio| + 0.3·|ln Feret ratio| + 0.1·relative
aspect-ratio difference) that add less than 1.5 µm even for a tenfold size
difference: **matching is effectively by distance**, with size and shape only
breaking near-ties. For the paper (Methods):

> Instruments were aligned by a similarity transform found from large
> particles (≥ 100 µm) and fibres, accepted when the landmark pairs exceeded
> the number expected by chance at the landmark density (Poisson p < 10⁻⁸);
> otherwise from the plastic particles passing each instrument's quality
> threshold, followed by iterative-closest-point refinement. Only particles
> passing their instrument's quality threshold (FTIR AAU ≥ 0.70, Raman HQI
> ≥ 70, LDIR ≥ 0.85) and at least 20 µm in Feret diameter were paired,
> excluding pairs whose Feret diameters differed by more than a factor of 5, one-to-one with the Hungarian algorithm,
> minimising the summed centre-to-centre distance within a 50 µm gate (for
> LDIR, three times its positional error fitted to each alignment, 10–50 µm). An alignment was used only when its anchor pairs exceeded what the
> pose search finds by chance (p < 0.01). Size and shape enter the cost only as small
> tie-breaking terms (< 1.5 µm equivalent), so pairing is by position, not by
> identification. For each instrument pair, the number of pairs expected by
> chance was estimated from the density of the reference particles within the
> gate and reported with the results.

**Dataset A**

| Pair | matched | distance median / p90 / max (µm) |
|---|---|---|
| FTIR (PerkinElmer)–Raman | 31 | 8.6 / 20.6 / 61.5 |
| FTIR (Bruker)–Raman | 18 | 16.6 / 34.4 / 43.3 |
| LDIR–Raman | 28 | 7.3 / 20.7 / 58.4 |
| LDIR–FTIR (PerkinElmer) | 31 | 6.6 / 11.1 / 21.9 |
| FTIR (Bruker)–FTIR (PerkinElmer) | 16 | 20.6 / 32.9 / 49.7 |
| FTIR (Bruker)–LDIR | 14 | 16.9 / 27.7 / 42.1 |

(Maxima above the base gate are large particles: the gate widens to 0.15 ×
their major axis.) Before the field-sample changes (gates 100 / 150 µm, all
particles matched) the three Raman pairs were 32 / 29 / 34 pairs with 32 / 21
/ 31 scored; the quality filter removes only pairs that were never scored,
the tighter gate one scored pair per instrument. The 20 µm minimum and the
5× size limit then removed two FTIR (Bruker)–Raman pairs and one LDIR–Raman
pair (Raman partners of 6–7 µm; a Bruker particle sized 18 µm that the
other instruments size 98–120 µm). The 20 µm LDIR gate removed one more
LDIR–Raman pair (28 at every gate from 15 to 40 µm); FTIR (Bruker)–LDIR,
now at 50 µm instead of 75 µm, has 14.

**How many of these could be chance?** The aligned source cloud is displaced
so that true correspondences are destroyed while density and overlap are kept,
and re-matched (200 displacements each): a random rotation about the Raman
centroid (excluding rotations that move the median particle by less than two
gates), a toroidal shift of at least two gates, and the Poisson estimate
n·(1 − e^(−ρπr²)).

| Pair | observed | chance: rotation / shift / Poisson | chance share |
|---|---|---|---|
| FTIR (PerkinElmer)–Raman, 50 µm | 31 | 0.5 / 0.4 / 0.5 | ≤ 2% |
| FTIR (Bruker)–Raman, 50 µm | 18 | 0.3 / 0.3 / 0.4 | ≤ 2% |
| LDIR–Raman, 20 µm | 28 | 0.2 / 0.2 / 0.1 | ≤ 1% |

The chance share is the largest null estimate over the observed count: an
upper bound on the fraction of pairs that a random arrangement would also
produce at this gate. It grows with particle density and gate, so it must be
re-estimated for denser samples -- every run now does (see Field samples).

**Gate sensitivity** (alignment fixed), from half to twice the gate: FTIR
(PerkinElmer) 30–32 pairs, LDIR 28–31 (≥ 96% of the baseline pairs kept at
every gate); FTIR (Bruker) 14–21 (78% kept at half the gate: its pairs sit
further apart, median 16.6 µm).

**Pair confidence and detection rates.** The positional-error fit also gives
each pair its probability of being the true partner (`p_true`): the 2-D
normal density at the pair's distance against the chance nearest-neighbour
density at the reference density, weighted by the fitted share of true pairs.
On the synthetic field sample (`test-pipeline-field.R`) the expected number of
wrong pairs, Σ(1 − p_true), was 2.5 against 2 actual. Detection rates
(`07_diagnostics/detection_rates.csv`) give, per instrument pair and
direction, the share of one instrument's particles that the other also found,
by size class, with Wilson 95% intervals; dataset A, all sizes: LDIR finds
78% (62–88%) of the Raman particles in its area, Raman 90% (75–97%) of the
LDIR particles.

## 3. LDIR image join

LDIR exports no coordinates; the pipeline finds the particles in the LDIR
image and attaches each Excel row to an image particle by size and shape.
Each run checks the join against Raman (`R/ldir_join_check.R`): for each LDIR
particle matched to Raman, the LDIR table's size and family are compared with
the Raman partner's, against the same with the LDIR rows randomly re-assigned
(1000 permutations).

Dataset A (28 pairs): size Spearman **0.97** (random join 0.00, 95th
percentile 0.31; p < 0.001; 0.79 with 30 pairs before the 5× size limit
removed a 318 µm LDIR particle paired with a 6 µm Raman one). The family statistic is not informative here
(100% vs 91% for a random join — the sample is almost all PET). The size
statistic is biased upward because the matcher slightly prefers similar sizes;
the null keeps the matched positions, so the bias does not make a random join
look correct. The manual swap file of dataset A corrects 10 of 40 joins.

Synthetic check (`test-validation-tools.R`): a correct join is significant on
both statistics, a scrambled one on neither.

## 4. Material agreement

Scored only for pairs where both identifications pass their instrument's
cutoff (FTIR AAU ≥ 0.70, Raman HQI ≥ 70, LDIR ≥ 0.85) and neither is unknown.
Primary = family agreement, secondary = exact polymer.

| Pair | matched | excluded (quality) | scored | Exact / Family / Filler / Unresolved / Disagree |
|---|---|---|---|---|
| FTIR (PerkinElmer)–Raman | 31 | 0 | 31 | 31 / 0 / 0 / 0 / 0 |
| FTIR (Bruker)–Raman | 18 | 0 (16 FTIR particles below 0.70 not matched) | 18 | 18 / 0 / 0 / 0 / 0 |
| LDIR–Raman | 28 | 0 (8 LDIR particles below 0.85 not matched) | 28 | 28 / 0 / 0 / 0 / 0 |

Triplets (Raman particles matched by FTIR PerkinElmer and LDIR): 28, all
full agreement.

These rates say little about the method's discrimination: dataset A is one
polymer. The synthetic run checks the scoring itself: designed family-only
differences, disagreements, low-quality and unidentified particles are all
counted exactly (`test-pipeline-synthetic.R`).

## 5. Field samples

A field-sample filter is dense: the Raman file of Faanui SW1L R4 B has
14 562 particles on a 13 mm filter (2 960 with HQI ≥ 70). With all particles,
every LDIR particle has a chance Raman partner within 150 µm (> 99.9%), every
pose of a coarse search scores alike, and pairing measures nothing; a first
four-instrument run of sample A (12 513 Raman particles) gave a physically
implausible LDIR pose (scale 0.79), 0% family agreement on every pair and an
LDIR join check no better than random. Hence the quality filter before
matching, the tighter gates, Tier 2 anchors on quality-passing plastics, the
significance criterion for landmarks, and the chance check in every run.

Test with known truth: an LDIR table derived from the sample B Raman
particles ≥ 50 µm (1 243, rotated −90°, scale 1.03, 5 µm noise). Landmarks:
227 pairs against 4.7 expected by chance. LDIR–Raman at the 75 µm gate: 579
pairs, of which 345 (60%) wrong -- the LDIR particles whose true partner is
below HQI 70 pair with another particle -- and the run's chance estimate said
61%. Wrong pairs sit at a median 50 µm, true ones at 8 µm: within 25 µm 9%
of the pairs are wrong, within 15 µm 2%. On field samples the LDIR gate
therefore has to follow the positional accuracy actually achieved; the
chance estimate says per run whether the pairs can be interpreted.

Field sample A (Faanui SW1L R4 A; 12 513 Raman particles, 2 387 with
HQI ≥ 70; LDIR 1 409, 463 with quality ≥ 0.85). LDIR landmarks: 26 pairs
against 0.9 expected by chance (p = 2e-28); ICP residual 10.5 µm. At the
75 µm gate 170 LDIR–Raman pairs: 68 within 10 µm, 92 within 15 µm, then a
flat tail. The positional-error fit gives σ = 7.3 µm and ~100 true pairs;
expected share of wrong pairs 41% at 75 µm, 23% at 50 µm, 5% at 20 µm, 3% at
15 µm. Hence an LDIR gate of about 3σ, now set per run (with the analysed
LDIR image σ = 3.0 µm and the gate 10 µm). Within 15 µm the sizes agree (Spearman
0.56, tail 0.28) but the identities do not: LDIR says cellulose / natural
polyamide / chitin, Raman (HQI 70–98) tin(II) chloride, p-diaminobenzene,
rutin and similar -- a question for the libraries, not the matching.
Both FTIR exports had been filtered to plastics in the instrument software
(15 and 5 particles in a 3.5 mm patch); against 2 387 Raman particles that
cannot fix a pose (PerkinElmer came out mirrored, the two Bruker searches
disagreed by 140°), and the pose check now reports such an alignment as not
determined.

## 6. Images

- **FTIR image** placement by registering the particles visible in the image
  to the table: dataset A, 2.0 µm/px.
- **Raman image** from the WITec placement file (`<image>.txt`): on the field
  sample Faanui SW1L R4 B (14,562 particles), with X/Y read as the top-left
  corner in the stage frame, all particles lie on the image and the image
  brightness at the particle positions peaks at zero offset and exactly scale
  1.00 (it halves at ±1% scale); every other reading puts 0–22% of particles
  on the image. The auto-calibration (without the file) found no confident fit
  on this sample.
- **LDIR**: see Known limitations.

## 7. Reproducibility

Two runs of the same input give byte-identical output folders apart from
timestamps, the run's own path and PDF creation dates
(`test-pipeline-synthetic.R`). The aligners use fixed seeds and restore the
caller's random-number state. Each run records its full configuration, R
session, package versions, input files (with MD5) and git commit
(`00_manifest/`). The dataset A results are pinned by a regression test.

## Known limitations

- **One real lab dataset**, essentially single-polymer; field samples are
  dense; each run reports its chance-pair share and positional error σ, and
  warns when the gate is well below 3σ.
- **LDIR scan circle**: the detector overestimates the circle radius (~12% on
  a synthetic disk) and only works when the circle covers < 20% of the
  image; the radius sets µm/px on optical LDIR images. Dataset A uses the
  analyzed overlay instead, so its results are unaffected. Not yet tested on
  a real image with a visible circle.
- **Raman spectral-library index** (trade-name resolution) is used when the
  licensed library listings are present and recorded in the manifest; not
  part of the public deposit.
