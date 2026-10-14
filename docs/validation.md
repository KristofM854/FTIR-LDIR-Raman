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
(rotation, uniform scale, translation): landmarks (large particles and
fibres), else RANSAC / global registration, then ICP.

**Dataset A**

| Pair | rotation | scale | tier used |
|---|---|---|---|
| FTIR (PerkinElmer) → Raman | −178.69° | 0.9965 | landmarks (18 / 24 inliers) |
| FTIR (Bruker) → Raman | +1.62° | 1.0021 | landmarks (17 / 20) |
| LDIR → Raman | −178.45° | 1.0384 | landmarks (16 / 23) |

Without a transform no FTIR particle matches its Raman partner (0 pairs): the
instruments' frames differ by ~180° and millimetres of offset. ICP, on top of
the landmark transform, leaves the matched pairs unchanged on this sparse
filter but tightens them (PerkinElmer median distance 12.9 → 8.7 µm, Bruker
19.0 → 15.8 µm).

The LDIR scale (+3.8%) reflects the LDIR overlay's scale calibration from the
Excel areas; a < 5% size discrepancy is accepted for LDIR.

**Stability** (`tools/validation.R`, one parameter changed at a time, every
pair re-aligned): landmark tier on/off (forcing RANSAC vs global
registration), landmark size 50–200 µm, Raman size filter 10–40 µm, TPS
on/off, RANSAC inlier distance 100–300 µm with Tier 2 forced. **Every variant
keeps 100% of the baseline pairs.** Different coarse starts (landmark vs
RANSAC; coarse rotations differing by up to 0.4°) converge to the same final
transform after ICP.

**Synthetic** (`tests/testthat/test-pipeline-synthetic.R`): 120 particles in
the Raman frame; FTIR (PerkinElmer) rotated 178° and scaled 0.997, written in
millimetres, with 4 particles 20 mm outside the field; FTIR (Bruker) rotated
3° with image-style (downward) y; LDIR rotated −90° and scaled 1.03. Recovered:
rotations within 0.2°, scales within 0.2% (LDIR 0.5%); every matched pair is a
true correspondence (100 of 100 per instrument).

## 2. Matching and chance matches

Matching is one-to-one (Hungarian) within a distance gate: 100 µm, 150 µm for
pairs involving LDIR.

**Dataset A**

| Pair | matched | distance median / p90 / max (µm) |
|---|---|---|
| FTIR (PerkinElmer)–Raman | 32 | 8.7 / 21.1 / 61.5 |
| FTIR (Bruker)–Raman | 29 | 15.8 / 41.3 / 56.8 |
| LDIR–Raman | 34 | 8.5 / 40.9 / 113.0 |
| LDIR–FTIR (PerkinElmer) | 33 | 6.5 / 11.0 / 21.9 |
| FTIR (Bruker)–FTIR (PerkinElmer) | 25 | 20.5 / 44.5 / 54.2 |
| FTIR (Bruker)–LDIR | 29 | 20.3 / 56.2 / 125.8 |

**How many of these could be chance?** The aligned source cloud is displaced
so that true correspondences are destroyed while density and overlap are kept,
and re-matched (200 displacements each): a random rotation about the Raman
centroid (excluding rotations that move the median particle by less than two
gates), a toroidal shift of at least two gates, and the Poisson estimate
n·(1 − e^(−ρπr²)).

| Pair | observed | chance: rotation / shift / Poisson | chance share |
|---|---|---|---|
| FTIR (PerkinElmer)–Raman, 100 µm | 32 | 2.9 / 1.7 / 2.1 | ≤ 9% |
| FTIR (Bruker)–Raman, 100 µm | 29 | 3.1 / 2.0 / 2.6 | ≤ 11% |
| LDIR–Raman, 150 µm | 34 | 5.8 / 4.1 / 5.4 | ≤ 17% |

The chance share is the largest null estimate over the observed count: an
upper bound on the fraction of pairs that a random arrangement would also
produce at this gate. It grows with particle density and gate, so it must be
re-estimated for denser samples.

**Gate sensitivity** (alignment fixed): from half to twice the gate, the
matched count changes by at most 2 (31–33, 27–31, 33–34) and ≥ 93% of the
baseline pairs are kept at every gate.

## 3. LDIR image join

LDIR exports no coordinates; the pipeline finds the particles in the LDIR
image and attaches each Excel row to an image particle by size and shape.
Each run checks the join against Raman (`R/ldir_join_check.R`): for each LDIR
particle matched to Raman, the LDIR table's size and family are compared with
the Raman partner's, against the same with the LDIR rows randomly re-assigned
(1000 permutations).

Dataset A (34 pairs): size Spearman **0.78** (random join 0.00, 95th
percentile 0.28; p < 0.001). The family statistic is not informative here
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
| FTIR (PerkinElmer)–Raman | 32 | 0 | 32 | 32 / 0 / 0 / 0 / 0 |
| FTIR (Bruker)–Raman | 29 | 8 (FTIR) | 21 | 21 / 0 / 0 / 0 / 0 |
| LDIR–Raman | 34 | 3 (LDIR) | 31 | 31 / 0 / 0 / 0 / 0 |

Triplets (Raman particles matched by FTIR PerkinElmer and LDIR): 32; 30 full
agreement, 2 FTIR+Raman agree (the LDIR side below 0.85, relabelled unknown).

These rates say little about the method's discrimination: dataset A is one
polymer. The synthetic run checks the scoring itself: designed family-only
differences, disagreements, low-quality and unidentified particles are all
counted exactly (`test-pipeline-synthetic.R`).

## 5. Images

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

## 6. Reproducibility

Two runs of the same input give byte-identical output folders apart from
timestamps, the run's own path and PDF creation dates
(`test-pipeline-synthetic.R`). The aligners use fixed seeds and restore the
caller's random-number state. Each run records its full configuration, R
session, package versions, input files (with MD5) and git commit
(`00_manifest/`). The dataset A results are pinned by a regression test.

## Known limitations

- **One real dataset**, essentially single-polymer; the chance-match share and
  the agreement rates must be re-established on each new sample type.
- **LDIR scan circle**: the detector overestimates the circle radius (~12% on
  a synthetic disk) and only works when the circle covers < 20% of the
  image; the radius sets µm/px on optical LDIR images. Dataset A uses the
  analyzed overlay instead, so its results are unaffected. Not yet tested on
  a real image with a visible circle.
- **Raman spectral-library index** (trade-name resolution) is used when the
  licensed library listings are present and recorded in the manifest; not
  part of the public deposit.
