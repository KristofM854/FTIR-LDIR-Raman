# Roadmap and open questions

## Validation on more samples

The method is validated on one reference filter measured by all four
instruments (dataset A), a polycarbonate reference filter and the first field
samples (see [validation.md](validation.md)). Next:

- Field samples with the LDIR analysed (labelled) image, so that every LDIR
  table row is tied to its image particle without the size-and-shape join.
- Field samples with unfiltered FTIR exports. A list already reduced to
  plastics in the instrument software is too short to fix the alignment
  against a dense Raman field (field sample A: 15 and 5 particles; both
  alignments reported as not determined).
- Further comparison filters, to repeat the chance-match null model, the
  parameter sensitivity and the detection rates per size class on each.
- A sweep of the quality cutoffs (FTIR AAU, Raman HQI, LDIR quality) and its
  effect on agreement and detection rates.

## Open questions on identification

To be settled with the analysts who run the instruments, not in code:

- **Natural particles in Raman.** On field sample A, LDIR identifies the
  closest LDIR-Raman pairs (within 15 µm, sizes agreeing) mostly as cellulose,
  natural polyamide or chitin, while Raman reports library hits such as tin(II)
  chloride, p-diaminobenzene or rutin at HQI 70-98. Does the Raman library hold
  cellulose and natural-fibre references? If not, should such hits count as
  unresolved rather than as a disagreement?
- **PET vs PC.** On the polycarbonate reference filter, FTIR (PerkinElmer)
  calls 8 particles PET where Raman and LDIR call them PC. To check: the raw
  identification labels (`identification_raw`), since the classification
  counts any label containing "polyester" as PET.
- **File naming** across the instruments, so that the files of one filter are
  recognisable as such.

## Temporary

- `tools/stitch_ftir_segments.R` joins the six-segment PerkinElmer exports.
  Remove it once the Spotlight software update exports one table and one
  image per scan.

## Release

- Zenodo deposit with a DOI; complete the citation metadata (`CITATION.cff`,
  `.zenodo.json`) and add the method paper once published.
- Re-record the package versions (`renv.lock`) and the Python versions
  (`inst/python/requirements.txt`) of the analysis environment.
