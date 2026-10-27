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

- **Natural-material references for Raman (for the Raman analysts).** On
  field sample A, FTIR and LDIR identify most particles as cellulose, protein,
  natural polyamide or chitin; Raman reports reference chemicals instead
  (79 of 108 FTIR-Raman pairs: p-diaminobenzene 17, rutin 11,
  1,5-pentanediol 10, essential oils, esters, paraffins) at a median HQI of 84,
  higher than its cellulose hits, so the HQI cutoff does not remove them. The
  nine S.T. Japan libraries (21 105 entries, ~10 700 substances) hold almost no
  natural-fibre references: one cotton, ~12 native celluloses (the other 120
  cellulose entries are derivatives), two chitin, no wool, silk, keratin,
  linen, hemp, jute or lignin. Where a reference exists Raman finds it (all 5
  Raman cellulose hits agree with FTIR). Since this change such hits count as
  unresolved, not as disagreements (field sample A: family agreement 26%
  (13-45%) of 27 resolved pairs instead of 7% of 107; 81 unresolved). The fix
  belongs on the instrument side: search only the relevant libraries
  (L60035 microplastics and related compounds plus natural materials) and/or
  add references measured in-house -- cotton, wool, viscose, native cellulose,
  chitin, and biofilm or plankton material from the sampling sites.
  Also to check: tin(II) chloride (from the inorganics library) is the Raman
  hit for 5 cellulose/protein particles and still counts as a disagreement.
- **PET vs PC.** On the polycarbonate reference filter, FTIR (PerkinElmer)
  calls 8 particles PET where Raman and LDIR call them PC. To check: the raw
  identification labels (`identification_raw`), since the classification
  counts any label containing "polyester" as PET.
- **File naming** across the instruments, so that the files of one filter are
  recognisable as such.

## Temporary

- Segmented PerkinElmer exports are stitched automatically
  (`R/ftir_segments.R`, called at the start of `run_pipeline()`;
  `tools/stitch_ftir_segments.R`; `tests/testthat/test-ftir-segments.R`).
  Remove these once the Spotlight software update exports one table and one
  image per scan; the placement file reader (`read_ftir_image_metadata()`)
  can stay.

## Release

- Zenodo deposit with a DOI; complete the citation metadata (`CITATION.cff`,
  `.zenodo.json`) and add the method paper once published.
- Re-record the package versions (`renv.lock`) and the Python versions
  (`inst/python/requirements.txt`) of the analysis environment.
