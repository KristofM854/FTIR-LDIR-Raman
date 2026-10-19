# Material tables

One table per instrument role (`ftir`, `ftir_bruker`, `raman`, `ldir`). Each row
translates one identification label, exactly as that instrument exports it,
into the name the pipeline should classify:

```csv
label,material
Nylon 6/6 copolymer,PA66
Hostalen GM 7040,Polyethylene
```

- Matching ignores case and surrounding spaces; everything else must match.
- The exported label is kept in the `identification_raw` column; `material`
  holds the translated name, from which `material_polymer` and
  `material_family` are derived (`R/08b_material_map.R`).
- Labels not in the table are used unchanged. The built-in rules recognise
  the vocabulary of the bundled example dataset; `ftir.csv` adds the
  PerkinElmer library spelling "PolyEthylen" (field samples, 2026).
- To use your own table without editing these files, set
  `config$<role>_material_table` to a CSV path (or a data frame), or
  `config$material_table_dir` to a folder holding `<role>.csv` files.

Use a table for labels the rules cannot know: trade names, library-specific
spellings, or a lab's own abbreviations. Check which labels end up
"Unknown" or "Other polymer" with `tools/material_vocabulary.R`.
