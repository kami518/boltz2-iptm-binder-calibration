# Release v1.2.1

This release adds the CSBJ revision to the public reproducibility record without
changing the earlier `v1.1.0` DOI. It includes the final figure-free manuscript,
landscape tables, Supporting Information, 21 figure sets in PDF and SVG with
600-dpi PNG previews, panel-level traceability, and the accepted target-as-new
run.

`v1.2.1` supersedes the first `v1.2.0` upload after the final vector audit
removed remaining default Matplotlib colors and embedded raster elements from
the requested figure scope. Numerical results and manuscript conclusions are
unchanged.

## Figure delivery

- Main Figures 1-8 and Supplementary Figures S1-S10 use a maximum width of
  180 mm, DejaVu Sans, editable SVG text, and the Okabe-Ito colorblind-safe
  palette.
- Supplementary Figures S11-S13 are retained as separate vector assets.
  S12-S13 were outside the requested S1-S10 unification pass and keep their
  original physical dimensions.
- Error bars, legends, panel labels, and axis text were visually checked in the
  generated PDF, SVG, and 600-dpi PNG files.

## Traceability

`manuscript/csbj_revision_2026-11-26/FIGURE_TRACEABILITY_2026-11-26.csv`
contains 30 panel-level rows. Each row names the accepted source file, run
directory, hash manifest, generation script, and resulting figure files. The
directly traced CSV sources, run metadata, QC records, and scripts are included
in this archive.

## Manuscript files

The submission DOCX files contain captions and tables but no embedded figures,
DrawingML, VML, OLE objects, or `word/media` entries. Matching PDFs contain no
image XObjects. Figures are submitted as separate PDF and SVG files.

## Verification

Run:

```bash
python analysis/verify_csbj_release_v1_2_1_2026-11-26.py
```

The command verifies the manifest, hashes, ZIP, document structure, figure
formats and widths, and every path recorded in the traceability table.
