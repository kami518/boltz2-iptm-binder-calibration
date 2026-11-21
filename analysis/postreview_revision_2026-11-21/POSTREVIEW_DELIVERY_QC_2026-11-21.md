# Post-review delivery QC

Status: **passed**

## Final files

| File | SHA-256 |
|---|---|
| `publication/manuscript/main_article.md` | `3999d58ea30ec8c62108dfbdd154b00c48c0fae7dfa2edc4edc759e49a11cf00` |
| `publication/manuscript/supplementary_materials.md` | `0d6a00c634b11d6edce8ea0030de1dad5090ce2dc4a067ab7965ec95271e00b2` |
| `publication/manuscript/PEDS_Main_and_Supplement_Review_Draft_v1.1.0.pdf` | `aaf638d0c1c162ee005351bf0bc8540de4698f62a6c79f281f963fa9574eae6c` |

The submission-candidate PDF is A4, contains 57 pages, and is 5,708,787 bytes.

## Scientific and numerical checks

- C01 official ipSAE recomputation passed for 2,074/2,074 pairs.
- Released ipSAE-field agreement is reported descriptively: 1,048/2,074 for the directional minimum and 821/2,074 for the other field.
- C01B combination outputs contain the complete six-target U/I/F table under both weightings.
- C01B did not run a random forest or add another feature.
- C02 contains 2,000 valid paired cluster-bootstrap replicates for each of the primary-six, primary-five, and external datasets, with zero failed replicates.
- The key relative reductions are locked as 45.1% (39.5%-49.8%), 42.9% (38.4%-46.8%), and 25.9% (19.3%-31.9%).
- The probability estimand is explicitly limited to experimental binding among designs already selected and submitted for assay.
- The external transport analysis is described as changing target/campaign and score-generation conditions together.
- Author identity, ORCID, affiliation, corresponding email, CRediT roles, and the approved acknowledgement are present.

## Manuscript checks

- Title is identical in the main manuscript, supplement, and PDF metadata.
- Abstract is 9 sentences and 240 words.
- The canonical `post hoc two-feature logistic model` wording is present.
- `prespecified two-feature` and `universal selection threshold` are absent.
- The original Anthropic technical report and the dataset release are cited separately.
- 41 references are alphabetically ordered; 35 unique DOI-bearing records are present.
- The eight methodological sources added after the reference-density audit were verified through Crossref or the official JMLR page.
- S1-S16 retain their original order; Table S17 and Supplementary Figure S11 are present.
- Cohort accounting now reconciles 4,943 main-source pairs and the 1,320 external interpretable pairs without unexplained rows.

## Figure checks

- Main figures retain numbering 1-7.
- Main Figures 1-7 and Supplementary Figure S11 each exist as PNG, SVG, and PDF.
- Figure 2 uses the corrected `Target prevalence` label; Figure 5 and Supplementary Figure S7 remove the retired `pilot value` wording.
- All 24 production files exist and match the hash manifest.
- No `RQ1`-`RQ6` internal work labels occur in the generated SVG panels.
- C02 intervals appear in Figure 4 and Table 3.
- C01, C01B, and C02 results appear in Supplementary Figure S11 and Table S17.

## PDF checks

- All 39 main/supplementary figure and table titles are present.
- All 18 figure images are embedded on their corresponding title pages.
- Page numbering is continuous from Page 2 through Page 57.
- All pages are A4.
- Bookmarks and internal destinations resolve.
- Internal and external link annotations were checked; none were invalid or unclassified.
- The PDF contains no `file:` URI, local absolute path, retired pilot wording, or review-draft footer.
- Key sample sizes, interval endpoints, and C01B differences are present in extracted text.
- No raw Markdown, replacement glyphs, missing figure image, or text-empty body page was detected.
- All figure and table pages were rendered and visually inspected for clipping, overlap, or unreadable small text.

## Remaining human tasks

- Apply the current Protein Science journal-specific file format and length limits.
- Obtain mainline approval, then push the clean source snapshot and create the `v1.1.0` GitHub Release.
- Upload the verified release pieces and archive to Zenodo; activate DOI `10.5281/zenodo.24400262` only after final approval.
