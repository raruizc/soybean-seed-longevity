# Raw data (not versioned)

All inputs are public or come from the curation of the two GWAS papers. They
are not redistributed in this repository: download each file and place it in
the folder shown, with the name shown.

| Folder | File | Source |
|---|---|---|
| `gwas/` | `gwas_candidates_curated.xlsx` | Curated candidate list from Wu et al. (2024) and Thenveettil et al. (2025), with IDs converted to Wm82.a4. Sheets `Master_Candidates` (one row per source gene) and `Final_Loci` (one row per Wm82.a4 gene) |
| `gwas/` | `candidate_aliases.csv` | Older IDs of each candidate, from the Legume Information System pangene table (`Glycine.pan5.MKRS.table_ref_lines.tsv.gz`) |
| `annotation/` | `glyma.Wm82.gnm4.ann1.T8TQ.info_annot.txt.gz` | SoyBase/LIS Data Store, `Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ/` |
| `annotation/` | `glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz` | Same folder |
| `GSE98199/` | `GSE98199_processed.txt.gz` | NCBI GEO series GSE98199 (Pereira Lima et al., 2017), supplementary count matrix |
| `GSE98199/` | `metadata_GSE98199.csv` | Built from the GEO series page: columns `GEO_accession, description, stage, treatment, replicate, DAF, cultivar` (11 libraries) |
| `atlas/` | `SEA_39genes_TPM.tsv`, `SEA_39genes_metadata.tsv` | Soybean Expression Atlas v2 (Almeida-Silva et al., 2023), "Search gene list" with the 39 IDs written by `scripts/03_expression_atlas/01`, all body parts, TPM |
| `network/` | `SEA_<part>_TPM.tsv.gz`, `SEA_<part>_metadata.tsv` for `seed`, `embryo`, `cotyledon`, `endosperm`, `seed_coat`, `suspensor` | Soybean Expression Atlas v2, "Download by body part", TPM, gene level |
| `external/pereira_lima_2017/` | `S3Table.xlsx`, `S4Table.xlsx`, `S5Table.xlsx` | Supplementary material of Pereira Lima et al. (2017), PLoS ONE |
| `external/righetti_2015/` | `Righetti2015_DS1_MatNet.xls`, `..._DS2_TF_MatNet.xls`, `..._DS3_longevity_module.xls`, `..._DS4_desiccation_tolerance.xls`, `..._DS5_ABI3_regulation.xls` | Supplemental Data Sets 1 to 5 of Righetti et al. (2015), Plant Cell (PMC4682330), renamed as shown |
| `external/qtl/` | `storability_qtl_curated.csv` | Curated by hand from published soybean QTL and GWAS studies of seed storability, longevity, viability, germination and vigor (one row per QTL, Wm82.a4 coordinates; template and column rules in `docs/qtl_table_template.csv`). Optional: `run_all.R` skips stage 7c while it is missing |

## References

- Almeida-Silva F, Pedrosa-Silva F, Venancio TM (2023). The Soybean Expression Atlas v2: a comprehensive database of over 5000 RNA-seq samples. The Plant Journal. doi:10.1111/tpj.16459
- Pereira Lima JJ et al. (2017). Molecular characterization of the acquisition of longevity during seed maturation in soybean. PLoS ONE 12:e0180282. doi:10.1371/journal.pone.0180282
- Righetti K et al. (2015). Inference of longevity-related genes from a robust coexpression network of seed maturation identifies regulators linking seed storability to biotic defense-related pathways. The Plant Cell 27:2692-2708. doi:10.1105/tpc.15.00632
- Thenveettil N et al. (2025). Genome-wide association study reveals the QTLs and candidate genes associated with seed longevity in soybean (Glycine max (L.) Merrill). BMC Plant Biology 25:829. doi:10.1186/s12870-025-06822-1
- Wu et al. (2024). Identification of candidate genes for soybean storability via GWAS and WGCNA approaches. Agronomy 14:2457. doi:10.3390/agronomy14112457
