# Correspondence with the original analysis folders

The analysis was first run in one folder per stage (`EtapaN_...`, Portuguese
names). This repository reorganizes the same code into one project. Every
migrated script was re-run on the original data and checked against the
original outputs (identical tables, or numerical differences < 1e-12).

| Original script(s) | Repository script | Key output (original -> repository) |
|---|---|---|
| Etapa 1 Excel; Etapa6 `08_Verificacoes` (V1); `09_Verificacao_V1b` | `01_candidates/01_candidates_loci_positions.R` | `6B2/06_locus_map_39` -> `01_candidates/02_locus_map`; `V/06_V1b_snp_verdict` -> `06_snp_verdict` |
| Etapa2 `01`-`04` | `02_maturation_timecourse/01_qc_profiles_GSE98199.R`, `02_deseq2_contrasts.R` | `09_stage2_candidate_evidence_final` -> `09_candidate_temporal_evidence`; `05_DESeq2_all_genes_two_contrasts` -> `07_deseq2_all_genes` |
| Etapa3 `01`-`04` | `03_expression_atlas/01_atlas_bodypart_summary.R`, `02_seed_context_study_support.R` | `15_STAGE3_FINAL_evidence_39_candidates` -> `10_atlas_evidence_39` |
| Etapa4 `00b` | `04_coexpression_network/01_combine_atlas_downloads.R` | combined TPM matrix -> `data/processed/SEA_seedrelated_gene_TPM.tsv.gz` |
| Etapa4 `01`, `02` | `02_sample_audit_curation.R` | `14_sample_curation_decisions` -> `03_sample_curation_decisions` |
| Etapa4 `03` (4C) | `03_connectivity_qc_confounding.R` | `21_sample_connectivity_ZK_global` -> `09_sample_connectivity_ZK` |
| Etapa4 `04`, `05` (4D, 4D2) | `04_gene_filter_soft_threshold.R` | `36_SFT_variance_grid_summary` -> `17_soft_threshold_grid`; `04D2_genes_top75_variance` -> `20_network_genes` |
| Etapa4 `06` (4E) | `05_network_construction.R` | `40_module_assignments_all_genes` -> `21_module_assignments`; `42` -> `23_module_eigengenes`; `44` -> `25_module_bodypart_correlations` |
| Etapa4 `07` (4F) | `06_centrality_neighbors_projection.R` | `54_candidate_top30_intramodular_neighbors` -> `31_candidate_top30_neighbors`; `57_stage4F_candidate_network_evidence_39` -> `33_candidate_network_evidence_39` |
| Etapa6 `08_Verificacoes` (V2) | `07_organellar_artifact_check.R` | `V/03_V2_organelle_fraction_by_module` -> `34_organelle_fraction_by_module` |
| Etapa5 `01` (5A) | `05_functional_enrichment/01_go_map_module_enrichment.R` | `02_GO_module_enrichment_all` (same name) |
| Etapa5 `01` (neighbourhoods) + `03_V3` | `02_neighborhood_enrichment.R` | `04_GO_..._within_module` -> `05_GO_neighborhood_within_module`; `V3/02_summary_by_candidate` -> `08_BH_family_summary_by_candidate` |
| Etapa5 `02` (5B) | `03_redundancy_functional_profile.R` | `14_candidate_functional_profile_39` -> `12_candidate_functional_profile_39` |
| Etapa6 `01` (6A) | `06_evidence_integration/01_integrated_matrix_primary_ranking.R` | `6A/02_integrated_evidence_matrix_39` (same name, folder `01_integrated_matrix`) |
| Etapa6 `02` (6B) | `02_ranking_sensitivity.R` | `6B/05_stage6B_complete_case_layer_scores_30` -> `01_complete_case_layer_scores_30` |
| Etapa6 `03` (6B2) | `03_convergence_weights_loci.R` | `03_convergence_permutation_test`, `04_weakest_link_2domains` (same names) |
| Etapa6 `04` (6C) | `04_network_maturation_anchor.R` | `05_candidate_module_anchor_39` (same name) |
| Etapa6 `05` (6C2) | `05_module_tissue_within_project.R` | `02_module_tissue_within_vs_original` (same name) |
| Etapa6 `06` (6D) | `06_final_evidence_classes.R` | `01_final_evidence_table_39_genes` (same name) |
| Etapa6 `10_..._sensibilidade_V1b` | `07_sensitivity_displaced_genes.R` | `6S_V1b/01_convergence_by_scenario` (same name) |
| Etapa6 `11_Etapa6E` | `08_a4_genes_at_displaced_snps.R` | `6E/01_genes_a4_vs_autores_temporal_rede` -> `01_authors_vs_a4_genes` |
| Etapa7 `01` (7A) | `07_external_validation/01_pereira_lima_longevity_lists.R` | `7A/02_candidate_enrichment_S5` (same name) |
| Etapa7 `02` (7B) | `02_medicago_longevity_network.R` | `7B/04_module_enrichment_by_unique_AGI` (same name) |

Not migrated: `Etapa6/07_Etapa6_figuras_banner.R` (figures made only for the
IX SBV poster) and `go_bh_check.R` (an audit check now covered by
`05_functional_enrichment/02`).

Labels changed to English in the outputs: SNP verdict `consistente/borda/deslocado`
-> `consistent/borderline/displaced`; locus label `semPosicao` -> `noPosition`.
