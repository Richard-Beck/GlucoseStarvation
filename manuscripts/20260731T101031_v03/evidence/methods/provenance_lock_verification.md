# Provenance Lock Verification

- Table: `agent-dev/manuscript_methods/20260730_v03_methods_provenance_reconstruction/locked_provenance_table.md`
- Root: `/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/GlucoseStarvation`
- Parsed rows: 271

## Status Summary

| status | n_rows |
|---|---:|
| ok | 234 |
| changed | 0 |
| missing | 0 |
| ambiguous | 0 |
| unresolved | 0 |
| not_checked | 37 |

## Non-OK Rows

| row | line | status | reason | id | lock_target | stored_sha256 | current_sha256 |
|---:|---:|---|---|---|---|---|---|
| 123 | 129 | not_checked | no_stored_sha256 | `raw_source:run_20260721_163410/original_shard_measurements` | `NA` | `NA` | `NA` |
| 124 | 130 | not_checked | no_stored_sha256 | `raw_source:run_20260721_163410/reprocessed_35_fields` | `NA` | `NA` | `NA` |
| 127 | 133 | not_checked | no_stored_sha256 | `unresolved:source_plate_map_pngs` | `NA` | `NA` | `NA` |
| 128 | 134 | not_checked | no_stored_sha256 | `raw_source:filename_encoded_acquisition_identity` | `NA` | `NA` | `NA` |
| 133 | 139 | not_checked | no_stored_sha256 | `named_object:segmentation_classifier_90frame/20260722#validation_bundle` | `NA` | `NA` | `NA` |
| 137 | 143 | not_checked | no_stored_sha256 | `raw_source:validation90_phase_alive_dead_tiffs` | `NA` | `NA` | `NA` |
| 139 | 145 | not_checked | no_stored_sha256 | `external:resnet18-f37072fd.pth` | `NA` | `NA` | `NA` |
| 140 | 146 | not_checked | no_stored_sha256 | `unresolved:validation_ObjectClassifier_execution_source` | `NA` | `NA` | `NA` |
| 141 | 147 | not_checked | no_stored_sha256 | `unresolved:validation_CPSAM_model_weights` | `NA` | `NA` | `NA` |
| 143 | 149 | not_checked | no_stored_sha256 | `named_object:data/train/composites#matched_training_tiff_collection` | `NA` | `NA` | `NA` |
| 144 | 150 | not_checked | no_stored_sha256 | `named_object:data/train/cellpose_inputs#matched_training_tiff_collection` | `NA` | `NA` | `NA` |
| 147 | 153 | not_checked | no_stored_sha256 | `raw_source:all_raw_acquisition_tiffs` | `NA` | `NA` | `NA` |
| 150 | 156 | not_checked | no_stored_sha256 | `external:imutils.image_utils#training_image_construction` | `NA` | `NA` | `NA` |
| 153 | 159 | not_checked | no_stored_sha256 | `external:imutils.curation_controller.CurationController` | `NA` | `NA` | `NA` |
| 154 | 160 | not_checked | no_stored_sha256 | `external:imutils.object_classification.ObjectClassifier` | `NA` | `NA` | `NA` |
| 155 | 161 | not_checked | no_stored_sha256 | `named_object:run_20260324_233122#current_selected_montage_sources` | `NA` | `NA` | `NA` |
| 156 | 162 | not_checked | no_stored_sha256 | `raw_source:selected_montage_phase_alive_dead_tiffs_and_masks` | `NA` | `NA` | `NA` |
| 170 | 176 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#single_line_stan_data_collection` | `NA` | `NA` | `NA` |
| 177 | 183 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#stan_data_by_optimization_fit` | `NA` | `NA` | `NA` |
| 178 | 184 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#per_start_optimization_result_triplets` | `NA` | `NA` | `NA` |
| 180 | 186 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#combined_optimization_result_triplets` | `NA` | `NA` | `NA` |
| 184 | 190 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#current_enabled_nuts_draw_collection` | `NA` | `NA` | `NA` |
| 187 | 193 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#current_nuts_config_collection` | `NA` | `NA` | `NA` |
| 201 | 207 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#all_lines_posterior_prediction_collection` | `NA` | `NA` | `NA` |
| 204 | 210 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#no_sum_posterior_prediction_collection` | `NA` | `NA` | `NA` |
| 206 | 212 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#selected_best_map_no_sum_models` | `NA` | `NA` | `NA` |
| 207 | 213 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#selected_best_map_transfer_cases` | `NA` | `NA` | `NA` |
| 209 | 215 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#strategy_simulation_shard_collection` | `NA` | `NA` | `NA` |
| 230 | 236 | not_checked | no_stored_sha256 | `named_object:red_a30_counts_20260722#rank_matched_morphology_initializers` | `NA` | `NA` | `NA` |
| 231 | 237 | not_checked | no_stored_sha256 | `named_object:morphology_metrics_a30_20260729#per_start_optimization_results` | `NA` | `NA` | `NA` |
| 233 | 239 | not_checked | no_stored_sha256 | `named_object:morphology_metrics_a30_20260729#selected_optimization_result_collections` | `NA` | `NA` | `NA` |
| 251 | 257 | not_checked | no_stored_sha256 | `raw_source:sum159_competition_phase_green_alive_dead_images` | `NA` | `NA` | `NA` |
| 261 | 267 | not_checked | no_stored_sha256 | `unresolved:run_20260324_233122/cellpose_default_model` | `NA` | `NA` | `NA` |
| 265 | 271 | not_checked | no_stored_sha256 | `raw_source:SUM159/targeted_nuclear_phase_alive_dead_tiffs` | `NA` | `NA` | `NA` |
| 266 | 272 | not_checked | no_stored_sha256 | `named_object:run_20260721_163410#selected_classifier_object_rows` | `NA` | `NA` | `NA` |
| 267 | 273 | not_checked | no_stored_sha256 | `raw_source:SUM159/selected_red_cell_nuclear_image_collection` | `NA` | `NA` | `NA` |
| 271 | 277 | not_checked | no_stored_sha256 | `raw_source:SUM159/selected_multimodal_channels_and_masks` | `NA` | `NA` | `NA` |
