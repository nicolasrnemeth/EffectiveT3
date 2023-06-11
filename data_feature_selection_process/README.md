This folder contains saved numpy arrays containing the balanced_accuracy_scores
for 6 different models as trained on 6 different training set splits -> bacc_scores_d{num_train_split}.npy

feat_names.npy -> feature names

optimized_parameters_d{num_train_split}.json -> optimized hyperparameters as obtained by one-by-one and heuristic optimization

shap_value_mean_d1Tod6.npy -> mean of the shap-value obtained from the 6 different test splits

model_d{num_train_split}.bin -> the trained models that can be loaded into memory

sorted_idxs_highToLowImp.npy -> sorted indexes of the shap-values

selected_features_85dim_obtained_by_gbdt_with_regularization.npy -> feature names of the features that have been selected by the feature selection process

´protein_sequences_train_test_splits_1to6´ -> this folder contains the training and test splits of the protein sequences
that were used to obtain the data in ´data_feature_selection_process´
