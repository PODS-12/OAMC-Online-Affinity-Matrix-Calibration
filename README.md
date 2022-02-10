# Online Affinity Matrix Calibration (OAMC)

Supplementary materials for "Online Affinity Matrix Calibration with Incomplete Observations"

(SIGKDD'2022 submission #730)

The code has been tested on MATLAB R2016b and R2020a. It should be able to run on other recent versions.

## Main files:

- example_main.m: demo of missing data processing in Section 5.3
- example_scale.m: demo of scalable extention in Section 5.4

## Baseline_MVI (missing value imputation):

- impute_zero.m: ZERO Imputation 
- impute_mean.m: MEAN Imputation
- impute_knn.m: kNN Imputation
- impute_lr.m: Linear Regression-based Imputation
- impute_rf.m: Random Forest-based Imputation
- impute_grouse.m: GROUSE Imputation
- impute_kfmc.m: KFMC Imputation

## Baseline_MC (matrix calibration):

- calibrate_dmc.m: DMC Calibration
- calibrate_cmc.m: CMC Calibration
- calibrate_oamc.m: OAMC Calibration
- calibrate_oamc_block.m: Block OAMC Calibration
- calibrate_soamc.m: Scalable OAMC Calibration

## Other files:

- demo_data.mat: data file for demo_missing.m
- similarity.m: approximate similarity matrix with incomplete samples


Feb. 9th, 2022
