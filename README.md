# NucleosomeProteomeAnalysis
This was tested using R version 4.2.2

## To install dependencies, 
#### run `installation.R` in RStudio

## Step 1: run `preproccess.R` with 
```
RScript preprocess.R [experiments options: grouping_all | v1 | walk | ap] [/path/to/output/dir]
```
## Step 2: run analysis script for profile plot outputs, PCA, and other protein summarization analysis
### Run `analysis.R`for multiple mixture experiment
```
RScript analysis.R [/path/to/output/dir/with/preprocess/outputs]
```

### Run `analysis_single_run.R`for single-run experiment
```
RScript analysis_single_run.R [/path/to/output/dir/with/preprocess/outputs]
```
## Step 3: Pairwise analysis
### Run `change.R` to output volcano plots for any experiment
```
RScript change.R [experiments options: grouping_all | v1 | walk | ap] [/path/to/output/dir/with/preprocess/outputs]
```
## Step 4: Heatmap analysis
### Run `heatmap.R` to output heatmap plots for any experiment
```
RScript heatmap.R [experiments options: grouping_all | v1 | walk | ap] [/path/to/output/dir/with/preprocess/outputs]
```
