# NucleosomeProteomeAnalysis

## First run `preproccess.R` with 
```
RScript preprocess.R [experiments options: grouping_all | v1 | walk | ap] [/path/to/output/dir]
```

## Run `analysis.R` for profile plot outputs, PCA, and other protein summarization analysis for multiple mixture experiment
```
RScript analysis.R [/path/to/output/dir/with/preprocess/outputs]
```

## Run `analysis.R` for profile plot outputs, PCA, and other protein summarization analysis for multiple mixture experiment
```
RScript analysis_single_run.R [/path/to/output/dir/with/preprocess/outputs]
```

## Run `heatmap.R` to output heatmap plot for multiple mixture experiment
```
RScript heatmap.R [/path/to/output/dir/with/preprocess/outputs]
```
