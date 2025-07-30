# Transcritpomic outlier QC 

This pipeline is used to identify outlier samples for eQTL calling using WGCNA. Most approaches for identifying outlier samples tend to rely on using PCA, in which case outliers may only be detected by standardizing first one or two principal components. This may not accurately capture all samples that are outliers in the data 

Using WGCNA  we first calculate a correlation matrix (biweight midcorrelation) across all samples after filtering out lowly expressed genes and then performing a rank normal transformation. From this correlation matrix a network is then derived where nodes are samples and the edges describe the weighted correlation between samples. Finally, using the network outlier samples can be identified by standardizing the connectivity score of each sample. Outlier samples are then identified as a sample have a connectivity Z score < -3. 
