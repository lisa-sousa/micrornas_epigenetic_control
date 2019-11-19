# A statistical model for epigenetic control of miRNAs

Published in German Conference on Bioinformatics 2016 Collection
doi: 10.7287/peerj.preprints.2423v1

### Abstract
microRNAs are small, non-coding RNAs involved in post-transcriptional gene regulation. Since the dysregulation of only a few miRNAs can affect many biological pathways, miRNAs are thought to play a key role in cancer development and can be used as biomarkers for cancer diagnosis and prognosis. In order to understand how miRNA dysregulation leads to a cancer phenotype it is important to determine the basic regulatory mechanisms that drive miRNA expression. Although much is known about miRNA-mediated post-transcriptional regulation, little is known about the epigenetic control of miRNAs. Here, we performed cell-line specific miRNA promoter predictions and built a classification model for expressed and non-expressed miRNAs based on several epigenetic features, e.g. histone marks and DNA methylation at both, miRNA promoters and around miRNA hairpins. We were able to classify intragenic and intergenic miRNAs with an accuracy of 79% and 85%, respectively, and identified the most important features for classification via feature selection. Surprisingly, we found that DNA methylation seems to have a dual role in regulating miRNA expression at transcriptional level: at promoters, high levels of DNA methylation correlate with transcriptional repression, while around miRNA hairpins high levels of DNA methylation have a positive impact on the expression level of the mature miRNA.

###Content

This repository contains all data sources and scripts that were used to perform all analysis for the paper. 

data: the folder contains the called miRNA promoters (with PROmiRNA) for expressed and non-expressed miRNAs in two cell lines: IMR90 and HELA. Furthermore, it contains the expression status (RPKM values), calculated from small RNA-seq for all miRNAs in both cell lines. In addition, it includes processed data for DNA methylation and histone marks on both cell lines.

scripts: this folder contains all scripts needed to calculate expression status of miRNAs as well as to preprocess DNA methylation and histone marks data sets. Furthermore it contains the scripts for the logistic regression model with additional regularization.