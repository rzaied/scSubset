# scSubset

scSubset: R package to evaluate optimal cell count in single cell transcriptomics

## About 
scSubset is designed to help users identify the sufficient number of cells to use in their scRNA-seq experiments. The package interactively interrogates deposited single cell datasets and down-samples them into smaller subsets having 20%, 40%, 60% and 80% of the parent dataset, respectively. Clustering projections of each subset is compared to that of the reference using the adjusted Rand index (ARI) and normalized mutual information (NMI) scores. The degree of overlap of differentially expressed genes and conserved marker genes between subsets and the reference dataset will also be computed.

![overview figure](https://github.com/rzaied/scSubset/blob/master/figures/overview.png)


## Usage Scenario

10K peripheral mononuclear cells (PBMCs) from a healthy donor obtained from 10X Genomics was used for this example. At a resolution of 0.3, scSubset analysis shows that using 60% of the dataset drops the adjusted rand index (ARI) and normalized mutual information scores (NMI) from 1.0 to to 0.8 and 0.85, respectively:

![ARI/NMI plot](https://github.com/rzaied/scSubset/blob/master/figures/PBMCs_ari_nmi.png)


For each cluster in the full dataset, the top 10 marker genes (MGs) were compared with the marker genes of each subset and an UpSet plot is used to demonstrate the degree of overlap. In this dataset, 70 of 150 MGs were resolved across all subsets. 20 MGs were unique to the reference data set and another 20 were uniquely shared between the 40%, 60%, 80% and the full dataset. 


![UUpSet plot](https://github.com/rzaied/scSubset/blob/master/figures/PBMCs_upset_plot.png)

Considering the slight improvement in results when increasing the dataset from 60% to 80%, we reasoned that 60% coverage would have been sufficient for this dataset and had the capacidt of saving ~£1700 of sequencing costs: 

Number of cells	 | number of clusters | % Overlapping markers	| ARI | NMI | Sequencing cost
------------ | ------------- | ------------- | ------------- | ------------- | -------------
2102 (20%) |	9	60.00 | 0.602 |	0.701	| 840.8
4204 (40%) | 11	73.33 | 0.790 |	0.825 |	1,681.6
6306 (60%) | 12	80.00 |	0.799 |	0.857 |	2,522.4
8408 (80%) | 14	80.00	| 0.919 |	0.916 |	3,363.2
10510 |	15 | 100.00	| 1.000	| 1.000	| 4,204.0

scSubset further provides tables showing the resolved MGs. This allows researchers to also consider subsets that resolve genes of biological interest:

![UUpSet plot](https://github.com/rzaied/scSubset/blob/master/figures/MGstats.png)

#### *To help in the identification of a suitable subset size, the following criteria can be used:*
  1) Subsets resulting in a high marginal increase of the NMI/ARI scores. 
  2) Subsets whose identified conserved marker genes and/or differentially expressed genes have a high degree of overlap with the full dataset.
  3) Subsets that can sufficiently resolve genes of specific biological interest in a given dataset.


## Usage 

Users can upload single or paired datasets; for each dataset, an .h5 file or a matrix.mtx, genes.tsv/features.tsv and barcodes.tsv files will be accepted (output from a Cell Ranger run). The pattern of the mitochondrial genes should be specified e.g. use "^MT-" for human datasets
and "^mt-" for mouse datasets, etc. The desired number of genes to resolve (default is top 10) from the reference dataset should be selected. The resolution should also be selected (default is 0.5). For a single dataset consisting of 10K cells, the analysis could take ~20 minutes. The same is true for an
integrated dataset of 10K cells. Computation of conserved marker genes for integrated datasets is optional and could add around an extra hour to the analysis.



## Dependancies

shiny (>= 1.5.0), shinythemes (>= 1.1.2), Seurat (>= 4.0.0), DT (>= 0.17), shinycssloaders (>= 0.3), shinydashboard (>= 0.7.1), shinyjs (>= 2.0.0), shinybusy (>= 0.2.0), aricode (>= 0.1.2), ggplot2 (>= 3.3.1), reshape2 (>= 1.4.4), UpSetR (>= 1.4.0), hdf5r (>= 1.3.2), tidyverse (>= 1.3.0), metap (>= 1.3), shinyWidgets (>= 0.5.4), cowplot (>= 1.0.0), patchwork (>= 1.0.0), shinyalert (>= 2.0.0), multtest (>= 2.42.0), stringr (>= 1.4.0), MAST (>= 1.16.0)

## Installation
```
library("devtools")
install_github("rzaied/scSubset")

# Run the application
library(scSubset)
scSubsetGo()
```

## Licence
This project is licensed under the MIT License.
