# Single cell sequencing of radial glia progeny reveals diversity of newborn neurons in the adult zebrafish brain

## Description of files in repository

Contains the source code for data analysis performed for Lange et al., 2019. The sequencing data have been deposited in NCBI's Gene Expression Omnibus and are accessible through GEO Series accession number GSE137525 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137525). For a description of the computational methods please see Lange et al., 2019.

* `envrionement.yml`: anaconda environment for executing `zebrafish_neurogenesis_smartseq_prepare_data.R` and `zebrafish_neurogenesis_smartseq_manuscript.ipynb`
* `zebrafish_neurogenesis_smartseq_prepare_data.R`: reading the count data, calculation of quality metrics
* `zebrafish_neurogenesis_smartseq_manuscript.ipynb`: contains all the analysis of the single-cell seq data except for the cell type homology analysis
* `environment_2.yml`: anaconda environment for executing `Hochgerner_2018.ipynb`
* `mart_export.txt`: the list of orthogolous genes used for the cell type homology analysis
* `Hochgerner_2018.ipynb`: contains the cell type homology analysis

## Citation

Lange C, Rost F, Machate A, Reinhardt S, Lesche M, Weber A, Kuscha V, Dahl A, Rulands S, Brand M. Single cell sequencing of radial glia progeny reveals diversity of newborn neurons in the adult zebrafish brain, under review
