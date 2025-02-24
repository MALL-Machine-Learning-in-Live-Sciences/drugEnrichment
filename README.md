## Drug enrichment based on Gene Regulatory Networks
These scripts are designed to study Gene Regulatory Networks (GRNs) associated with drug resistance. They include the inference of transcription factors and metabolic pathways, and their relationship with drug response across various cancer cell lines.


### How to download data

- Data come from DepMap repository. All matrices can be downloaded from this [link](https://depmap.org/portal/data_page/?tab=allData)
    - DepMAp 23Q2 Data:
    - Model.csv
    - OmicsExpressionProteinCodingGenesTPMLogp1.csv
    - OmicsSomaticMutations.csv
- PRISM 19Q4 Data:
    - secondary-screen-dose-response-curve-parameters.csv
- Sanger GDSC1 and GDSC2: 
    - sanger-dose-response.csv
    - screened_compounds_rel_8.5.csv
- GSE162104 avalible at [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162104)

### How to use

 - Run the scripts sequentially, ensuring that the paths at the beginning of each script are updated accordingly.
     - 01: Preprocessing the 3 ds
     - 02: Transcriptión factor activity inference and Drug Enrichment
     - 03: Non-negative matrix factorization and Drug Enrichment
     - 04-06: Cohort study
 - It is possible to select which tumor type to analyze; the list is available at the end of the scripts.
 - The code for the plots is located in the plots folder. Ensure that all intermediate data is available.
 - The code is designed to plot specific MoA/biological features. It is important to check which features are being used for the plots.
