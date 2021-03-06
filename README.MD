# Drop-Seurat
This repository contains the code used to analyse the Drop-seq data from the McCarrol lab as well as the Shiny-seq app to query the data. Fork this repo, populate the raw folder with the Drop-seq data, and run the process_drop-seq.Rmd to generate a local version of the Drop-seq data processed with Seurat and viewable with Shiny-seq.

The reports for each of the Seurat object derived from the dataset can be found in this repo here:

[Reports](https://wtoddfarmer.github.io/drop-seurat/)  

Because of their size, the Seurat objects referred to here can be found elsewhere in the lab.




Folders:  
- report_generation contains the scripts for generating basic and QC reports for each object
- raw contains the raw data from Drop-viz  
- ALL contain the calculated object and the report of the cells that have passed QC for all regions
- TISSUE contains the reports of the cells that have passed QC for each region in the dataset  
- CLASS contains the reports of the cells that have passed QC for each class in the dataset   
- the annotated_* folders contain reports all cells with no filtering
- shiny-seq contains the code for the Shiny web app





