# Drop-Seurat
This repository contains the code used to analyse the Drop-seq data from the McCarrol lab.  

Becuase of their size, the Seurat objects referred to here can be found elsewhere in the lab.

Folders:

- raw contains the raw data from Drop-viz  
- ALL contain the calculated object and the report of the cells that have passed QC for all regions
- TISSUE contains the reports of the cells that have passed QC for each region in the dataset  
- CLASS contains the reports of the cells that have passed QC for each class in the dataset   
- the annotated_* folders contain reports all cells with no filtering
- shiny-seq contains the code for the Shiny web app



TODO:

make main data ALL
  make downsampled versions of all seurat objects

select level 1 region or class

dropdown reflects l1 choice

select level 2 region, class, or none

dropdown reflects l2 choice

button to create report


