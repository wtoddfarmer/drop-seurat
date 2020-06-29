library(rmarkdown)

folder <- 'CLASS/'

SOs <- list.files(folder, pattern='*\\.RDS', recursive=FALSE, full.names=TRUE)

for (file in SOs){
  print(SO)
  rmarkdown::render(input = "drop_QC_report.Rmd",
                    output_format = 'html_document',
                    output_file = paste0(gsub(".RDS","", file), "_QC.html"),
                    output_dir = folder, 
                    envir = new.env()
                     
                    
  )
}





