library(rmarkdown)

folder <- 'CLASS'

SOs <- list.files(folder, pattern='*\\.RDS', recursive=FALSE, full.names=TRUE)

for (file in SOs){
  print(file)
  rmarkdown::render(input = "report_generation/drop_report_class.Rmd",
                    output_format = 'html_document',
                    output_file = paste0(gsub(".RDS","", file), ".html"),
                    output_dir = folder, 
                    envir = new.env()
                     
                    
  )
}




