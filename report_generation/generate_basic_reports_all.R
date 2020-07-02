library(rmarkdown)

folder <-"ALL"

SOs <- list.files(folder, pattern='ALL_calc.RDS', recursive=FALSE, full.names=TRUE)
print(SOs)
for (file in SOs){
  print(file)
  rmarkdown::render(input = "report_generation/drop_report_all.Rmd",
                    output_format = 'html_document',
                    output_file = paste0(gsub(".RDS","", file), ".html"),
                    output_dir = folder, 
                    envir = new.env()
                     
                    
  )
}




