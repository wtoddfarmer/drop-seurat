# This script will create and save a downsampled vesion of all Seurat objects in the path.

# Note: non Seurat objects will cause the script to fail


folders <- c("CLASS", "TISSUE", "ALL")

SOs <- list.files(folders, pattern='*\\.RDS', recursive=TRUE, full.names=TRUE)

for (file in SOs){
  print(paste0("Reading ", file, "..."))
  SO <- readRDS(file)
  print(SO)
  
  print("Downsampling...")
  # need a reasnable value for downsample (max number of cells/ident...needs to vary some 
  SO <- SetIdent(SO, value = SO@meta.data$common_name)
  subset <- subset(SO, cells = WhichCells(object = SO, downsample = 1000))
  print("created:")
  print(subset)
  name_agg <- subset@meta.data %>% group_by(common_name) %>% count(common_name) 
  print(name_agg)
  print("Writing to disk at:")
  new_file_name <- gsub(".RDS", "_downSample.RDS", file)
  print(new_file_name)
  saveRDS(subset, file= new_file_name)       
  
}
