# R Functions for Drop-Seq and Seurat

# Convert the data provided from Drop-Viz website into annotated Seurat object
# Function to generate R annotated Seurat object from the drop-seq (Saunders et al) sparse matix and annotation files. 
# Note that this function will not calculate any transformations or dimensional reductions.
extractDGE <- function(rawDGE, outcomes, idents, prefix){
  dataDir <- dirname(rawDGE)
  baseName <- strapplyc(as.character(rawDGE), "P60(.*).raw", simplify = TRUE)
  
  RDS <- readRDS(outcomes)
  RDS[] <- lapply(RDS, as.character)
  RDS$tissue_subcluster <- paste0(prefix, "_", RDS$subcluster)
  RDS$names <- rownames(RDS)
  RDS <- merge(RDS, idents, by.x="tissue_subcluster", by.y="tissue_subcluster", all.x = TRUE)
  rownames(RDS) <- RDS$names
  DGE <- loadSparseDge(rawDGE)
  SO <- CreateSeuratObject(counts = DGE, 
                           min.cells = 3,
                           min.features = 200,
                           project = "drop-seq")
  SO <- PercentageFeatureSet(object = SO, pattern = "^mt-", col.name = "percent.mito")
  SO <- AddMetaData(object = SO, metadata = prefix, col.name = "tissue")
  SO <- AddMetaData(object = SO, metadata = select(RDS, subcluster.x), col.name = "subcluster")
  SO <- AddMetaData(object = SO, metadata = select(RDS, tissue_subcluster), col.name = "tissue_subcluster")
  SO <- AddMetaData(object = SO, metadata = select(RDS, reason), col.name = "reason")
  SO <- AddMetaData(object = SO, metadata = select(RDS, class), col.name = "class")
  SO <- AddMetaData(object = SO, metadata = select(RDS, full_name), col.name = "full_name")
  SO <- AddMetaData(object = SO, metadata = select(RDS, common_name), col.name = "common_name")
  SO@meta.data$class[is.na(SO@meta.data$class)] <- "UNKNOWN"
  SO@meta.data$common_name[is.na(SO@meta.data$common_name)] <- "unknown"
  SO@meta.data$full_name[is.na(SO@meta.data$full_name)] <- "unknown"
  SO@meta.data$reason[is.na(SO@meta.data$reason)] <- "passed"
  print(paste0("Saving annotated Seurat object for ", prefix))
  saveRDS(SO, file=paste0("annotated/", prefix,"_annotated.RDS"))
  return(SO)
  
}
# Filter a Seurat object on nFeatures and "passed" in reason metadata field
dropPass <- function(SO){
  filteredSO <- subset(x = SO, subset = nFeature_RNA > 500)
  filteredSO <- subset(x = filteredSO, subset = reason == "passed")
  return(filteredSO)
}



# Perform standard calculations on a Seurat object
calcDrop <- function(SO){
  # SO <- dropPass(SO)
  #stopifnot(length(SO@meta.data$orig.ident) >100)
  message("Performing Normalization and Variance Stablization")
  SO <- suppressWarnings(SCTransform(SO,
                                     vars.to.regress = "percent.mito" ,
                                     verbose = TRUE,
                                     return.only.var.genes = FALSE
                                      ))
  
  message("Running PCA")

  SO <- try(suppressWarnings(RunPCA(SO, verbose = FALSE)))
  # write test to extract max number of PCs
  # set dims 
  message("Testing Dimensionality")
  try(if (dim(SO@reductions[["pca"]])[2] < 30){
    dims=dim(SO@reductions[["pca"]])[2]
  }
  else{dims=30})
  # SO <- JackStraw(SO, num.replicate = 100, dims = 50)
  # SO <- ScoreJackStraw(SO, dims = 1:50)
  message("Running Dimension Reduction")
  SO <- try(RunUMAP(SO, dims = 1:30, verbose = FALSE))
  SO <- try(RunTSNE(SO, dims = 1:30, verbose = FALSE))
  # clustering
  message("Finding Nearest Neighbors")
  SO <- try(FindNeighbors(SO, dims = 1:30, verbose = FALSE))
  message("Finding Clusters")
  SO <- try(FindClusters(SO, verbose = FALSE))
  #saveRDS(SO, file=paste0("calculated/", prefix,"_ann.RDS")) #
  return(SO)
}

# function to load downsampled Seurat objects
loadDownsampleData <- function (file_list){
  
  for (file in file_list){
    print(paste0("loading ", file, "..."))
    var_name <- strsplit(basename(file), "_downSample")[[1]][1]
    print(paste0("Assigning as ", var_name))
    assign(var_name, readRDS(file), envir = .GlobalEnv)
  }
  
}


