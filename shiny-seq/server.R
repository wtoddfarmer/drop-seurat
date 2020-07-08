library(shiny)
library(ggplot2)
library(Seurat)
source("/Volumes/one/DATA/drop-seq/DropSeuratFunctions.R")
# load data
# load all of the data
# load all tissues

folder <- "/Volumes/one/DATA/drop-seq/"
SOs <- list.files(folder, pattern='*\\_downSample.RDS', recursive=TRUE, full.names=TRUE)
loadDownsampleData(SOs)
#ALL_calc <- readRDS("../ALL/ALL_calc_downSample.RDS")

tissue_list = list("CB", "ENT", "FC", "GP", "HC", "PC", "SN", "STR", "TH")
class_list = list("ASTROCYTE",
                  "ENDOTHELIAL_STALK",
                  "ENDOTHELIAL_TIP",
                  "EPENDYMAL",
                  "MACROPHAGE",
                  "MICROGLIA",
                  "MURAL",
                  "NEUROGENESIS",
                  "NEURON",
                  "OLIGODENDROCYTE",
                  "POLYDENDROCYTE")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  l1 <<- reactive(input$l1)
  l2 <<- reactive(input$l2)
     # Image of cell classes on UMAP
    output$classUmap <- renderImage({
        filename = "class.png"
        list(src = filename,
             contentType = 'image/png',
             width = 300,
             height = 300,
             alt = "UMAP projection showing classes")
    }, deleteFile = FALSE)
    
    # Level 1 Variable
    output$L1 <- renderText({paste0("Level 1: ", input$l1)})
    
    
    # LEVEL 0--------------------------------------------------------------
    # UMAP with Class
    # add if statement for 'none'
    output$classPlot <- renderPlot({
        DimPlot(ALL_calc, group.by = "class", pt.size = 0.25, label = TRUE, repel =TRUE)  + theme(legend.position = "none")
    }, width = 600, height = 400)
    
    # UMAP with Expressiom
    output$umapPlot <- renderPlot({
       FeaturePlot(ALL_calc, features = input$GOI) 
    }, width = 600, height = 400)
    
    # Violin Class
    output$vlnClass <- renderPlot({
        VlnPlot(ALL_calc, features = input$GOI, group.by = "class" ) + ggtitle(NULL) + theme(legend.position = "none")
    }, width = 600, height = 400) 
    
    # Violin Tissue
    output$vlnTissue <- renderPlot({
        VlnPlot(ALL_calc, features = input$GOI, group.by = "tissue" ) + ggtitle(NULL) + theme(legend.position = "none")
    }, width = 600, height = 400) 
   
    
    # LEVEL 1--------------------------------------------------------------
    # Level 1 UMAP plot
    # add conditionals for plots
    
    output$L1overview <- renderPlot({
        if (input$l1 %in% tissue_list){
                DimPlot(get(input$l1), group.by = "class", pt.size = 0.25, label = TRUE, repel =TRUE)  + ggtitle(NULL) + theme(legend.position = "none")}
        else  if (input$l1 %in% class_list){
                DimPlot(get(input$l1), group.by = "tissue", pt.size = 0.25, label = TRUE, repel =TRUE)  + ggtitle(NULL) + theme(legend.position = "none")}    
    }, width = 600, height = 400)     
    
    
    output$umapL1 <- renderPlot({
        FeaturePlot(get(input$l1), features = input$GOI)  + ggtitle(NULL) + theme(legend.position = "none")
    }, width = 600, height = 400)
   
    output$vlnL1 <- renderPlot({
      
      if (input$l1 %in% tissue_list){
        VlnPlot(get(input$l1),
                features = input$GOI, group.by = "class") + ggtitle(NULL) + theme(legend.position = "none")
      }
      else  if (input$l1 %in% class_list){
        VlnPlot(get(input$l1),
                features = input$GOI, group.by = "tissue" ) + ggtitle(NULL) + theme(legend.position = "none")}
    }, width = 600, height = 400)
    
    # LEVEL 2--------------------------------------------------------------
    output$L2 <- renderText({paste0("Level 2: ", input$l2)})
    
    # UMAP of all of the common names for L1 and L2 selections
    output$overviewL2 <- renderPlot({
      
      if (input$l2 %in% tissue_list){
        DimPlot(subset(x = get(input$l1),
                       subset = tissue == l2()),
                group.by = "common_name") +
          theme(legend.position='bottom', legend.text=element_text(size=8), legend.key.size = unit(0.00, 'cm'))
      }
      else  if (input$l2 %in% class_list){
        DimPlot(subset(x = get(input$l1),
                       subset = class == l2()), 
                group.by = "common_name" ) +
          theme(legend.position='bottom', legend.text=element_text(size=8), legend.key.size = unit(0.00, 'cm'))
        }
    }, width = 800, height = 600)
    
    # UMAP of all of the common names for L1 and L2 selections
    output$umapL2 <- renderPlot({
      
      if (input$l2 %in% tissue_list){
        FeaturePlot(subset(x = get(input$l1),
                       subset = tissue == l2()),
                features = input$GOI) + ggtitle(NULL) + theme(legend.position = "none")
      }
      else  if (input$l2 %in% class_list){
        FeaturePlot(subset(x = get(input$l1),
                       subset = class == l2()), 
                features = input$GOI) + ggtitle(NULL) + theme(legend.position = "none")}
    }, width = 800, height = 600)
    
    # output violin where the cells are limited to the L1 and L2 selections
    output$vlnL2 <- renderPlot({
      
      if (input$l2 %in% tissue_list){
        VlnPlot(subset(x = get(input$l1),
                         subset = tissue == l2()),
              features = input$GOI, group.by = "common_name") + 
          ggtitle(NULL) +
          theme(legend.position = "none") 
          }
      else  if (input$l2 %in% class_list){
        VlnPlot(subset(x = get(input$l1),
                           subset = class == l2()), 
                features = input$GOI, group.by = "common_name" ) +
          ggtitle(NULL) +
          theme(legend.position = "none") 
          }
        }, width = 800, height = 800)
    
    # REPORT --------------------------------------------------------------   
    # output report
    output$report <- downloadHandler(
        # For PDF output, change this to "report.pdf"
        filename = "report.html",
        content = function(file) {
            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).
            tempReport <- file.path(tempdir(), "report.Rmd")
            file.copy("report.Rmd", tempReport, overwrite = TRUE)
            
            # Set up parameters to pass to Rmd document
            params <- list(n = input$slider)
            
            # Knit the document, passing in the `params` list, and eval it in a
            # child of the global environment (this isolates the code in the document
            # from the code in this app).
            rmarkdown::render(tempReport, output_file = file,
                              params = params,
                              envir = new.env(parent = globalenv())
            )
        })
    outputOptions(output, 'L1', suspendWhenHidden = FALSE)
})
