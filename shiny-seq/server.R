library(shiny)
library(ggplot2)
library(Seurat)

# load data
SO <- readRDS("GP_downsample.RDS")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    # Image of cell classes on UMAP
    output$classUmap <- renderImage({
        filename = "class.png"
        list(src = filename,
             contentType = 'image/png',
             width = 300,
             height = 300,
             alt = "UMAP projection showing classes")
    }, deleteFile = FALSE)

    # UMAP
    output$umapPlot <- renderPlot({
       FeaturePlot(SO, features = input$GOI) 
    }, width = 600, height = 400)
    
    # Violin
    output$vlnPlot <- renderPlot({
        VlnPlot(SO, features = input$GOI, group.by = "class", ) + ggtitle(NULL) + theme(legend.position = "none")
    }, width = 600, height = 400) 
    
    
        
          
})
