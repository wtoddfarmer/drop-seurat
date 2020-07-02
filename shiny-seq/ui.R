
library(shiny)

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

shinyUI(fluidPage(

    # Application title
    titlePanel("Shiny-seq: Single Cell Gene Expression Report"),
    tags$h4("This app shows the expression of a single gene. The data used is from the publication", 
            tags$a(href="https://www.cell.com/cell/fulltext/S0092-8674(18)30955-3?", "Saunders et. al, 2018"),
            "and was processed using the R Package",
            tags$a(href="https://satijalab.org/seurat/", "Seurat.")),
    
    # Sidebar with gene entry field
    sidebarLayout(
        sidebarPanel(
            tags$h1("Gene of Interest"),
            textInput("GOI", "Input a mouse gene symbol without quotes (i.e. Slc1a3) to see which cell types express it :", value = "none", width = NULL,
                      placeholder = NULL),
            selectInput("l1", "Select Level 1:",
                        list('none',
                            `Tissue` = tissue_list,
                             `Class` = class_list)
            ),
            selectInput("l2", "Select Level 2:",
                        list('none',
                            `Tissue` = tissue_list,
                             `Class` = class_list)
            ),
           # imageOutput("classUmap"),
            
            #downloadButton("report", "Generate report")
            
        ),

        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Level 0", 
                                 tags$h1("Level 0: All Cells/Regions/Classes"),
                                 tags$h2("Cell Class on UMAP"),
                                 plotOutput("classPlot"),
                                 tags$h2("Expression of GOI on UMAP"),
                                 plotOutput("umapPlot"),
                                 tags$h2("Violin Plot GOI in Cell Classes"),
                                 plotOutput("vlnClass"),
                                 tags$h2("Violin Plot GOI by Tissue"),
                                 plotOutput("vlnTissue"),
                                 ),
                        tabPanel("Level 1", 
                                 textOutput("L1"),
                                 tags$head(tags$style("#L1{font-size: 36px; font-weight: 500; line-height: 1.1}")),
                                 plotOutput("L1overview"),
                                 tags$h2("Expression of GOI on L1 UMAP"),
                                 plotOutput("umapL1")),
                        tabPanel("Level 2", 
                                 textOutput("L2"),
                                 tags$head(tags$style("#L2{font-size: 36px; font-weight: 500; line-height: 1.1}")),
                                 #plotOutput("L2overview"),
                                 tags$h2("Expression of GOI on L2 UMAP"),
                                 #plotOutput("umapL1")),
                                 )
            )
            
        
            
            #includeHTML("../ALL/ALL_calc.html")
        )
    )
))
