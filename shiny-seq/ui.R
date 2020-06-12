
library(shiny)

shinyUI(fluidPage(

    # Application title
    titlePanel("Shiny-seq: Single Cell Gene Expression Report"),
    tags$h4("This app shows the expression of a single gene in the mouse globus palidus. The data used is from the publication", 
            tags$a(href="https://www.cell.com/cell/fulltext/S0092-8674(18)30955-3?", "Saunders et. al, 2018"),
            "and was processed using the R Package",
            tags$a(href="https://satijalab.org/seurat/", "Seurat.")),
    tags$h2("Plots shown:"),
    tags$h4("UMAP: the level of expression on the dimensional reduction plot."),
    tags$h4("Violin Plot: the level of expression in each cell class."),
    # Sidebar with gene entry field
    sidebarLayout(
        sidebarPanel(
            tags$h1("Gene of Interest"),
            textInput("GOI", "Input a mouse gene symbol without quotes (i.e. Slc1a3) to see which cell types express it :", value = "Slc1a3", width = NULL,
                      placeholder = NULL),
            tags$h4("Some Ideas:"),
            tags$h6("Slc1a2 is expressed in astrocytes"),
            tags$h6("Slc17a6 is expressed in some excitory neurons"),
            tags$h6("Gad1 is expressed in some inhibitory neurons"),
            tags$h6("Olig2 is expressed in oligodendrocytes"),
            tags$h6("Tek is expressed in endotheilal cells"),
            tags$h1("Cell Classes on UMAP"),
            imageOutput("classUmap")
              
        ),

        mainPanel(
            # Show a plots
            tags$h1("Expression of GOI on UMAP"),
            
            plotOutput("umapPlot"),
            tags$h1("Violin Plot GOI in Cell Classes"),
            plotOutput("vlnPlot")
        )
    )
))
