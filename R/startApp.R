scSubsetGo <- function() {

  if (interactive()){
    #load packages
    #for app.R
    library(shiny)
    library(shinythemes)
    library(Seurat)
    library(DT)
    library(shinycssloaders)
    library(shinydashboard)
    library(shinyjs)
    library(shinybusy)

    #for autoClustering.R
    library(aricode)
    library(ggplot2)
    library(reshape2)
    library(dplyr)

    #for exampleDatar.R
    #to plot bar plot
    library(UpSetR)

    #for integration.R
    library(metap)
    #for NMI/ARI
    library(aricode)
    #to plot bar plot
    library(ggplot2)
    #for melt function
    library(reshape2)


    #for filter() function
    library(dplyr)
    library(shinyalert)
    library(cowplot)
    library(patchwork)

    options(shiny.maxRequestSize=500*1024^2)
    app <- shinyApp(ui = shinyUI(ui),
                    server = shinyServer(server))
    runApp(app)

   # runApp("R")
  }
}


