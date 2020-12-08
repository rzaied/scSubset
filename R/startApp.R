scSubsetGo <- function() {

  if (interactive()){
    #reset tmp dir
    unixtools::set.tempdir("/data/Roan/tmp/")
    #load packages
    #for app.R
    library(shiny)
    library(shinythemes)
    library(unixtools)
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
    #set tmp directory
    library(unixtools)

    #for exampleDatar.R
    #to plot bar plot
    library(UpSetR)

    #for integration.R
    #for NMI/ARI
    library(aricode)
    #to plot bar plot
    library(ggplot2)
    #for melt function
    library(reshape2)


    #for filter() function
    library(dplyr)
    library(SeuratData)
    library(shinyalert)
    library(unixtools)
    library(shinyalert)

    #BiocManager::install("SeuratData")
    #devtools::install_github("satijalab/seurat-data", ref = 'develop')
    library(cowplot)
    library(patchwork)
   # ulimit::memory_limit(80000)
    options(shiny.maxRequestSize=500*1024^2)
    app <- shinyApp(ui = shinyUI(ui),
                    server = shinyServer(server))
    runApp(app)

   # runApp("R")
  }
}


