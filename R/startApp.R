#' Function to load required dependanies and start Shiny session
#'
#' @export

scSubsetGo <- function() {

  if (interactive()){
    #load packages

    library(shiny)
    library(shinythemes)
    library(Seurat)
    library(DT)
    library(shinycssloaders)
    library(shinydashboard)
    library(shinyjs)
    library(shinybusy)
    library(aricode)
    library(ggplot2)
    library(reshape2)
    library(tidyverse)
    library(UpSetR)
    library(metap)
    library(aricode)
    library(ggplot2)
    library(reshape2)
    library(shinyalert)
    library(cowplot)
    library(patchwork)
    library(stringr)
    library(MAST)

    options(shiny.maxRequestSize=500*1024^2)
    app <- shinyApp(ui = shinyUI(ui),
                    server = shinyServer(server))
    runApp(app)

  }
}


