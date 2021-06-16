#' Single Cell marker genes Tab UI
#'
#' @export
#' @return None
#'
findMarkerGenesUI <- function(id) {
  ns <- NS(id)


  tagList(fluidRow(column(
    8,
    wellPanel(
      style = "background-color: #fff; border-color: #2c3e50; height: 515px;",
      p(
        tags$b('UpSet plot of intersecting differentially expressed genes (cluster biomarkers)', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
      ),
      hr(),
      plotOutput(ns("UpsetMarker"), height = "360px") %>% withSpinner(color =
                                                                    "#0dc5c1"),
      textOutput(ns("upsetLegend"))
    )
  )),
  fluidRow(column(
    12,
    wellPanel(
      style = "background-color: #fff; border-color: #2c3e50; height: 650px;",
      p(
        tags$b('differentially expressed genes (cluster biomarkers) present per subset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
      ),
      downloadLink(ns('downloadData'), 'Download table as .csv'),
      hr(),
      DT::dataTableOutput(ns("markerTable")) %>% withSpinner(color =
                                                "#0dc5c1")
    )
  )),
  fluidRow(column(
    12,
    wellPanel(
      style = "background-color: #fff; border-color: #2c3e50; height: 650px;",
      p(
        tags$b('differentially expressed genes (cluster biomarkers) summary statistics per subset ', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
      ),
      downloadLink(ns('downloadData2'), 'Download table as .csv'),
      hr(),
      DT::dataTableOutput(ns("sumStatMarkerTable")) %>% withSpinner(color =
                                                       "#0dc5c1")
    )
  )))

}

#' Single Cell marker genes Server
#'
#' @param seuratObjectsList Reactive value containing list of downsampled seurat objects
#' @param selectedNumGenes Reactive value containing user input for number of top genes to test for per cluster per subset
#' @export
#' @return Returns a Reactive value containing a table listing marker genes present in each subset
#'

findMarkerGenes <-
  function(input, output, session, seuratObjectsList, selectedNumGenes) {
    TablesList = list()
    sumStatMarkerTable1 = c()
    #for each seurat object, find marker genes
    for (i in 1:(length(seuratObjectsList))) {
      subsetSize = paste(nrow(seuratObjectsList[[i]]@meta.data) / 1000, "K", sep =
                           "")
      subset.markers = FindAllMarkers(
        seuratObjectsList[[i]],
        only.pos = TRUE,
        min.pct = 0.3,
       # test.use = "MAST",
        logfc.threshold = 0.3
      )
      subset_top_markers <- subset.markers

       if (i == length(seuratObjectsList)) {
        #to only look at top 10 of
        subset_top_markers = subset.markers %>% group_by(cluster) %>% top_n(n = selectedNumGenes, wt = avg_log2FC)


      }

      subset_top_markers$subsetSize <- subsetSize
      #concatenate markers from all subsets together
      sumStatMarkerTable1 = rbind(sumStatMarkerTable1, subset_top_markers)


      #update progress bar, following on from 0.5.
      update_modal_progress((i + 5) / 12)
    }
    sumStatMarkerTable1 = data.frame(sumStatMarkerTable1)
    sumStatMarkerTable1 <-
      sumStatMarkerTable1[, c(7, 1, 2, 3, 4, 5, 6, 8)]


    #parse table for use in upsetplot
    combinedMarkersTable <- tableParserMG(sumStatMarkerTable1)

    #plot upsetPlot via combined table
    upsetPlotMG = upset(
      combinedMarkersTable,
      sets = names(combinedMarkersTable)[ncol(combinedMarkersTable):2],
      nsets = ncol(combinedMarkersTable)-1,
      number.angles = 30,
      point.size = 2.5,
      line.size = 1,
      mainbar.y.label = "Marker Intersections",
      sets.x.label = "Marker genes per subset",
      order.by = "freq",
      keep.order = TRUE,
      text.scale = c(1.5, 1.5, 1.2, 1.2, 1.75, 1.3)
    )
    #remove not needed objects to save mem
   # rm(seuratObjectsList)
    rm(TablesList)

    output$UpsetMarker <- renderPlot({
      upsetPlotMG
    })
    output$upsetLegend <- renderText({
      " UpSet plot showing the number of cell marker genes that are shared
                   across subsets."
    })
    output$sumStatMarkerTable <-
      renderDT({
        sumStatMarkerTable1
      },
      options = list(
        info = F,
        paging = F,
        searching = T,
        stripeClasses = F,
        lengthChange = F,
        scrollY = '445px',
        scrollCollapse = T
      ),
      rownames = F)

    output$markerTable <-
      renderDT({
        combinedMarkersTable
      },
      options = list(
        info = F,
        paging = F,
        searching = T,
        stripeClasses = F,
        lengthChange = F,
        scrollY = '445px',
        scrollCollapse = T
      ),
      rownames = F)
    #Allow users to download tables
    output$downloadData <- downloadHandler(
      filename = function() {
        paste('MarkersTable_', Sys.Date(), '.csv', sep = '')
      },
      content = function(con) {
        write.csv(combinedMarkersTable, con)
      }
    )
    output$downloadData2 <- downloadHandler(
      filename = function() {
        paste('SummaryStatistics_', Sys.Date(), '.csv', sep = '')
      },
      content = function(con) {
        write.csv(sumStatMarkerTable1, con)
      }
    )

    return(combinedMarkersTable)



  }

#' Function that record the presence/absence of markers per subset.
#'
#' @param sumStatMarkerTable1 Reactive value containing table that lists the marker genes present in each subset
#' @export
#' @return Returns a Reactive value listing marker genes resolved per cluster per subset

tableParserMG <- function(sumStatMarkerTable1) {
  #PARSING TABLE FOR USE IN UPSETR**************************************
  #this code record its precesnce/absence of each marker in each of the other subsets
    combinedMarkersTable <- sumStatMarkerTable1

    combinedMarkersTable$Cluster_Marker = paste(combinedMarkersTable$cluster,
                                                combinedMarkersTable$gene)

  combinedMarkersTable <- combinedMarkersTable %>% ungroup() %>%
    select(Cluster_Marker, subsetSize) %>% mutate(presence=1) %>%
    pivot_wider(names_from=subsetSize , values_from=presence) %>%
    replace(is.na(.), 0) %>% as.data.frame()


  #to only look at markers present in reference
  combinedMarkersTable <- combinedMarkersTable %>%
    filter(!(combinedMarkersTable[, ncol(combinedMarkersTable)] == "0"))

  return(combinedMarkersTable)
}
