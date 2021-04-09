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
        tags$b('UpSet plot of intersecting marker genes', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
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
        tags$b('Marker genes present per subset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
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
        tags$b('Marker genes summary statistics per subset ', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
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
    #carrying on from Autoclustering_V2.R script,
  #  combinedMarkersTableHeader = c("Cluster_Marker")
    TablesList = list()
    sumStatMarkerTable1 = c()
    #for each seurat object, find marker genes
    for (i in 1:(length(seuratObjectsList))) {
      subsetSize = paste(nrow(seuratObjectsList[[i]]@meta.data) / 1000, "K", sep =
                           "")
      print(subsetSize)
      #logfc.threshold is by default 0.25, increasing to 0.3 speeds up the function but could miss weaker signals
      subset.markers = FindAllMarkers(
        seuratObjectsList[[i]],
        only.pos = TRUE,
        min.pct = 0.3,
        logfc.threshold = 0.3
      )
     # change avg_logFC to avg_log2FC
      markers_col_names<- names(subset.markers)
      subset_top5 = subset.markers %>% group_by(cluster) %>% top_n(n = selectedNumGenes, wt = avg_log2FC)
  #  x<-  subset.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_logFC)
      #to create variables dynamically and assign them with respective
      #top 10 table "subset_2000top10

      #this is just so that th e tables dont get overwritten and are saved in env
      assign(paste("subset_", subsetSize, "_top5", sep = ""), subset_top5)

      #append to list
     # TablesList[[i]] = subset_top5
      #add subsetsize coloumn to table
      subset_top5$subsetSize <- subsetSize
      #concatenate markers from all subsets together
      sumStatMarkerTable1 = rbind(sumStatMarkerTable1, subset_top5)

   #  combinedMarkersTableHeader = c(combinedMarkersTableHeader, subsetSize)

      #update progress bar, following on from 0.5.
      update_modal_progress((i + 5) / 12)
    }
    sumStatMarkerTable1 = data.frame(sumStatMarkerTable1)
    sumStatMarkerTable1 <-
      sumStatMarkerTable1[, c(7, 1, 2, 3, 4, 5, 6, 8)]
    #rownames(sumStatMarkerTable1)<-c(1:nrow(sumStatMarkerTable1))


    #parse table for use in upsetplot
    combinedMarkersTable <- tableParserMG(sumStatMarkerTable1)


    #finally, rename the coloumns in combinedMarkersTable
  #  names(combinedMarkersTable) = combinedMarkersTableHeader
    print("line 100 findMG")
    #change factor to numeric
    # combinedMarkersTable[, 2:(length(TablesList) + 1)] <-
    #   sapply(combinedMarkersTable[, 2:(length(TablesList) + 1)], as.character)
    # combinedMarkersTable[, 2:(length(TablesList) + 1)] <-
    #   sapply(combinedMarkersTable[, 2:(length(TablesList) + 1)], as.numeric)
    print("line 106 findMG")
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

    print("line 121 findMG")
    #remove not needed objects to save mem
    rm(seuratObjectsList)
    rm(TablesList)

    output$UpsetMarker <- renderPlot({
      upsetPlotMG
    })
    print("line 130 findMG")
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

    print("line 150 findMG")
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
    print("line 162 findMG")
    #Allow users to download tables
    output$downloadData <- downloadHandler(
      filename = function() {
        paste('MarkersTable_', Sys.Date(), '.csv', sep = '')
      },
      content = function(con) {
        write.csv(combinedMarkersTable, con)
      }
    )
    print("line 172 findMG")
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

  print("line 256 findMG")
  return(combinedMarkersTable)
}
