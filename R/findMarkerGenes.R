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
    combinedMarkersTableHeader = c("Cluster_Marker")
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
      #change avg_logFC to avg_log2FC
      subset_top5 = subset.markers %>% group_by(cluster) %>% top_n(n = selectedNumGenes, wt = avg_log2FC)

      #to create variables dynamically and assign them with respective
      #top 10 table "subset_2000top10

      #this is just so that the tables dont get overwritten and are saved in env
      assign(paste("subset_", subsetSize, "_top5", sep = ""), subset_top5)

      #append to list
      TablesList[[i]] = subset_top5
      #add subsetsize coloumn to table
      subset_top5$subsetSize <- subsetSize
      #concatenate markers from all subsets together
      sumStatMarkerTable1 = rbind(sumStatMarkerTable1, subset_top5)

      combinedMarkersTableHeader = c(combinedMarkersTableHeader, subsetSize)

      #update progress bar, following on from 0.5.
      update_modal_progress((i + 5) / 12)
    }
    sumStatMarkerTable1 = data.frame(sumStatMarkerTable1)
    sumStatMarkerTable1 <-
      sumStatMarkerTable1[, c(7, 1, 2, 3, 4, 5, 6, 8)]
    #rownames(sumStatMarkerTable1)<-c(1:nrow(sumStatMarkerTable1))


    #parse table for use in upsetplot
    combinedMarkersTable <- tableParserMG(TablesList)

    combinedMarkersTable = data.frame(combinedMarkersTable)
    #finally, rename the coloumns in combinedMarkersTable
    names(combinedMarkersTable) = combinedMarkersTableHeader
    print("line 100 findMG")
    #change factor to numeric
    combinedMarkersTable[, 2:(length(TablesList) + 1)] <-
      sapply(combinedMarkersTable[, 2:(length(TablesList) + 1)], as.character)
    combinedMarkersTable[, 2:(length(TablesList) + 1)] <-
      sapply(combinedMarkersTable[, 2:(length(TablesList) + 1)], as.numeric)
    print("line 106 findMG")
    #plot upsetPlot via combined table
    upsetPlotMG = upset(
      combinedMarkersTable,
      sets = combinedMarkersTableHeader[6:2],
      nsets = 5,
      number.angles = 30,
      point.size = 2.5,
      line.size = 1,
      mainbar.y.label = "Marker Intersections",
      sets.x.label = "Marker genes per subset",
      order.by = "freq",
      keep.order = TRUE,
      text.scale = c(1.5, 1.5, 1, 1.2, 1.75, 1.3)
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

#' Function that iterates through cell markers and record their presence/absence per subset.
#'
#' @param TablesList Reactive value containing a list of tables that list the marker genes present prer cluster per
#'   subset
#' @export
#' @return Returns a Reactive value listing marker genes resolved per cluster per subset


tableParserMG <- function(TablesList) {
  #table will be used to plot upsetplot
  combinedMarkersTable = c()
  #row that temporarily stores info about each marker gene from each subset
  #if a marker is present, the subset will have a score of 1. else, it will score a 0
  tmpRow = c()
  #PARSING TABLE FOR USE IN UPSETR**************************************
  #this code iterates through each cell marker and record its precesnce/absence in each of
  #the other subsets.

  #for each table, add the 8th coloumn (STARTING FROM smallest subset TABLE)
  for (i in 1:length(TablesList)) {
    #1) make a cluster_marker column in each table
    TablesList[[i]]$Cluster_Marker = paste(TablesList[[i]]$cluster,
                                           TablesList[[i]]$gene)
  }

  #for each table, starting at smallest subset
  for (i in 1:length(TablesList)) {
    #update progress bar, following on from 0.8.
    update_modal_progress((i + 8) / 13.5)

    #to know progress
    print(paste("parsing table:", i, "of 5"))

    #2) iterate through each row in that table
    for (row in 1:nrow(TablesList[[i]])) {
      #check if marker is already present in the combined table
      #col 8 holds the cluster and marker name
      clusterMarker = paste(TablesList[[i]][row, 8])
      coloumnWithPattern = which(grepl(clusterMarker, combinedMarkersTable))

      #to know progress thus far
      print(paste("row #", row, " of ", nrow(TablesList[[i]]), sep = ""))

      #if marker isn't present in combinedMarkersTable
      if (length(coloumnWithPattern) == 0) {
        #1)add it to tmpROW table
        #clusterMarker=paste(MyList[[i]][row,8])
        tmpRow = rbind(tmpRow, clusterMarker)
        print("line 227 findMG")
        for (j in 1:length(TablesList)) {
          #search if that marker exist in myTables
          x = which(grepl(clusterMarker, TablesList[[j]]))
          print("line 231, findMG ", i)
          #if it isn't present,
          if (length(x) == 0) {
            #add 0 to tmpRow
            tmpRow = cbind(tmpRow, 0)
            print("line 236 findMG")
          }
          else
            #add1 to tmp row
            tmpRow = cbind(tmpRow, 1)
          print("line 241 findMG")

        } #then search marker in next table
        #once done crosscomparing for selected marker,
        #add tmpRow to the combinedMarkersTable
        combinedMarkersTable = rbind(combinedMarkersTable, tmpRow)
        print("line 247 findMG")
        #reinitiate tmpRow to 0
        tmpRow = c()
      } #if it is present already it means you've already cross-compared it
      #so move on to next marker
    }
    #once done will ll rows of a table, move on to next table

  }
  print("line 256 findMG")
  return(combinedMarkersTable)


}

