conservedMarkersUI <- function(id) {
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
      DTOutput(ns("markerTable")) %>% withSpinner(color =
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
      DTOutput(ns("sumStatMarkerTable")) %>% withSpinner(color =
                                                       "#0dc5c1")
    )
  )))
}




findConservedMarkers <-
  function (input, output, session, seuratObjectsList, selectedNumGenes) {
    #carrying on from findDEgenes_integration.R script,
    combinedMarkerGenesTableHeader = c("Cluster_Marker")
    #list to store conserve marker genes from each cluster
    conservedMarkers_Tables_List = list()
    #combinedtop conserved markers table holds marker genes for one subset at a time
    combinedMarkers_top = c()
    #table holds all marker genes of all subsets
    sumStatMarkerTable1 = c()
    #for each seurat object, add data cloumns
    for (i in 1:(length(seuratObjectsList))) {
      DefaultAssay(seuratObjectsList[[i]]) <- "RNA"
      subsetSize = paste(nrow(seuratObjectsList[[i]]@meta.data) / 1000, "K", sep =
                           "")
      print(subsetSize)
      #update progress bar, following on from 0.5
      update_modal_progress((i + 5) / 10.5)
      #list to easily acess the clusters in each subset
      levelsList = levels(seuratObjectsList[[i]])
      #for each cluster in list
      for (j in 1:length(levelsList)) {
        #logfc.threshold is by default 0.25, increasing to 0.4 speeds up the function but could miss weaker signals
        #Find conserved marker genes (marker genes shared b/w both conditions)
        try(conserved.markers <-
              FindConservedMarkers(
                seuratObjectsList[[i]],
                ident.1 = levelsList[[j]],
                grouping.var = "orig.ident",
                verbose = TRUE,
                min.pct = 0.7,
                logfc.threshold = 0.3
              ))
        print("line 81")

        #in case a cluster doesn't exist in a condition, only one condition will be used (i.e. only 5 cols)

        #if tested cluster is present in both conditions (i.e. there are 12 cols) carry on, otherwise skip that cluster
        if (ncol(conserved.markers) == 12)  {

          #only keep markers that are positivly expressed in both conditions
          top_conservedMarkers = conserved.markers[conserved.markers[, 2] >= 0 &
                                                     conserved.markers[, 7] >= 0,]

          print("line 85")
          #use minimump_p_value, which is a combinded p_val to only keep significant markers
          top_conservedMarkers = top_conservedMarkers[top_conservedMarkers$minimump_p_val <=
                                                        0.01,]
          print("line 89")
          #select top5 (using many genes significantly increases computation time)
          top_conservedMarkers = top_conservedMarkers %>% top_n(n = -selectedNumGenes, wt = minimump_p_val)
          #add a coloumn with the identity of the cluster (in try block incase no top consv markers )
          try(top_conservedMarkers$cluster <-
                levelsList[[j]],
              silent = T)
          print("line 96")
          #append those to a combined_top table
          combinedMarkers_top = rbind(combinedMarkers_top, top_conservedMarkers)
          print("line 99")
        }


      }
      combinedMarkers_top = data.frame(combinedMarkers_top)

      #append tables to list
      conservedMarkers_Tables_List[[i]] = combinedMarkers_top
      #to create variables dynamically and assign them with respective
      #top 10 table "subset_2000top10
      # this is just so that the tables dont get overwritten and are saved in env for me
      assign(paste("top_marker_", subsetSize, sep = ""),
             combinedMarkers_top)
      #add subsetsize col
      combinedMarkers_top$subsetSize <- subsetSize
      combinedMarkers_top$gene = rownames(combinedMarkers_top)
      sumStatMarkerTable1 = rbind(sumStatMarkerTable1, combinedMarkers_top)

      #once we compiled the top marker genes from each cluster of a subset, move on to next subset
      combinedMarkers_top = c()
      #rename column to substsize
      combinedMarkerGenesTableHeader = c(combinedMarkerGenesTableHeader, subsetSize)

    }

    #add numbered row names to table
    rownames(sumStatMarkerTable1) <- c(1:nrow(sumStatMarkerTable1))
    combinedMarkersTable <- tableParserCMG(conservedMarkers_Tables_List, combinedMarkerGenesTableHeader)
    print("line 135 findCMG")
    #plot upsetPlot via combined table
    upsetPlotMG = upset(
      combinedMarkersTable,
      sets = names(combinedMarkersTable)[6:2],
      nsets = 5,
      number.angles = 30,
      point.size = 2.5,
      line.size = 1,
      mainbar.y.label = "Marker Intersections",
      sets.x.label = "Conserved markers/subset",
      order.by = "freq",
      keep.order = TRUE,
      text.scale = c(1.5, 1.5, 1.2, 1.2, 1.75, 1.3)
    )
    print("line 150 findCMG")
    #remove unwanted objects

    rm(conservedMarkers_Tables_List)

    output$UpsetMarker <- renderPlot({
      upsetPlotMG
    })
    output$upsetLegend <- renderText({
      " UpSet plot showing the number of conserved marker genes that are shared
                   across subsets."
    })
    #show tables
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
    #Allow user to download tables
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
    # update progress bar value
    update_modal_progress(1)
  }

#PARSING TABLE FOR USE IN UPSETR**************************************
#this code iterates through each cell marker and record its precesnce/absence in each of
#the other subsets.
#for each table, add the 8th coloumn (STARTING FROM smallest TABLE)


tableParserCMG <- function(conservedMarkers_Tables_List, combinedMarkerGenesTableHeader) {
  #tables will be used to plot upsetplot
  combinedMarkersTable = c()

  #row that temporarily stores info about each marker gene from each subset
  #if a marker is present, the subset will have a score of 1. else, it will score a 0
  tmpMarkerRow = c()

  for (i in 1:length(conservedMarkers_Tables_List)) {
    #1) make a gene name column in each table
    conservedMarkers_Tables_List[[i]]$gene = rownames(conservedMarkers_Tables_List[[i]])
    #make a cluster_marker column in each table
    conservedMarkers_Tables_List[[i]]$Cluster_Marker = paste(conservedMarkers_Tables_List[[i]]$cluster,
                                                             conservedMarkers_Tables_List[[i]]$gene)
    #note that these changes are only visible from within the list
    print("line 234 cmg")
  }

  #for each table, starting at smallest
  for (i in 1:length(conservedMarkers_Tables_List)) {
    #2) iterate through each row in that table
    for (row in 1:nrow(conservedMarkers_Tables_List[[i]])) {
      #check if marker is already present in the combined table
      #col 15 holds cluster_marker info
      clusterMarker = paste(conservedMarkers_Tables_List[[i]][row, 15])
      coloumnWithPattern = which(grepl(clusterMarker, combinedMarkersTable))
      print("line 245 cmg")
      #if marker isn't present in combinedMarkersTable
      if (length(coloumnWithPattern) == 0) {
        #1)add it to tmpMarkerRow table
        tmpMarkerRow = rbind(tmpMarkerRow, clusterMarker)
        print("line 250 cmg")
        for (j in 1:length(conservedMarkers_Tables_List)) {
          #search if that marker exist in myTables
          x = which(grepl(clusterMarker, conservedMarkers_Tables_List[[j]]))
          print("line 254 cmg")
          #if it isn't present,
          if (length(x) == 0) {
            #add 0 to tmpMarkerRow
            tmpMarkerRow = cbind(tmpMarkerRow, 0)
          }
          else
            #add1 to tmp row
            tmpMarkerRow = cbind(tmpMarkerRow, 1)
          print("line 263 cmg")
        } #then search marker in next table
        #once done crosscomparing for selected marker,
        #add tmpMarkerRow to the combinedMarkersTable
        combinedMarkersTable = rbind(combinedMarkersTable, tmpMarkerRow)
        #reinitiate tmpMarkerRow to 0
        tmpMarkerRow = c()
        print("line 270 cmg")
      } #if it is present already it means you've already cross-compared it
      #so move on to next marker
    }
    #once done will ll rows of a table, move on to next table

  }
  print("line 277 cmg")
  combinedMarkersTable = data.frame(combinedMarkersTable)
  #finally, rename the coloumns in combinedMarkersTable
  names(combinedMarkersTable) = combinedMarkerGenesTableHeader

  #change all subsets from factor to numeric
  combinedMarkersTable[, 2:(length(conservedMarkers_Tables_List) + 1)] <-
    sapply(combinedMarkersTable[, 2:(length(conservedMarkers_Tables_List) +
                                       1)], as.character)
  combinedMarkersTable[, 2:(length(conservedMarkers_Tables_List) + 1)] <-
    sapply(combinedMarkersTable[, 2:(length(conservedMarkers_Tables_List) +
                                       1)], as.numeric)
  return(combinedMarkersTable)
}

