computeCostUI <- function(id) {
  ns <- NS(id)
  tagList(
  fluidRow(column(
    12,
    wellPanel(
      style = "background-color: #fff; border-color: #2c3e50; height: 650px;",
      p(
        tags$b('Sequencing cost per subset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
      ),
      downloadLink(ns('downloadCostData'), 'Download table as .csv'),
      hr(),
      DT::dataTableOutput(ns("costTable")) %>% withSpinner(color =
                                                               "#0dc5c1")
    )
  )))

}

#server function
computeCost <-
  function(input, output, session, analysis_type,combinedConsMarkersTable,
           combinedMarkersTable, seuratObjectsList, costPerMil, depthPerCell) {

    scCostTable=c()
    #for each subset
    for (i in 1:(length(seuratObjectsList))) {
      print(i)
      #find number of cells in dataset
      numCells = nrow(seuratObjectsList[[i]]@meta.data)
      print(numCells)
      #find number of clusters
      numClusters=max(as.numeric(levels
                                 (seuratObjectsList[[i]]@meta.data$seurat_clusters)[seuratObjectsList[[i]]@meta.data$seurat_clusters]))


      #find total number of marker genes in the reference dataset
      reference_num_markers=nrow(combinedMarkersTable[!(combinedMarkersTable[,6]==0),])


      #remove rows where marker doesn't exist (marker==0) in both reference and query datasets
      combinedMarkersTable2<-combinedMarkersTable[!(combinedMarkersTable[,(i+1)]==0 & combinedMarkersTable[,6]==0),]

      #calculate the overlapb between query and reference datasets
      markers_overlap= colSums(combinedMarkersTable2[,(i+1), drop=FALSE]==combinedMarkersTable2[,6, drop=FALSE])/reference_num_markers*100
      markers_overlap<-round(markers_overlap, digits = 2)



      #find total depth per subset
      depth=depthPerCell*numCells

      #find total cost per subset
      cost=((costPerMil*depth)/1000000)
      cost<-round(cost, digits = 3)

      #repeat for conserved marker genes, if user chose to compute it
      if (!combinedConsMarkersTable==0) {
        reference_num_consv_markers=nrow(combinedConsMarkersTable[!(combinedConsMarkersTable[,6]==0),])
        combinedConsMarkersTable2<-combinedConsMarkersTable[!(combinedConsMarkersTable[,(i+1)]==0 & combinedConsMarkersTable[,6]==0),]

          consv_markers_overlap= colSums(combinedConsMarkersTable2[,(i+1), drop=FALSE]==combinedConsMarkersTable2[,6, drop=FALSE])/reference_num_consv_markers*100
          consv_markers_overlap<-round(consv_markers_overlap, digits = 2)

          scCostTable=rbind(scCostTable, c(numCells, numClusters,consv_markers_overlap, markers_overlap, cost))

      }
      #append to table
      if (combinedConsMarkersTable==0) {
        scCostTable=rbind(scCostTable, c(numCells, numClusters, markers_overlap, cost))
    }
}
    scCostTable=data.frame(scCostTable)

    if (analysis_type=="integration" && combinedConsMarkersTable==0) {
      print(analysis_type)
    names(scCostTable)=c("Number of cells", "number of clusters", " % Overlapping DE genes", "cost estimate")
    }

    if (analysis_type=="integration" && !combinedConsMarkersTable==0) {
      print(analysis_type)
      names(scCostTable)=c("Number of cells", "number of clusters", " % Overlapping conserved markers", " % Overlapping DE genes", "cost estimate")
    }
    if (analysis_type=="Single dataset") {
      print(analysis_type)
      names(scCostTable)=c("Number of cells", "number of clusters", " % Overlapping markers", "cost estimate")
    }


    #table output
    output$costTable <-
      renderDT({
        scCostTable
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



    output$downloadCostData <- downloadHandler(
      filename = function() {
        paste('costTable_', Sys.Date(), '.csv', sep = '')
      },
      content = function(con) {
        write.csv(scCostTable, con)
      }
    )

    }