#' Single Cell compute cost Tab UI
#'
#' @export
#' @return None
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


#' Single Cell compute cost Tab Server
#' @param analysis_type Reactive value containing type of analysis (integrated vs single analysis)
#' @param combinedConsMarkersTable Reactive value containing a table that lists conserved markers per cluster per
#'  subset if the "find conserved markers" option was selected
#' @param seuratObjectsList Reactive value containing list of downsampled seurat objects
#'  with reduced dimensions (PCA data), scaled counts, and cluster labels that corresponds across subsets
#' @param costPerMil Reactive value containing cost of sequencing one million reads
#' @param depthPerCell Reactive value containing the depth of sequencing per cell selected by the user
#'
#' @export
#'

#server function
computeCost <-
  function(input, output, session, analysis_type, combinedConsMarkersTable,
           markersTable, seuratObjectsList, costPerMil, depthPerCell) {

    scCostTable=c()
    #for each subset
    for (i in 1:(length(seuratObjectsList))) {

      #find number of cells in dataset
      numCells = nrow(seuratObjectsList[[i]]@meta.data)
      #find number of clusters
      numClusters=max(as.numeric(levels
                                 (seuratObjectsList[[i]]@meta.data$seurat_clusters)[seuratObjectsList[[i]]@meta.data$seurat_clusters]))

      numClusters = numClusters+1 #since count starts at 0
      #find total number of marker genes in the reference dataset
      reference_num_markers=nrow(markersTable[!(markersTable[,ncol(markersTable)]==0),])


      #remove rows where marker doesn't exist (marker==0) in both reference and query datasets
      markersTable2<-markersTable[!(markersTable[,(i+1)]==0 & markersTable[,ncol(markersTable)]==0),]

      #calculate the overlapb between query and reference datasets
      markers_overlap= colSums(markersTable2[,(i+1), drop=FALSE]==markersTable2[,ncol(markersTable2), drop=FALSE])/reference_num_markers*100
      markers_overlap<-round(markers_overlap, digits = 2)

      #find total depth per subset
      depth=depthPerCell*numCells

      #find total cost per subset
      cost=((costPerMil*depth)/1000000)
      cost<-round(cost, digits = 3)

      #repeat for conserved marker genes, if user chose to compute it
      if (!combinedConsMarkersTable==0) {
        reference_num_consv_markers=nrow(combinedConsMarkersTable[!(combinedConsMarkersTable[,ncol(combinedConsMarkersTable)]==0),])
        combinedConsMarkersTable2<-combinedConsMarkersTable[!(combinedConsMarkersTable[,(i+1)]==0 & combinedConsMarkersTable[,ncol(combinedConsMarkersTable)]==0),]

        consv_markers_overlap= colSums(combinedConsMarkersTable2[,(i+1), drop=FALSE]==combinedConsMarkersTable2[,ncol(combinedConsMarkersTable), drop=FALSE])/reference_num_consv_markers*100
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

      names(scCostTable)=c("Number of cells", "number of clusters", " % Overlapping DE genes", "cost estimate")
    }

    if (analysis_type=="integration" && !combinedConsMarkersTable==0) {

      names(scCostTable)=c("Number of cells", "number of clusters", " % Overlapping conserved markers", " % Overlapping DE genes", "cost estimate")
    }
    if (analysis_type=="Single dataset") {
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

