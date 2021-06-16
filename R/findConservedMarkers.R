#' Single Cell conserved markers Tab UI
#'
#' @export
#' @return None

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


#' Single Cell conserved markers Server
#'
#' @param seuratObjectsList Reactive value containing list of downsampled seurat objects
#' @param selectedNumGenes Reactive value containing user input for number of top genes to test for per cluster per subset
#' @export
#' @return Returns a reactive value containing a table that lists conserved markers per cluster per subset
#'

findConservedMarkers <-
  function (input,
            output,
            session,
            seuratObjectsList,
            selectedNumGenes) {

    #combinedtop conserved markers table holds marker genes for one subset at a time
    top_conservedMarkers_combined = c()
    #table holds all marker genes of all subsets
    sumStatMarkerTable1 = c()
    #for each seurat object, add data cloumns
    for (i in 1:(length(seuratObjectsList))) {
      DefaultAssay(seuratObjectsList[[i]]) <- "RNA"
      subsetSize = paste(nrow(seuratObjectsList[[i]]@meta.data) / 1000, "K", sep =
                           "")

      #update progress bar, following on from 0.5
      update_modal_progress((i + 5) / 10.5)
      #list to easily access the clusters in each subset
      levelsList = levels(seuratObjectsList[[i]])
      #for each cluster in list
      for (j in 1:length(levelsList)) {

        tryCatch({
          conserved.markers <-
            FindConservedMarkers(
              seuratObjectsList[[i]],
              ident.1 = levelsList[[j]],
              grouping.var = "orig.ident",
              verbose = TRUE,
              min.pct = 0.3,
              logfc.threshold = 0.3
            )

          #in case a cluster doesn't exist in a condition, only one condition will be used (i.e. only 5 cols)
          #if tested cluster is present in both conditions (i.e. there are 12 cols) carry on, otherwise skip that cluster
          if (ncol(conserved.markers) == 12)  {
            #only keep markers that are positively expressed in both conditions and use minimump_p_value, which is a combinded p_val to only keep significant markers
            conserved.markers <- conserved.markers %>%
              mutate("cluster" = levelsList[[j]]) %>%
              mutate("gene" = rownames(conserved.markers))

            conserved.markers<- conserved.markers[conserved.markers$dataset2_avg_log2FC >= 0.3 &
                                                    conserved.markers$dataset1_avg_log2FC >= 0.3 &
                                                    conserved.markers$minimump_p_val <= 0.01, ]  #select top5 (using many genes significantly increases computation time) #add a coloumn with the identity of the cluster

            top_conservedMarkers<- conserved.markers



            if (i == length(seuratObjectsList)) {
              #to only look at top 10 of
              top_conservedMarkers <- conserved.markers %>%
                slice_head(n = selectedNumGenes)

            }


            #append those to a top_conservedMarkers_combined table
            top_conservedMarkers_combined = rbind(top_conservedMarkers_combined,
                                                  top_conservedMarkers)
          }

        },
        error = function(cond) {

          # Choose a return value in case of error
          return(NA)
        },
        warning = function(cond) {

          # Choose a return value in case of warning
          return(NULL)
        },
        finally = {

        })
      }

      #add subsetsize col
      top_conservedMarkers_combined = top_conservedMarkers_combined %>%
        as.data.frame() %>%
        mutate("subsetSize" = subsetSize)

      sumStatMarkerTable1 = rbind(sumStatMarkerTable1, top_conservedMarkers_combined)

      #once we compiled the top marker genes from each cluster of a subset, move on to next subset
      top_conservedMarkers_combined = c()

    }

    #add numbered row names to table
    rownames(sumStatMarkerTable1) <- c(1:nrow(sumStatMarkerTable1))

    combinedMarkersTable <- tableParserCMG(sumStatMarkerTable1)
    #plot upsetPlot via combined table
    upsetPlotMG = upset(
      combinedMarkersTable,
      sets = names(combinedMarkersTable)[ncol(combinedMarkersTable):2],
      nsets = ncol(combinedMarkersTable)-1,
      number.angles = 30,
      point.size = 2.5,
      line.size = 1,
      mainbar.y.label = "Marker Intersections",
      sets.x.label = "Conserved markers/subset",
      order.by = "freq",
      keep.order = TRUE,
      text.scale = c(1.5, 1.5, 1.2, 1.2, 1.75, 1.3)
    )

    #remove unwanted objects
    #rm(conservedMarkers_Tables_List)

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
   # update_modal_progress(1)
    return(combinedMarkersTable)
  }

#PARSING TABLE FOR USE IN UPSETR**************************************

#' Function that iterates through cell markers and record their presence/absence per subset.
#'
#' @param sumStatMarkerTable1 Reactive value containing table that lists the conserved markers present in each subset
#'
#' @export
#' @return Returns a Reactive value listing conserved marker genes resolved per subset

tableParserCMG <- function(sumStatMarkerTable1) {
  #PARSING TABLE FOR USE IN UPSETR**************************************
  #this code record its presence/absence of each marker in each of the other subsets
  combinedMarkersTable <- sumStatMarkerTable1

  combinedMarkersTable$Cluster_Marker = paste(combinedMarkersTable$cluster,
                                              combinedMarkersTable$gene)

  combinedMarkersTable <- combinedMarkersTable %>% ungroup() %>%
    select(Cluster_Marker, subsetSize) %>% mutate(presence = 1) %>%
    pivot_wider(names_from = subsetSize , values_from = presence) %>%
    replace(is.na(.), 0) %>% as.data.frame()

  #to only look at
  combinedMarkersTable <- combinedMarkersTable %>%
    filter(!(combinedMarkersTable[, ncol(combinedMarkersTable)] == "0"))


  return(combinedMarkersTable)
}
