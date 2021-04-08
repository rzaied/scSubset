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
    #carrying on from findDEgenes_integration.R script,
    #combinedtop conserved markers table holds marker genes for one subset at a time
    top_conservedMarkers_combined = c()
    #table holds all marker genes of all subsets
    sumStatMarkerTable1 = c()
    #for each seurat object, add data cloumns
    for (i in 1:(length(seuratObjectsList))) {
      DefaultAssay(seuratObjectsList[[i]]) <- "RNA"
      subsetSize = paste(nrow(seuratObjectsList[[i]]@meta.data) / 1000, "K", sep =
                           "")
      print(subsetSize)
      #update progress bar, following on from 0.5
      #update_modal_progress((i + 5) / 10.5)
      #list to easily acess the clusters in each subset
      levelsList = levels(seuratObjectsList[[i]])
      #for each cluster in list
      for (j in 1:length(levelsList)) {
        #logfc.threshold is by default 0.25, increasing to 0.4 speeds up the function but could miss weaker signals
        #Find conserved marker genes (marker genes shared b/w both conditions)
        tryCatch({
          conserved.markers <-
            FindConservedMarkers(
              seuratObjectsList[[i]],
              ident.1 = levelsList[[j]],
              grouping.var = "orig.ident",
              verbose = TRUE,
              min.pct = 0.7,
              logfc.threshold = 0.3
            )
          print("line 81")

          #in case a cluster doesn't exist in a condition, only one condition will be used (i.e. only 5 cols)

          #if tested cluster is present in both conditions (i.e. there are 12 cols) carry on, otherwise skip that cluster
          if (ncol(conserved.markers) == 12)  {
            #only keep markers that are positively expressed in both conditions and use minimump_p_value, which is a combinded p_val to only keep significant markers
            conserved.markers = conserved.markers[conserved.markers$dataset2_avg_log2FC >= 0 &
                                                    conserved.markers$dataset1_avg_log2FC >= 0 &
                                                    conserved.markers$minimump_p_val <= 0.01, ] %>%  #select top5 (using many genes significantly increases computation time)
              mutate(cluster = levelsList[[j]]) #add a coloumn with the identity of the cluster

            top_conservedMarkers = sort_by_minimump_p_val(selectedNumGenes, conserved.markers)
            #append those to a top_conservedMarkers_combined table
            top_conservedMarkers_combined = rbind(top_conservedMarkers_combined,
                                                  top_conservedMarkers)
            print("line 99")
          }

        },
        error = function(cond) {
          # message(cond)
          # Choose a return value in case of error
          return(NA)
        },
        warning = function(cond) {
          #message(cond)
          # Choose a return value in case of warning
          return(NULL)
        },
        finally = {
          #things to execute regardless
        })
      }

      #append tables to list
      #  conservedMarkers_Tables_List[[i]] = combinedMarkers_top

      # this is just so that the tables dont get overwritten and are saved in env for me
      assign(paste("top_marker_", subsetSize, sep = ""),
             top_conservedMarkers_combined)
      #add subsetsize col
      top_conservedMarkers_combined = top_conservedMarkers_combined %>%
        as.data.frame() %>%
        mutate("subsetSize" = subsetSize) %>%
        mutate("gene" = rownames(top_conservedMarkers_combined))

      sumStatMarkerTable1 = rbind(sumStatMarkerTable1, top_conservedMarkers_combined)

      #once we compiled the top marker genes from each cluster of a subset, move on to next subset
      top_conservedMarkers_combined = c()

    }

    #add numbered row names to table
    rownames(sumStatMarkerTable1) <- c(1:nrow(sumStatMarkerTable1))

    #Round Pval

    combinedMarkersTable <- tableParserCMG(sumStatMarkerTable1)
    #finally, rename the coloumns in combinedMarkersTable
    #  names(combinedMarkersTable) = combinedMarkerGenesTableHeader
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
  #this code record its precesnce/absence of each marker in each of the other subsets
  combinedMarkersTable <- sumStatMarkerTable1

  combinedMarkersTable$Cluster_Marker = paste(combinedMarkersTable$cluster,
                                              combinedMarkersTable$gene)

  combinedMarkersTable <- combinedMarkersTable %>% ungroup() %>%
    select(Cluster_Marker, subsetSize) %>% mutate(presence = 1) %>%
    pivot_wider(names_from = subsetSize , values_from = presence) %>%
    replace(is.na(.), 0) %>% as.data.frame()

  print("line 256 findMG")
  return(combinedMarkersTable)
}


#' Function that sorts and select top n genes
#'
#' @param conserved.markers Reactive value containing conserved marker genes of a cluster in a given subset
#' @param selectedNumGenes Reactive value containing number of Deferentially expressed genes to be tested for
#' as chosen by the user
#' @export
#' @return Returns a Reactive value of the top conserved markers after sorting


sort_by_minimump_p_val <-
  function(selectedNumGenes, conserved.markers) {
    conserved.markers <-
      conserved.markers %>% #seperate coefficient and exponent to sort
      mutate("coefficient" = as.numeric(as.character(
        str_extract(conserved.markers$minimump_p_val, regex("([^e]+)"))
      ))) %>%
      mutate("exponent" = as.numeric(as.character(
        str_extract(conserved.markers$minimump_p_val, regex("[^e]*$"))
      ))) %>%
      arrange(exponent, coefficient) %>% select(-exponent,-coefficient) %>%
      slice_head(n = selectedNumGenes) #select top 5 rows (most statistically significant)


    return(conserved.markers)

  }
