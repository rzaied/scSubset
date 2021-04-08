#' Single Cell differential expression Tab UI
#'
#' @export
#' @return None
findDEgenesUI <- function(id) {
  #namespaces make it so that ids only need to be unique within a namespace and not across the app
  ns <- NS(id)
  tagList(fluidRow(column(
    12,
    wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 50px;",
              p(
                tags$b(
                  'Differential expression is only computed for integrated datasets',
                  style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;"
                )
              ), )
  )),
  fluidRow(column(
    8,
    wellPanel(
      style = "background-color: #fff; border-color: #2c3e50; height: 500px;",
      p(
        tags$b('UpSet plot of intersecting DE genes', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
      ),
      hr(),
      plotOutput(ns("upsetDE"), height = "360px") %>% withSpinner(color =
                                                                    "#0dc5c1"),
      textOutput(ns("upsetDELegend"))
    )
  )),
  fluidRow(column(
    12,
    wellPanel(
      style = "background-color: #fff; border-color: #2c3e50; height: 650px;",
      p(
        tags$b('DE genes present per subset ', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
      ),
      hr(),
      downloadLink(ns('downloadDataDE'), 'Download table as .csv'),
      DT::dataTableOutput(ns("DEtable")) %>% withSpinner(color =
                                                           "#0dc5c1")
    )
  )),
  fluidRow(column(
    12,
    wellPanel(
      style = "background-color: #fff; border-color: #2c3e50; height: 650px;",
      p(
        tags$b('DE genes and summary statistics per subset ', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
      ),
      downloadLink(ns('downloadDataDE2'), 'Download table as .csv'),
      hr(),
      DT::dataTableOutput(ns("sumStatDEtable")) %>% withSpinner(color =
                                                                  "#0dc5c1")
    )
  )))

}

#' Single Cell differential expression Server
#'
#' @param seuratObjectsList Reactive value containing list of downsampled seurat objects
#' @param dataset1_name Reactive value holding project name for first dataset
#' @param dataset1_name Reactive value holding project name for second dataset
#' @param selectedNumGenes Reactive value containing user input for number of top genes to test for per cluster per subset
#' @export
#' @return Returns a Reactive value containing a table that lists deferentially expressed genes present in each subset
#'

findDEgenes <-
  function(input,
           output,
           session,
           seuratObjectsList,
           dataset1_name,
           dataset2_name,
           selectedNumGenes) {
    #UPSET PLOT FOR DE GENES
    print("DE 59")
    #carrying on from Autoclustering_V2.R script,
  #  combinedDEgenesTableHeader = c("DE_gene")

    #list to store DE genes from each cell type
  #  DE_Tables_List = list()
    #names of clusters before reassigning which condition they belong to for findCMG
    #origlevelsList = list()
    #combinedtop DE genes table. this table will hold all DE genes per subset
    combinedDE_top = c()
    #this table will hold collective DE genes (of all subsets)
    sumStatDEtable1 = c()

    #for each seurat object, add data cloumns
    for (i in 1:(length(seuratObjectsList))) {
      DefaultAssay(seuratObjectsList[[i]]) <- "RNA"
      subsetSize = paste(nrow(seuratObjectsList[[i]]@meta.data) / 1000, "K", sep =
                           "")
      print(subsetSize)

      #update progress bar, following on from 0.2.
      update_modal_progress((i + 5) / 20)

      levelsList = levels(seuratObjectsList[[i]])
      #origlevelsList[[i]] = levelsList
      seuratObjectsList[[i]]$celltype.cond <-
        paste(Idents(seuratObjectsList[[i]]),
              seuratObjectsList[[i]]$orig.ident,
              sep = "_")
      seuratObjectsList[[i]]$celltype <-
        Idents(seuratObjectsList[[i]])
      #to know which sample each cluster comes from
      Idents(seuratObjectsList[[i]]) <- "celltype.cond"
      print(levels(seuratObjectsList[[i]]))
      print("Line 95, find DE genes")
      #for each cluster in list
      for (j in 1:length(levelsList)) {
        print(levelsList)

        #in try catch incase 1: no genes meet filtering criteria
        #2: some clusters exist in one dataset but not the other
        tryCatch({
          DE.response <-
            FindMarkers(
              seuratObjectsList[[i]],
              ident.1 = paste(levelsList[j], "_", dataset1_name, sep = ""),
              ident.2 = paste(levelsList[j], "_", dataset2_name, sep = ""),
              min.pct = 0.5,
              logfc.threshold = 0.25,
              verbose = TRUE
            )

          print("Line 112, find DE genes")
          #find the top 5 positive DE genes
          #or use all those that have pvalue 0.01, table with 0 if none
          #get the top 5 most sig genes (or -selectedNumGenes) via pvalue so lowest pvalues
          DE.response = DE.response[DE.response$p_val_adj <= 0.01, ] %>%  #add a coloumn with the identity of the cell
            mutate(cluster=  levelsList[[j]])

          combinedDE_top=sort_by_p_val_adj(selectedNumGenes, DE.response)


          print("Line 116, find DE genes")

          #DE_top holds DE genes per CLUSTER per SUBSET
          #append those to a combinedDE_top table
          combinedDE_top = rbind(combinedDE_top, DE.response)
          #end of inner loop
          print("Line 130, find DE genes")

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
      print("line 134 DEG")
      combinedDE_top = data.frame(combinedDE_top)
      print("Line 135, find DE genes")
      #append tables to list
  #    DE_Tables_List[[i]] = combinedDE_top
      print("Line 142, find DE genes")
      #add subsetsize col
      combinedDE_top$subsetSize <- subsetSize
      #add gene name as seperate coloumn
      combinedDE_top$gene = rownames(combinedDE_top)
      sumStatDEtable1 = rbind(sumStatDEtable1, combinedDE_top)
      #once we compiled the top DE genes from each cluster of a subset, move on to next subset
      combinedDE_top = c()
   #   combinedDEgenesTableHeader = c(combinedDEgenesTableHeader, subsetSize)
      print("Line 152, find DE genes")
    }
    print("Line 153, find DE genes")
    #add row names to table
    rownames(sumStatDEtable1) <- c(1:nrow(sumStatDEtable1))
    sumStatDEtable1 = data.frame(sumStatDEtable1)
    #make gene names first coloumn
    sumStatDEtable1 <- sumStatDEtable1[, c(7, 1, 2, 3, 4, 5, 6, 8)]
    print("Line 159, find DE genes")


    combinedDEgenesTable <- tableParserDEG(sumStatDEtable1)
    #finally, rename the coloumns in combinedDEgenesTable
    names(combinedDEgenesTable) = combinedDEgenesTableHeader

    print("LINE 165 findDEG")
    #plot upsetPlot via combined table
    upsetPlotDE = upset(
      combinedDEgenesTable,
      sets = names(combinedDEgenesTable)[6:2],
      nsets = 5,
      number.angles = 20,
      point.size = 2.5,
      line.size = 1,
      mainbar.y.label = "DE genes intersections",
      sets.x.label = "DE genes per subset",
      order.by = "freq",
      keep.order = TRUE,
      text.scale = c(1.5, 1.5, 1.2, 1.2, 1.75, 1.3)
    )

    #remove not needed obj to save mem
    rm(DE.response)
    rm(DE_Tables_List)

    print("Line 181, find DE genes")

    output$upsetDE <- renderPlot({
      upsetPlotDE
    })
    output$upsetDELegend <- renderText({
      " UpSet plot showing the number of differentially expressed genes that are shared across subsets."
    })
    output$sumStatDEtable <-
      renderDT({
        sumStatDEtable1
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

    output$DEtable <-
      renderDT({
        combinedDEgenesTable
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
    #ALLOW users to download tables
    output$downloadDataDE <- downloadHandler(
      filename = function() {
        paste('DEgenesTable_', Sys.Date(), '.csv', sep = '')
      },
      content = function(con) {
        write.csv(combinedDEgenesTable, con)
      }
    )

    output$downloadDataDE2 <- downloadHandler(
      filename = function() {
        paste('DEsumStat_', Sys.Date(), '.csv', sep = '')
      },
      content = function(con) {
        write.csv(sumStatDEtable1, con)
      }
    )



    return(combinedDEgenesTable)
  }



#' Function that iterates through differentially expressed genes and records their presence/absence per cluster per subset.
#'
#' @param sumStatDEtable1 Reactive value containing table that lists the differentially expressed genes present in each subset
#' @export
#' @return Returns a Reactive value listing conserved marker genes resolved per cluster per subset

tableParserDEG <- function(sumStatDEtable1) {
  #PARSING TABLE FOR USE IN UPSETR**************************************
  #this code record its precesnce/absence of each marker in each of the other subsets
  combinedDEgenesTable <- sumStatDEtable1

  combinedDEgenesTable$Cluster_Marker = paste(combinedDEgenesTable$cluster,
                                              combinedDEgenesTable$gene)

  combinedDEgenesTable <- combinedDEgenesTable %>% ungroup() %>%
    select(Cluster_Marker, subsetSize) %>% mutate(presence=1) %>%
    pivot_wider(names_from=subsetSize , values_from=presence) %>%
    replace(is.na(.), 0) %>% as.data.frame()

  print("line 256 findMG")

  return(combinedDEgenesTable)
}


#' Function that sorts and select top n genes
#'
#' @param DE.response Reactive value containing deferentially expressed genes of a cluster in a given subset
#' @param selectedNumGenes Reactive value containing number of Deferentially expressed genes to be tested for
#' as chosen by the user
#' @export
#' @return Returns a Reactive value of the top deferentially expressed genes after sorting


sort_by_p_val_adj<-function(selectedNumGenes, DE.response) {

  DE.response<- DE.response %>% #seperate coefficient and exponent to sort
    mutate("coefficient"=as.numeric(as.character(str_extract(DE.response$p_val_adj, regex("([^e]+)"))))) %>%
    mutate("exponent"=as.numeric(as.character(str_extract(DE.response$p_val_adj, regex("[^e]*$"))))) %>%
    arrange(exponent, coefficient) %>% select(-exponent, -coefficient) %>%
    slice_head(n= selectedNumGenes) #select top 5 rows (most statistically significant)


  return(DE.response)

}
