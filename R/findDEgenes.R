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

    #combinedtop DE genes table. this table will hold all DE genes per subset
    combinedDE_top = c()
    #this table will hold collective DE genes (of all subsets)
    sumStatDEtable1 = c()

    #for each seurat object, add data cloumns
    for (i in 1:(length(seuratObjectsList))) {
      DefaultAssay(seuratObjectsList[[i]]) <- "RNA"
      subsetSize = paste(nrow(seuratObjectsList[[i]]@meta.data) / 1000, "K", sep =
                           "")


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

      #for each cluster in list
      for (j in 1:length(levelsList)) {


        #in try catch in case 1: no genes meet filtering criteria
        #2: some clusters exist in one dataset but not the other
        tryCatch({
          DE.response <-
            FindMarkers(
              seuratObjectsList[[i]],
              ident.1 = paste(levelsList[j], "_", dataset1_name, sep = ""),
              ident.2 = paste(levelsList[j], "_", dataset2_name, sep = ""),
              min.pct = 0.3,
              logfc.threshold = 0.3,
        #      test.use = "MAST",
              verbose = TRUE
            )

          DE.response<- DE.response %>%
            mutate("cluster"=  levelsList[[j]]) %>%
            mutate("gene"= rownames(DE.response))

          DE.response = DE.response[DE.response$p_val_adj <= 0.01, ] #add a coloumn with the identity of the cell
          DE.response_top <-DE.response

          #DE_top holds DE genes per CLUSTER per SUBSET
          #append those to a combinedDE_top table
          combinedDE_top = rbind(combinedDE_top, DE.response_top)


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
          #things to execute regardless
        })
      }

      combinedDE_top = data.frame(combinedDE_top)
      #add subsetsize col
      combinedDE_top$subsetSize <- subsetSize
      sumStatDEtable1 = rbind(sumStatDEtable1, combinedDE_top)
      #once we compiled the top DE genes from each cluster of a subset, move on to next subset
      combinedDE_top = c()
    }

    #add row names to table
    rownames(sumStatDEtable1) <- c(1:nrow(sumStatDEtable1))
    sumStatDEtable1 = data.frame(sumStatDEtable1)
    #make gene names first coloumn
    sumStatDEtable1 <-sumStatDEtable1 %>%
      relocate(gene, .before= p_val)

    combinedDEgenesTable <- tableParserDEG(sumStatDEtable1)

    #plot upsetPlot via combined table
    upsetPlotDE = upset(
      combinedDEgenesTable,
      sets = names(combinedDEgenesTable)[ncol(combinedDEgenesTable):2],
      nsets = ncol(combinedDEgenesTable)-1,
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


  return(combinedDEgenesTable)
}


#' Function that sorts and select top n genes
#'
#' @param DE.response Reactive value containing deferentially expressed genes of a cluster in a given subset
#' @param selectedNumGenes Reactive value containing number of Deferentially expressed genes to be tested for
#' as chosen by the user
#' @export
#' @return Returns a Reactive value of the top deferentially expressed genes after sorting


# sort_by_p_val_adj<-function(selectedNumGenes, DE.response) {
#
#   DE.response<- DE.response %>% #seperate coefficient and exponent to sort
#     mutate("coefficient"=as.numeric(as.character(str_extract(DE.response$p_val_adj, regex("([^e]+)"))))) %>%
#     mutate("exponent"=as.numeric(as.character(str_extract(DE.response$p_val_adj, regex("[^e]*$"))))) %>%
#     arrange(exponent, coefficient) %>% select(-exponent, -coefficient) %>%
#     slice_head(n= selectedNumGenes) #select top 5 rows (most statistically significant)
#
#
#   return(DE.response)
#
# }
