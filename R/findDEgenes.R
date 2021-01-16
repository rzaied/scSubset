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
              ),)
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


findDEgenes <- function(input, output, session, seuratObjectsList, dataset, dataset2, selectedNumGenes) {
  #UPSET PLOT FOR DE GENES

  #carrying on from Autoclustering_V2.R script,
  combinedDEgenesTableHeader = c("DE_gene")

  #list to store DE genes from each cell type
  DE_Tables_List = list()
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
    #update_modal_progress((i + 5) / 20)

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
    levels(seuratObjectsList[[i]])
    print("Line 95, find DE genes")
    #for each cluster in list
    for (j in 1:length(levelsList)) {

      #in try catch incase 1: no genes meet filtering criteria
        #2: some clusters exist in one dataset but not the other
      tryCatch( {

        DE.response <-
              FindMarkers(
                seuratObjectsList[[i]],
                ident.1 = paste(levelsList[j], "_", dataset, sep = ""),
                ident.2 = paste(levelsList[j], "_", dataset2, sep = ""),
                min.pct = 0.5,
                logfc.threshold = 0.25,
                verbose = TRUE
              )

      print("Line 112, find DE genes")
        #find the top 5 positive DE genes
        #or use all those that have pvalue 0.01, table with 0 if none.
        DE.response = DE.response[DE.response$p_val_adj <= 0.01,]
        print("Line 116, find DE genes")
        #get the top 5 most sig genes (or selectedNumGenes)
        DE.response = DE.response %>% top_n(n = -selectedNumGenes, wt = p_val_adj)
        #bottom_num=DE.response %>% top_n(n = -5, wt = avg_logFC)
        #DE_top holds DE genes per CLUSTER per SUBSET
        #add a coloumn with the identity of the cell
        DE.response$cluster <- levelsList[[j]]
        print("Line 127, find DE genes")
        #append those to a combinedDE_top table
        combinedDE_top = rbind(combinedDE_top, DE.response)
        #end of inner loop
        print("Line 130, find DE genes")

        },
      error=function(cond) {
          message(cond)
          # Choose a return value in case of error
          return(NA)
        },
        warning=function(cond) {
          message(cond)
          # Choose a return value in case of warning
          return(NULL)
        },
        finally={
        #things to execute regardless
      })
    }
    print("line 134 DEG")
    combinedDE_top = data.frame(combinedDE_top)
    print("Line 135, find DE genes")
    #append tables to list
    DE_Tables_List[[i]] = combinedDE_top
    #to create variables dynamically and assign them with respective
    #top 10 table "subset_2000top10
    # this is just so that the tables dont get overwritten and are saved in env
    #assign(paste("top_DE_",subsetSize, sep = ""), combinedDE_top)
    print("Line 142, find DE genes")
    #add subsetsize col
    combinedDE_top$subsetSize <- subsetSize
    #add gene name as seperate coloumn
    combinedDE_top$gene = rownames(combinedDE_top)
    sumStatDEtable1 = rbind(sumStatDEtable1, combinedDE_top)
    #once we compiled the top DE genes from each cluster of a subset, move on to next subset
    combinedDE_top = c()
    combinedDEgenesTableHeader = c(combinedDEgenesTableHeader, subsetSize)
    print("Line 152, find DE genes")
  }
  print("Line 153, find DE genes")
  rownames(sumStatDEtable1) <- c(1:nrow(sumStatDEtable1))
  sumStatDEtable1 = data.frame(sumStatDEtable1)
  #make gene names first coloumn
  sumStatDEtable1 <- sumStatDEtable1[, c(7, 1, 2, 3, 4, 5, 6, 8)]
  print("Line 159, find DE genes")


  combinedDEgenesTable <- tableParserDEG(DE_Tables_List, combinedDEgenesTableHeader)
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




tableParserDEG <- function(DE_Tables_List, combinedDEgenesTableHeader) {
  #table will be used to plot upsetplot
  combinedDEgenesTable = c()
  #row that temporarily stores info about each marker gene from each subset
  #if a marker is present, the subset will have a score of 1. else, it will score a 0

  tmpDErow = c()
  #put tables in a list to iterate through them


  #PARSING TABLE FOR USE IN UPSETR**************************************
  #this code iterates through each cell marker and record its precesnce/absence in each of
  #the other subsets.
  #for each table, add the 8th coloumn (STARTING FROM smallest TABLE)

  #******************DE GENES***************************************************************
  for (i in 1:length(DE_Tables_List)) {
    #1) make a gene name column in each table
    DE_Tables_List[[i]]$gene = rownames(DE_Tables_List[[i]])
    #make a cluster_marker column in each table
    DE_Tables_List[[i]]$Cluster_DEgene = paste(DE_Tables_List[[i]]$cluster,
                                               DE_Tables_List[[i]]$gene)
    #note that these changes are only visible from within the list
    print("Line 269, find DE genes")
  }

  #for each table, starting at smallest
  for (i in 1:length(DE_Tables_List)) {
    #2) iterate through each row in that table
    for (row in 1:nrow(DE_Tables_List[[i]])) {
      #check if marker is already present in the combined table
      #col 8 has the cluster_gene name info
      print("Line 278, find DE genes")
      clusterMarker = paste(DE_Tables_List[[i]][row, 8])
      coloumnWithPattern = which(grepl(clusterMarker, combinedDEgenesTable))

      #if marker isn't present in combinedDEgenesTable
      if (length(coloumnWithPattern) == 0) {
        #1)add it to tmpDErow table
        tmpDErow = rbind(tmpDErow, clusterMarker)
        print("Line 287, find DE genes")
        for (j in 1:length(DE_Tables_List)) {
          #search if that marker exists in topDE table of each subset
          x = which(grepl(clusterMarker, DE_Tables_List[[j]]))
          print("Line 290, find DE genes")
          #if it isn't present,
          if (length(x) == 0) {
            #add 0 to tmpDErow
            tmpDErow = cbind(tmpDErow, 0)
          }
          else
            #add1 to tmp row
            tmpDErow = cbind(tmpDErow, 1)

        } #then search marker in next table
        #once done crosscomparing for selected marker,
        #add tmpDErow to the combinedDEgenesTable
        combinedDEgenesTable = rbind(combinedDEgenesTable, tmpDErow)
        #reinitiate tmpDErow to 0
        tmpDErow = c()
      } #if it is present already it means you've already cross-compared it
      #so move on to next marker
    }
    #once done will ll rows of a table, move on to next table

  }

  print("Line 313, find DE genes")
  combinedDEgenesTable = data.frame(combinedDEgenesTable)
  #finally, rename the coloumns in combinedDEgenesTable
  names(combinedDEgenesTable) = combinedDEgenesTableHeader
  print("Line 317, find DE genes")
  #change all subsets from factor to numeric
  combinedDEgenesTable[, 2:(length(DE_Tables_List) + 1)] <-
    sapply(combinedDEgenesTable[, 2:(length(DE_Tables_List) + 1)], as.character)
  combinedDEgenesTable[, 2:(length(DE_Tables_List) + 1)] <-
    sapply(combinedDEgenesTable[, 2:(length(DE_Tables_List) + 1)], as.numeric)

  print("Line 324, find DE genes")
  return(combinedDEgenesTable)

}







#to know # of cells per condition+cluster
#table(seuratObjectsList[[i]]@meta.data$celltype.cond)
#save plot
#pdf("./Plots/RA_integrated10k_UpsetDE_rep10.pdf", height=4, width=6, paper = "USr")
#print(upsetPlotDE)
#dev.off()
