#' Single Cell Cluster Tab UI
#'
#' @export
#' @return None


scClusteringUI <- function(id) {
  #namespaces make it so that ids only need to be unique within a namespace and not across the app
  ns <- NS(id)
  tagList(
    fluidRow(column(
      8,
      wellPanel(
        style = "background-color: #fff; border-color: #2c3e50; height: 515px;",
        p(
          tags$b('Comparing clustering solutions', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
        ),
        hr(),
        plotOutput(ns("barplot"), height = "360px") %>% withSpinner(color =
                                                                      "#0dc5c1"),
        textOutput(ns("barplotLegend"))
      )
    )),
    fluidRow(column(
      6,
      wellPanel(
        style = "background-color: #fff; border-color: #2c3e50; height: 500px;",
        p(
          tags$b('Projecting 20% of the dataset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
        ),
        hr(),
        plotOutput(ns("UMAP1"), height = "360px") %>% withSpinner(color =
                                                                    "#0dc5c1"),
      )
    ),
    column(6, fluidRow(
      column(
        12,
        wellPanel(
          style = "background-color: #fff; border-color: #2c3e50; height: 500px;",
          p(
            tags$b('Projecting 40% of the dataset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
          ),
          hr(),
          plotOutput(ns("UMAP2"), height = "360px") %>% withSpinner(color =
                                                                      "#0dc5c1")
        )
      )
    ))),
    fluidRow(column(
      6,
      wellPanel(
        style = "background-color: #fff; border-color: #2c3e50; height: 500px;",
        p(
          tags$b('Projecting 60% of the dataset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
        ),
        hr(),
        plotOutput(ns("UMAP3"), height = "360px") %>% withSpinner(color =
                                                                    "#0dc5c1")
      )
    ),
    column(6, fluidRow(
      column(
        12,
        wellPanel(
          style = "background-color: #fff; border-color: #2c3e50; height: 500px;",
          p(
            tags$b('Projecting 80% of the dataset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
          ),
          hr(),
          plotOutput(ns("UMAP4"), height = "360px") %>% withSpinner(color =
                                                                      "#0dc5c1")
        )
      )
    ))),
    fluidRow(column(
      6,
      wellPanel(
        style = "background-color: #fff; border-color: #2c3e50; height: 500px;",
        p(
          tags$b('Projecting 100% of the dataset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
        ),
        hr(),
        plotOutput(ns("UMAP5"), height = "360px") %>% withSpinner(color =
                                                                    "#0dc5c1")
      )
    ))
  )
}

#' Single Cell Variable Genes Tab UI
#'
#' @export
#' @return None

varGenesUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(column(
      12,
      wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 50px;",
                p(
                  tags$b('Variable genes are only computed for single datasets', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
                ), )
    )),
    fluidRow(column(
      6,
      wellPanel(
        style = "background-color: #fff; border-color: #2c3e50; height: 500px;",
        p(
          tags$b('Variable genes in 20% of the dataset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
        ),
        hr(),
        plotOutput(ns("variableGenes1"), height = "360px") %>% withSpinner(color =
                                                                             "#0dc5c1")
      )
    ),
    column(6, fluidRow(
      column(
        12,
        wellPanel(
          style = "background-color: #fff; border-color: #2c3e50; height: 500px;",
          p(
            tags$b('Variable genes in 40% of the dataset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
          ),
          hr(),
          plotOutput(ns("variableGenes2"), height = "360px") %>% withSpinner(color =
                                                                               "#0dc5c1")
        )
      )
    ))),
    fluidRow(column(
      6,
      wellPanel(
        style = "background-color: #fff; border-color: #2c3e50; height: 500px;",
        p(
          tags$b('Variable genes in 60% of the dataset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
        ),
        hr(),
        plotOutput(ns("variableGenes3"), height = "360px") %>% withSpinner(color =
                                                                             "#0dc5c1")
      )
    ),
    column(6, fluidRow(
      column(
        12,
        wellPanel(
          style = "background-color: #fff; border-color: #2c3e50; height: 500px;",
          p(
            tags$b('Variable genes in 80% of the dataset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
          ),
          hr(),
          plotOutput(ns("variableGenes4"), height = "360px") %>% withSpinner(color =
                                                                               "#0dc5c1")
        )
      )
    ))),
    fluidRow(column(
      6,
      wellPanel(
        style = "background-color: #fff; border-color: #2c3e50; height: 500px; display:center-align;",
        p(
          tags$b('Variable genes in 100% of the dataset', style = "font-size: 100%; font-family:Helvetica; color:#4c4c4c; text-align:left;")
        ),
        hr(),
        plotOutput(ns("variableGenes5"), height = "360px") %>% withSpinner(color =
                                                                             "#0dc5c1")
      )
    ))
  )
}


#' Single Cell Cluster Tab Server
#'
#' @param seurat Reactive value containing seurat object
#' @param mito Reactive value containing user input for mitochondrial genes pattern
#' @param res Reactive value containing user input for clustering resolution
#' @param dataset1_name Reactive value containing user input of uploaded dataset name
#'
#' @export
#' @return Returns a Reactive value containing list of downsampled seurat objects
#'  with reduced dimensions (PCA data) and scaled counts
#'
scClustering <-
  function (input,
            output,
            session,
            seurat,
            mito,
            res,
            dataset1_name) {
    seuratObjectsList <- reactiveValues()

    print(mito)
    print(res)

    print("autoclust line 211")

    #assign dataset name (for "project" in seurat)
    print(dataset1_name)
    seurat
    dim = 15
    print("line 193, autoclustering")
    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    #use ^MT- for human data and mt for mouse
    seurat[["percent.mt"]] = PercentageFeatureSet(seurat, pattern = mito)

    seurat = subset(seurat, subset = nFeature_RNA > 200 &
                      nFeature_RNA < 8000 & percent.mt < 14)

    print("line 201, autoclustering")
    #Normalization
    seurat = NormalizeData(seurat,
                           normalization.method = "LogNormalize",
                           scale.factor = 10000)
    seurat = NormalizeData(seurat)

    #scaling data
    all.genes = rownames(seurat)
    seurat = ScaleData(seurat, features = all.genes)

    #SUBSETTING OF CELLS

    #get number of cells in dataset
    numCells = nrow(seurat@meta.data)
    #least number of cells to test
    minSubset = round(numCells * 0.2 - 1)
    incrementation = minSubset
    #empty list to store generated subsets
    seuratObjectsList = c()
    print("line 221, autoclustering")
    x = 0
    #seq(starting nuber, ending number, incrementation)
    for (i in seq(from = minSubset, to = numCells, by = incrementation)) {
      print(i)
      x = x + 1
      #update progress bar
      update_modal_progress(x / 10)
      print(x)
      #subsetting
      subset = subset(seurat, cells = sample(Cells(seurat), i))

      #Plot variable features
      subset = FindVariableFeatures(subset,
                                    selection.method = "vst",
                                    nfeatures = 2000)

      #Run PCA on scaled data
      subset = RunPCA(subset, features = VariableFeatures(object = subset))
      #clustering
      subset = FindNeighbors(subset, dims = 1:dim)
      subset = FindClusters(subset, resolution = res)
      subset = RunUMAP(subset, dims = 1:dim)

      #Plot UMAP
      #p1 <- DimPlot(subset, reduction = "umap", group.by = "orig.ident")
      p2 <- DimPlot(subset, reduction = "umap", label = TRUE)

      #DimPlot(subset, reduction = "umap", split.by = "orig.ident")
      assign(paste("UMAP", x, sep = ""), p2)
      #DimPlot(subset, reduction = "umap", label=T, pt.size = 1) + NoLegend()

      #append subset to list
      seuratObjectsList = c(seuratObjectsList, subset)
      #Add subset to list
    }
    #call renaming clusters function
    seuratObjectsList <- renameClusters(seuratObjectsList)



    combinedBarplot <- projectClusters(seuratObjectsList)

    print("line 283 autoclust")

    output$barplot <- renderPlot({
      combinedBarplot
    })
    output$barplotLegend <- renderText({
      "Bar plot showing similarity score of subset and full dataset. Here, the cell labels
                   from each subset are compared against the cell labels of
                   the reference dataset (full dataset); values closest to 1 are most similar to the
                   reference. ARI: Adjusted Rand Index; NMI: Normalized Mutual Information."
    })
    output$UMAP1 <- renderPlot({
      UMAP1
    })
    output$UMAP2 <- renderPlot({
      UMAP2
    })
    output$UMAP3 <- renderPlot({
      UMAP3
    })
    output$UMAP4 <- renderPlot({
      UMAP4
    })
    output$UMAP5 <- renderPlot({
      UMAP5
    })

    #return renamed clusters
    return(seuratObjectsList)

  }



#' Single Cell Variable genes Tab Server
#'
#' @param seuratObjectsList Reactive value containg list of downsampled seurat objects
#'  with reduced dimensions (PCA data) and reduced dimensions
#'
#' @export
#'
scVarGenes <- function(input, output, session, seuratObjectsList) {
  #for each object
  for (i in 1:5) {
    print(i)

    top10 = head(VariableFeatures(seuratObjectsList[[i]]), 10)

    # plot variable genes
    plot1 <- VariableFeaturePlot(seuratObjectsList[[i]])
    plot2 <- LabelPoints(plot = plot1,
                         points = top10,
                         repel = TRUE)
    assign(paste("variableGenes", i, sep = ""), plot2)

  }

  output$variableGenes1 <- renderPlot({
    variableGenes1
  })
  output$variableGenes2 <- renderPlot({
    variableGenes2
  })
  output$variableGenes3 <- renderPlot({
    variableGenes3
  })
  output$variableGenes4 <- renderPlot({
    variableGenes4
  })
  output$variableGenes5 <- renderPlot({
    variableGenes5
  })


}

#' Function to find Mode
#'
#' @param x a vector of integer values
#'
#' @export

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#RENAMING CLUSTERS
#function to rename clusters after subsetting
#this program
#selects cells that belong to a cluster,x, from the reference dataset.
#finds these cells in the query dataset
#finds which cluster the majority of these cells belong to in the query dataset
#renames that cluster to x
#repeats process for other clusters in reference until all clusters in the query
#dataset are renamed.
#this insures clusters are not mislabeled and allows to project query and reference clusters
#avaiding mislabelling errors


#' Function to reassign cluster labels across downsampled subsets
#'
#' @param seuratObjectsList Reactive value containg list of downsampled seurat objects
#'  with reduced dimensions (PCA data) and scaled counts
#'
#' @export
#' @return Returns a Reactive value containing list of downsampled seurat objects
#'  with reassigned cluster labels that correspoond across downsampled subsets



renameClusters <- function(seuratObjectsList) {
  ref_subset=length(seuratObjectsList)
  #function to rename clusters after subsetting
  print("line 354 autoclustering")
  #for each object
  for (i in (length(seuratObjectsList)):1) {
    #make an empty list of size [max cluster#] to store new names of query subset
    newNames_list <-
      vector(mode = "character", length = max(as.numeric(seuratObjectsList[[i]]@active.ident)))
    #for each cluster in second object
    for (j in 0:max(as.numeric(as.character(seuratObjectsList[[i]]@active.ident)))) {
      #select all cell_ids where seurat_cluster = j from query subset
      clusterTable <-
        filter(
          seuratObjectsList[[i]]@meta.data,
          seuratObjectsList[[i]]@meta.data$seurat_clusters == j
        )
      print("line 369 autoclustering")
      #find corresponding cell ids in the reference dataset via merging
      overlappingCells = merge(clusterTable, seuratObjectsList[[ref_subset]]@meta.data, by = 0)

      #Which clusters do most of these cells belong to in reference?
      predominantCluster = as.numeric(as.character(Mode(
        overlappingCells$seurat_clusters.y
      )))
      print("line 378 autoclustering")
      #In position of predominantcluster, store cluster label, j (as per reference label)
      #list starts at 1 while clusters start at 0 so predominantCluster+1
      newNames_list[[j + 1]] = predominantCluster
      #repeat for next cluster in first table
      newNames_list
    }
    print("line 385 autoclustering")
    #Rename all clusters in first subset as per reference subset
    names(newNames_list) = levels(seuratObjectsList[[i]])
    seuratObjectsList[[i]] = RenameIdents(seuratObjectsList[[i]], newNames_list)
    levels(x = seuratObjectsList[[i]])
  }
  rm(clusterTable)
  return(seuratObjectsList)
}

#' Function to project downsampled subsets against the parent dataset and score projection quality
#'
#' @param seuratObjectsList Reactive value containg list of downsampled seurat objects
#'  with reduced dimensions (PCA data), scaled counts, and cluster labels that corresponds across subsets
#'
#' @export
#' @return Returns a bar plot showing projection quality scores across subsets
#'
#'
projectClusters <- function(seuratObjectsList) {
  #CALCULTE PROJECTION QUALITY************************************
  #table to store cluster similarity calculation
  ref_subset=length(seuratObjectsList)
  projectionQualityTable = c()
  for (i in 1:(length(seuratObjectsList))) {
    print("line 404 autoclustering")
    overlappingCells = merge(seuratObjectsList[[i]]@active.ident,
                             seuratObjectsList[[ref_subset]]@active.ident,
                             by =
                               0)
    print("line 409 autoclustering")


    overlappingCells$x<-as.numeric(as.character(overlappingCells$x))
    overlappingCells$y<-as.numeric(as.character(overlappingCells$y))

    #compute ARI and NMI
    ARI = ARI(overlappingCells$x, overlappingCells$y)
    NMI = NMI(overlappingCells$x, overlappingCells$y)
    print("439 autoclust")

    subsetSize = paste(nrow(seuratObjectsList[[i]]@meta.data) / 1000, "K", sep =
                         "")
    projectionQualityTable = rbind(projectionQualityTable, c(i, subsetSize, ARI, NMI))

    print("444 autoclust")
  }
  #num column is used to specify levels
  colnames(projectionQualityTable) = c("num","Subset", "ARI", "NMI")
  projectionQualityTable = data.frame(projectionQualityTable)
  projectionQualityTable$ARI = as.numeric(as.character(projectionQualityTable$ARI))
  projectionQualityTable$NMI = as.numeric(as.character(projectionQualityTable$NMI))

  print("454 autoclust")
  #to reorder the table by increasing size of subset
 projectionQualityTable$Subset <-
   factor(projectionQualityTable$Subset, levels = projectionQualityTable$Subset[order(projectionQualityTable$num)])
  #melting to plot multiple y values
  projectionQualityTable = melt(projectionQualityTable[,2:4], id.vars = 'Subset')
  colnames(projectionQualityTable)[2] = "Key"

  print("462 autoclust")
  #plot barplot
  combinedBarplot = ggplot(projectionQualityTable, aes(x = Subset, y = value, fill =
                                                         Key)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_minimal() + scale_fill_brewer(palette = "BuPu")

  print("472 autoclust")

  #**************************
  return(combinedBarplot)
}
