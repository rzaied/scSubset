
#'  #' Single Cell Integration Tab Server
#'
#' @param seurat Reactive value containing seurat object
#' @param seurat2 Reactive value containing seurat object
#' @param mito Reactive value containing user input for mitochondrial genes pattern
#' @param res Reactive value containing user input for clustering resolution
#' @param dataset Reactive value containg user input of uploaded dataset name
#'
#' @export
#' @return Returns a Reactive value containing list of downsampled seurat objects
#'  with reduced dimensions (PCA data) and scaled counts
#'
scIntegrate <- function(input, output, session, seurat, seurat2, mito, res) {
  seuratObjectsList <- reactiveValues()

  print(mito)
  print(res)
  dim = 15

  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  #use ^MT- for human data and mt for mouse
  seurat[["percent.mt"]] = PercentageFeatureSet(seurat, pattern = mito)
  seurat2[["percent.mt"]] = PercentageFeatureSet(seurat2, pattern = mito)

  # Visualize QC metrics as a violin plot

  seurat = subset(seurat, subset = nFeature_RNA > 200 &
                    nFeature_RNA < 8000 & percent.mt < 14)
  seurat2 = subset(seurat2, subset = nFeature_RNA > 200 &
                     nFeature_RNA < 8000 & percent.mt < 14)

  print("line 44 integration")
 # seurat = SCTransform(seurat, vars.to.regress = "percent.mt", verbose = TRUE)

  seurat <- NormalizeData(seurat, verbose = TRUE)
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
#  seurat2 = SCTransform(seurat2, vars.to.regress = "percent.mt", verbose = TRUE)

  seurat2 <- NormalizeData(seurat2, verbose = TRUE)
  seurat2 <- FindVariableFeatures(seurat2, selection.method = "vst", nfeatures = 2000)

  seurat.anchors <-
    FindIntegrationAnchors(object.list = list(seurat, seurat2),
                           dims = 1:dim)
  seurat.combined <-
    IntegrateData(anchorset = seurat.anchors, dims = 1:dim)

  #remove not needed obj to save mem
  rm(seurat.anchors)

  DefaultAssay(seurat.combined) <- "integrated"

  # Run the standard workflow for visualization and clustering
  seurat.combined <- ScaleData(seurat.combined, verbose = TRUE)

  #here subset

  #************************************SUBSETTING OF CELLS***************************

  #get number of cells in dataset
  numCells = nrow(seurat.combined@meta.data)
  #least number of cells to test
  minSubset = round(numCells * 0.2 - 1)
  incrementation = minSubset

  #empty list to store generated subsets
  seuratObjectsList = c()

  x = 0
  #seq(starting nuber, ending number, incrementation)
  for (i in seq(from = minSubset, to = numCells, by = incrementation)) {
    x = x + 1
    #update progress bar
    #update_modal_progress(x / 25)
    print(i)

    #subsetting
    subset = subset(seurat.combined, cells = sample(Cells(seurat.combined), i))
    subset <- RunPCA(subset, npcs = 30, verbose = FALSE)
    #t-SNE and Clustering
    #subset<- RunTSNE(subset, reduction="pca")
    #TSNEPlot(subset, reduction= "tsne", label=TRUE)
    #TSNEPlot(subset, reduction= "tsne", split.by="orig.ident")
    subset <- RunUMAP(subset, reduction = "pca", dims = 1:dim)
    subset <- FindNeighbors(subset, reduction = "pca", dims = 1:dim)
    subset <- FindClusters(subset, resolution = res)

    # Visualization
    p1 <-
      DimPlot(subset, reduction = "umap", group.by = "orig.ident")
    p2 <- DimPlot(subset, reduction = "umap", label = TRUE)

    #DimPlot(subset, reduction = "umap", split.by = "orig.ident")
    assign(paste("UMAP", x, sep = ""), p2)

    #append subset to list
    seuratObjectsList = c(seuratObjectsList, subset)
    #Add subset to list

  }

  #call renaming clusters function
  seuratObjectsList <- renameClusters(seuratObjectsList)

  print("line 100 integr")
  combinedBarplot <- projectClusters(seuratObjectsList)

  print("line 102 integration")
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


  return(seuratObjectsList)
}

