#'  #' Single Cell .h5 data reading Tab Server
#'
#' @param path Reactive value containing path to matrix files
#' @param dataset_name Reactive value containing dataset name
#'
#' @export
#' @return Returns a Reactive value containing seurat object
#'
read_h5 <- function(input, output, session, dataset_name, path) {

    print("line 22")
    filePath <- paste0(path, "/", "0.h5")
    print(filePath)
    tryCatch({
      seurat.data <- Read10X_h5(filePath)
      seurat<- make_seurat_obj(seurat.data, dataset_name)
      return(seurat)
    },
    error = function(cond) {
      print("line 18")
      sendSweetAlert(
        session = session,
        title = "Data format error",
        text = "Ensure that 10X Genomics data were appropriately chosen:
              accepted files are either .h5 file or matrix.mtx, genes.tsv/features.tsv and barcodes.tsv files ",
        type = "error"
      )
      # return()
    })
  # })

}

#'  #' Single Cell 10x data reading Tab Server
#'
#' @param dataset_name Reactive value containing dataset name
#' @param uploadsDF Reactive value containing dataframe having name, size, type and datapath of the uploaded files
#' @export
#' @return Returns a Reactive value containing seurat object
#'
read_10x <- function(input, output, session, dataset_name, uploadsDF) {

    print("line 44")
    #Name of the temporary directory where files are uploaded
    tempdirname <- dirname(uploadsDF$datapath[1])
    print("line 55")
    print(tempdirname)
    print(uploadsDF$datapath[1])

    # default file names are 0.tsv, 1.tsv and 0.mtx
    # Use the names coloumn in uploadsDF to rename the files to genes.tsv, barcodes.tsv and matrix.tsv
    # else, files wont be recognized by read10x function
    for (i in 1:nrow(uploadsDF)) {
      file.rename(uploadsDF$datapath[i],
                  paste0(tempdirname, "/", uploadsDF$name[i]))
    }

    #files can now be read
    seurat.data <- tryCatch({
      seurat.data <- Read10X(data.dir = tempdirname)
      seurat<- make_seurat_obj(seurat.data, dataset_name)
    },
    error = function(cond) {
      sendSweetAlert(
        session = session,
        title = "Data format error",
        text = "Ensure that 10X Genomics data were appropriately chosen:
              accepted files are either .h5 file or matrix.mtx, genes.tsv/features.tsv and barcodes.tsv files ",
        type = "error"
      )
      # return()
    })

  return(seurat)
}

#' Function to create seaurat object
#'
#' @param seurat.data Reactive value containing
#'
#' @export
#' @return Returns a Reactive value containing list of downsampled seurat objects
#'  with reassigned cluster labels that correspoond across downsampled subsets

make_seurat_obj <- function(seurat.data, dataset_name) {

  #Create first seurat object
  seurat = CreateSeuratObject(
    counts = seurat.data,
    project = dataset_name,
    min.cells = 3,
    min.features = 200
  )

  return(seurat)

}
