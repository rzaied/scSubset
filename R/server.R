#' Initialize Server
#'
#' @export
#' @return None

server <- function(input, output, session) {
  observeEvent(input$submitButton,
               {
                 #code to assign uploaded file to seurat.data object
                 print("here")

                 if (is.null(input$file1)) {
                   print("empy files")
                   shinyalert("Please upload at least one Chromium 10X scRNA-seq dataset")

                 }
                 else if (grepl(".h5", input$file1$name[1], fixed = TRUE)) {
                   seurat.data <- reactive({
                     print("line 13")
                     print(input$file1$name[1])
                     print(input$file1$datapath[1])
                     tempdirname <- dirname(input$file1$datapath[1])
                     print("line 22")
                     print(tempdirname)
                     tempdirname <- paste0(tempdirname, "/", "0.h5")
                     tempdirname
                     seurat.data <- tryCatch({
                       seurat.data <- Read10X_h5(tempdirname)
                     },
                     error = function(cond) {
                       sendSweetAlert(
                         session = session,
                         title = "Data format error",
                         text = "Ensure that 10X Genomics data were appropriately chosen:
              accepted files are either .h5 file or matrix.mtx, genes.tsv/features.tsv and barcodes.tsv files ",
                         type = "error"
                       )
                       return()
                     })
                   })
                 }

                 else if (grepl(".tsv", input$file1$name[1], fixed = TRUE) ||
                          grepl(".mtx", input$file1$name[1], fixed = TRUE)) {
                   seurat.data <- reactive({
                     req(input$file1)
                     print("line 44")
                     # uploadsDF is a data.frame where columns= name, size, type and datapath of the uploaded files
                     uploadsDF <- input$file1

                     # Name of the temporary directory where files are uploaded
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
                     },
                     error = function(cond) {
                       sendSweetAlert(
                         session = session,
                         title = "Data format error",
                         text = "Ensure that 10X Genomics data were appropriately chosen:
              accepted files are either .h5 file or matrix.mtx, genes.tsv/features.tsv and barcodes.tsv files ",
                         type = "error"
                       )
                       return()
                     })

                   })
                 }

                 else {
                   print("wrong file")
                   shinyalert(
                     "Please adhere to required file formats
                     (.h5 file ormatrix.mtx, genes.tsv/features.tsv and barcodes.tsv files )"
                   )
                 }


                 print("line 45")
                 #to restric tab creation

                 #require file input to proceed
                 shiny::req(input$file1)

                 # show the spinner
                 show_spinner() # show the spinner
                 show_modal_progress_line() # show the modal window

                 ##assign mito gene pattern
                 mito <- input$organism
                 print(mito)

                 #assign costPerMil
                 costPerMil <- as.numeric(input$cost)
                 print(costPerMil)

                 #assign costPerMil
                 depthPerCell <- as.numeric(input$depth)
                 print(depthPerCell)

                 #assign top # of genes to test for
                 selectedNumGenes <- as.numeric(input$numGenes)

                 #assign top # of genes to test for
                 libraryPrepCost <- as.numeric(input$libraryPrep)

                 #assign resolution
                 res <- input$resolution

                 #assign dataset name (for "project" in seurat)
                 dataset <- input$dataset1

                 #assigns custom name for second dataset, only if integration option was selected
                 if (input$integrationChoice) {
                   if (is.null(input$file2)) {
                     print("empy files")
                     shinyalert("Please upload a second Chromium 10X scRNA-seq dataset for integration analysis")

                   }

                   else if (grepl(".h5", input$file2$name[1], fixed = TRUE)) {
                     seurat.data2 <- reactive({
                       print("line 117")
                       tempdirname <-
                         paste0(dirname(input$file1$datapath[1]), "/", "0.h5")
                       tempdirname
                       seurat.data2 <- tryCatch({
                         seurat.data2 <- Read10X_h5(tempdirname)
                       },
                       error = function(cond) {
                         sendSweetAlert(
                           session = session,
                           title = "Data format error",
                           text = "Ensure that 10X Genomics data were appropriately chosen:
              accepted files are either .h5 file or matrix.mtx, genes.tsv/features.tsv and barcodes.tsv files ",
                           type = "error"
                         )
                         return()
                       })
                     })
                   }

                   else if (grepl(".tsv", input$file2$name[1], fixed = TRUE) ||
                            grepl(".mtx", input$file2$name[1], fixed = TRUE)) {
                     seurat.data2 <- reactive({
                       print("line 44")
                       # uploadsDF is a data.frame where columns= name, size, type and datapath of the uploaded files
                       uploadsDF <- input$file2
                       # Name of the temporary directory where files are uploaded
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
                       seurat.data2 <- tryCatch({
                         seurat.data2 <- Read10X(data.dir = tempdirname)
                       },
                       error = function(cond) {
                         sendSweetAlert(
                           session = session,
                           title = "Data format error",
                           text = "Ensure that 10X Genomics data were appropriately chosen:
              accepted files are either .h5 file or matrix.mtx, genes.tsv/features.tsv and barcodes.tsv files ",
                           type = "error"
                         )
                         return()
                       })

                     })
                   }

                   dataset2 <- input$dataset2
                   seurat2 = CreateSeuratObject(
                     counts = seurat.data2(),
                     project = dataset2,
                     min.cells = 3,
                     min.features = 200
                   )

                 }

                 # update progress bar value
                 print("line 77")
                 update_modal_progress(0.03)

                 print("line 105")
                 #Create first seurat object
                 seurat = CreateSeuratObject(
                   counts = seurat.data(),
                   project = dataset,
                   min.cells = 3,
                   min.features = 200
                 )

                 # update progress bar value

                 update_modal_progress(0.05)

                 # Append cost tab UI ----
                 appendTab(
                   inputId = "mainPage",
                   tabPanel(
                     id = "costTab",
                     title = "Cost estimate",
                     icon = icon("dollar-sign"),
                     computeCostUI("costTab"),
                     value = "costTab"
                   )
                 )

                 #only preform integration if option was selected from menu
                 if (input$integrationChoice) {
                   analysis_type="integration"
                   req(input$file2)

                   #callModule(scIntegrate)
                   # Append clustering tab UI ----
                   appendTab(
                     inputId = "mainPage",
                     tabPanel(
                       id = "clusteringTab",
                       title = "Clustering solutions",
                       icon = icon("braille"),
                       scClusteringUI("clusteringTab"),
                       value = "clusteringTab"
                     )
                   )

                   print("line 140")

                   #The callModule function is used within the Server function
                   #to call functions that creat reactive output
                   seuratObjectsList <-
                     callModule(scIntegrate,
                                "clusteringTab",
                                seurat,
                                seurat2,
                                mito,
                                res)

                   # update progress bar value
                   update_modal_progress(0.10)




                   if (input$conservedMarkersChoice) {
                     # update progress bar value

                     appendTab(
                       inputId = "mainPage",
                       tabPanel(
                         id = "CMGenesTab",
                         title = "Conserved Marker genes",
                         icon = icon("bullseye"),
                         conservedMarkersUI("CMGenesTab"),
                         value = "CMGenesTab"
                       )
                     )


                     # update progress bar value
                     update_modal_progress(0.20)
                     callModule(findConservedMarkers,
                                "CMGenesTab",
                                seuratObjectsList,
                                selectedNumGenes)

                     #  updateTabsetPanel(session, 'mainPage', 'CMGenesTab')
                     # update progress bar value
                     update_modal_progress(0.50)

                   }

                   #DIFFERENTIAL EXPRESSION
                   appendTab(
                     inputId = "mainPage",
                     tabPanel(
                       id = "DEgenesTab",
                       title = "Differential Expression",
                       icon = icon("arrows-alt-v"),
                       findDEgenesUI("DEgenesTab"),
                       value = "DEgenesTab"
                     )
                   )


                   print("line 162")
                   combinedDEgenesTable <- callModule(
                     findDEgenes,
                     "DEgenesTab",
                     seuratObjectsList,
                     dataset,
                     dataset2,
                     selectedNumGenes
                   )

                   #call server, pass it the UI,
                   callModule(computeCost,
                              "costTab", analysis_type, combinedDEgenesTable, seuratObjectsList, costPerMil, depthPerCell, libraryPrepCost
                   )

                   update_modal_progress(0.60)



                   # update progress bar value

                   #remove unwanted obj once integration done and no conserved markers chosen
                   #rm(seurat.combined)
                   rm(seurat.data)
                   rm(seurat2)
                   rm(seurat)
                   rm(seuratObjectsList)
                   update_modal_progress(1)


                 }

                 #if no integration:
                 else {
                   analysis_type="Single dataset"
                   # update progress bar value
                   print("line 218")
                   update_modal_progress(0.10)
                   #when only one dataset is selected

                   # Append clustering tab UI and swap ----
                   appendTab(
                     inputId = "mainPage",
                     tabPanel(
                       id = "clusteringTab",
                       title = "Clustering solutions",
                       icon = icon("braille"),
                       scClusteringUI("clusteringTab"),
                       value = "clusteringTab"
                     )
                   )

                   print("line 233")

                   appendTab(
                     inputId = "mainPage",
                     tabPanel(
                       id = "varGenesTab",
                       title = "Variable Genes",
                       icon = icon("chart-bar"),
                       varGenesUI("varGenesTab"),
                       value = "varGenesTab"
                     )
                   )

                   print("line 245")

                   print("line 248")
                   seuratObjectsList <-
                     callModule(scClustering,
                                "clusteringTab",
                                seurat,
                                mito,
                                res,
                                dataset)


                   callModule(scVarGenes, "varGenesTab", seuratObjectsList)



                   # update progress bar value

                   update_modal_progress(0.60)

                   #Find marker genes for single dataset
                   appendTab(
                     inputId = "mainPage",
                     tabPanel(
                       id = "markerGenesTab",
                       title = "Marker genes",
                       icon = icon("bullseye"),
                       findMarkerGenesUI("markerGenesTab"),
                       value = "markerGenesTab"
                     )
                   )



                   combinedMarkersTable <-
                     callModule(findMarkerGenes,
                                "markerGenesTab",
                                seuratObjectsList,
                                selectedNumGenes)

                   #call server, pass it the UI,
                   callModule(computeCost,
                              "costTab", analysis_type, combinedMarkersTable, seuratObjectsList, costPerMil, depthPerCell, libraryPrepCost
                   )



                   #update progress bar value
                   update_modal_progress(1)

                 }


                 removeTab("mainPage", "Data Upload", session = session)


                 hide_spinner() # hide the spinner

                 remove_modal_progress() # remove it when done
               })

  observeEvent(input$exampleButton,
               {
                 #example data

                 show_spinner() # show the spinner


                 show_modal_progress_line() # show the modal window


                 #setting up seurat object

                 sc_example_data <- paste0(
                   system.file("extdata",
                               "pbmc_filtered_feature_bc_matrix",
                               package = "scSubset"),
                   "/"
                 )
                 dataset = "PBMC"
                 seurat.data = Read10X(data.dir = sc_example_data)
                 seurat = CreateSeuratObject(
                   counts = seurat.data,
                   project = dataset,
                   min.cells = 3,
                   min.features = 200
                 )
                 mito = "^MT-"
                 res = 0.6
                 selectedNumGenes = 5
                 costPerMil=8
                 depthPerCell=50000
                 libraryPrepCost=2000
                 analysis_type="Single dataset"

                 update_modal_progress(0.10)
                 #running example data UI

                 # Append cost tab UI ----
                 appendTab(
                   inputId = "mainPage",
                   tabPanel(
                     id = "costTab",
                     title = "Cost estimate",
                     icon = icon("dollar-sign"),
                     computeCostUI("costTab"),
                     value = "costTab"
                   )
                 )

                 # Append clustering tab UI and swap ----
                 appendTab(
                   inputId = "mainPage",
                   tabPanel(
                     id = "clusteringTab",
                     title = "Clustering solutions",
                     icon = icon("braille"),
                     scClusteringUI("clusteringTab"),
                     value = "clusteringTab"
                   )
                 )

                 appendTab(
                   inputId = "mainPage",
                   tabPanel(
                     id = "varGenesTab",
                     title = "Variable Genes",
                     icon = icon("chart-bar"),
                     varGenesUI("varGenesTab"),
                     value = "varGenesTab"
                   )
                 )

                 seuratObjectsList <-
                   callModule(scClustering,
                              "clusteringTab",
                              seurat,
                              mito,
                              res,
                              dataset)
                 #example data functions
                 #objList<- callModule(exampleDataClust, "clusteringTab")



                 callModule(scVarGenes, "varGenesTab", seuratObjectsList)

                 update_modal_progress(0.60)
                 #Find marker genes for single dataset
                 appendTab(
                   inputId = "mainPage",
                   tabPanel(
                     id = "markerGenesTab",
                     title = "Marker genes",
                     icon = icon("bullseye"),
                     findMarkerGenesUI("markerGenesTab"),
                     value = "markerGenesTab"
                   )
                 )

                 combinedMarkersTable <-
                   callModule(findMarkerGenes,
                              "markerGenesTab",
                              seuratObjectsList,
                              selectedNumGenes)

                 #call server, pass it the UI,
                 callModule(computeCost,
                            "costTab", analysis_type, combinedMarkersTable, seuratObjectsList, costPerMil, depthPerCell, libraryPrepCost
                 )

                 #callModule(exampleDataMarkers, "markerGenesTab")


                 #callModule(exampleDataVAR, "varGenesTab", objList)


                 # updateTabsetPanel(session, 'mainPage', 'markerGenesTab')
                 #updateTabsetPanel(session, 'mainPage', 'clusteringTab')


                 #hide the uploads tab if example data selected
                 removeTab("mainPage", "Data Upload", session = session)
                 hide_spinner() # hide the spinner

                 update_modal_progress(1)
                 # remove it when done
                 remove_modal_progress()

               })


  #If reset button was selected
  observeEvent(input$resetButton, {
    session$reload()
  })


}
