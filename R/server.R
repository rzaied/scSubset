#' Initialize Server
#'
#' @export
#' @return None

server <- function(input, output, session) {
  observeEvent(input$submitButton,
               {
                 #assign uploaded file to seurat.data object
                 print("here")
                 dataset1_name<-"dataset1"
                 dataset2_name<-"dataset2"

                 if (is.null(input$file1)) {
                   print("empy files")
                   shinyalert("Please upload at least one Chromium 10X scRNA-seq dataset")

                 }
                 else if (grepl(".h5", input$file1$name[1], fixed = TRUE)) {

                   path=dirname(input$file1$datapath[1])
                   print(path)
                   #assign dataset name (for "project" in seurat)
                   #assign manually to avoid errors

                   seurat1 <-
                     callModule(read_h5,
                                "ui", dataset1_name, path)
                   print("server line 27")

                 }

                 else if (grepl(".tsv", input$file1$name[1], fixed = TRUE) ||
                          grepl(".mtx", input$file1$name[1], fixed = TRUE)) {
                   #uploadsDF is a data.frame where columns= name, size, type and datapath of the uploaded files
                   uploadsDF <- input$file1

                   seurat1 <-
                     callModule(read_10x,
                                "ui", dataset1_name, uploadsDF )
                   print("server line 40")

                 }

                 else {
                   print("wrong file")
                   shinyalert(
                     "Please use the required file formats
                     (.h5 file or matrix.mtx, genes.tsv/features.tsv and barcodes.tsv files )"
                   )
                 }


                 print("line 45")
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

                 #assign # of genes per cluster to test for
                 selectedNumGenes <- as.numeric(input$numGenes)

                 #assign resolution
                 res <- input$resolution


                 #assigns custom name for second dataset, only if integration option was selected
                 if (input$integrationChoice) {
                   if (is.null(input$file2)) {
                     print("empy files")
                     shinyalert("Please upload a second Chromium 10X scRNA-seq dataset for integration analysis")

                   }

                   #if the second file is a .h5 file
                   else if (grepl(".h5", input$file2$name[1], fixed = TRUE)) {

                     path=dirname(input$file2$datapath[1])
                     print(path)
                     #assign dataset name (for "project" in seurat)
                     #assign manually to avoid errors
                     seurat2 <-
                       callModule(read_h5,
                                  "ui", dataset2_name, path)
                     print("server line 27")

                    }

                   #if the second file is of .mtx, .genes, .barcodes formats
                   else if (grepl(".tsv", input$file2$name[1], fixed = TRUE) ||
                            grepl(".mtx", input$file2$name[1], fixed = TRUE)) {

                     uploadsDF <- input$file2

                     seurat2 <-
                       callModule(read_10x,
                                  "ui", dataset2_name, uploadsDF )
                     print("server line 40")

                   }

                 }

                 # update progress bar value
                 print("line 77")
                 update_modal_progress(0.03)

                 print("line 105")


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
                   analysis_type = "integration"
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
                                seurat1,
                                seurat2,
                                mito,
                                res)

                   # update progress bar value
                   update_modal_progress(0.10)

                   #place holder table
                   combinedConsMarkersTable <- 0
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
                     combinedConsMarkersTable <-
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
                     dataset1_name,
                     dataset2_name,
                     selectedNumGenes
                   )

                   #call server, pass it the UI,
                   callModule(
                     computeCost,
                     "costTab",
                     analysis_type,
                     combinedConsMarkersTable,
                     combinedDEgenesTable,
                     seuratObjectsList,
                     costPerMil,
                     depthPerCell
                   )

                   update_modal_progress(0.60)



                   # update progress bar value

                   rm(seurat.data)
                   rm(seurat2)
                   rm(seurat1)
                   rm(seuratObjectsList)
                   update_modal_progress(1)


                 }

                 #if no integration:
                 else {
                   analysis_type = "Single dataset"
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
                                seurat1,
                                mito,
                                res,
                                dataset1_name)


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

                  print("start of find markers")
                  print(seuratObjectsList)
                   combinedMarkersTable <-
                     callModule(findMarkerGenes,
                                "markerGenesTab",
                                seuratObjectsList,
                                selectedNumGenes)
                   #place holder table value for conserved marker genes (not actually used by function for single datasets)
                   combinedConsMarkersTable<-0
                   #call server, pass it the UI,
                   callModule(
                     computeCost,
                     "costTab",
                     analysis_type,
                     combinedConsMarkersTable,
                     combinedMarkersTable,
                     seuratObjectsList,
                     costPerMil,
                     depthPerCell
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
                 dataset1_name = "PBMC"
                 seurat.data = Read10X(data.dir = sc_example_data)
                 seurat1 = CreateSeuratObject(
                   counts = seurat.data,
                   project = dataset1_name,
                   min.cells = 3,
                   min.features = 200
                 )
                 mito = "^MT-"
                 res = 0.4
                 selectedNumGenes = 5
                 costPerMil = 8
                 depthPerCell = 50000
                 #libraryPrepCost=2000
                 analysis_type = "Single dataset"

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
                              seurat1,
                              mito,
                              res,
                              dataset1_name)

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
                 #place holder table
                 combinedConsMarkersTable <- 0
                 #call server, pass it the UI,
                 callModule(
                   computeCost,
                   "costTab",
                   analysis_type,
                   combinedConsMarkersTable,
                   combinedMarkersTable,
                   seuratObjectsList,
                   costPerMil,
                   depthPerCell
                 )




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
