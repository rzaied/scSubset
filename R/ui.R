#' Initiate GUI
#'
#'
#' @export
#' @return None

# Create user interface -------------------------------
ui <- tagList(
  navbarPage(
    id = "mainPage",
    title = "scSubset",
    theme = shinytheme('flatly'),
    tabPanel("About", icon = icon("bars"), fluidRow(column(
      12,
      wellPanel(style = "background-color: #fff; border-color: #2c3e50; height: 1000px;",
                includeHTML("inst/extdata/about.html"))
    ))),


    tabPanel(
      "Data Upload",
      icon = icon("upload"),

      useShinyjs(),
      useShinyalert(),

      #use_busy_spinner(spin = "fading-circle"),
      fluidRow(column(
        12,
        wellPanel(
          style = "background-color: #fff; border-color: #2c3e50; height: 900px;",
          #Load example data option
          actionButton(
            inputId = 'exampleButton',
            label = 'Load Example 10k PBMC Data',
            class = 'btn btn-primary btn-sm',
            style = "color: white;background-color: grey"
          ),

          br(),
          hr(),
          # home section
          checkboxInput(
            inputId = "integrationChoice",
            label = strong("Preform data integration")
          ),

          # Input: name dataset1 ----
          textInput(
            inputId = "dataset1",
            value = "PBMC",
            label = "Dataset Name"
          ),

          # Input: Select a file ----
          fileInput(
            "file1",
            "Upload genes.tsv, barcodes.tsv, and matrix.mtx files",
            multiple = TRUE,
            accept = c(
              "text/csv/application/zip",
              "text/comma-separated-values,text/plain",
              ".csv",
              ".tsv",
              ".mtx",
              ".gz"
            )
          ),
          # Display only if the smoother is checked
          conditionalPanel(
            condition = "input.integrationChoice == true",

            # Input: name dataset1 ----
            textInput(
              inputId = "dataset2",
              value =
                "PBMC2",
              label = "Dataset 2 Name"
            ),

            # Input: Select a file ----
            fileInput(
              "file2",
              "Upload genes, barcodes, and matrix files",
              multiple = TRUE,
              accept = c(
                "text/csv/application/zip",
                "text/comma-separated-values,text/plain",
                ".csv",
                ".tsv",
                ".mtx",
                ".gz"
              )
            ),

            checkboxInput(
              inputId = "conservedMarkersChoice",
              label = strong("Find conserved markers")
            ),
            HTML(
              "*Identifying conserved markers approximately doubles processing time"
            )
          ),
          br(),

          textInput(
            inputId = "organism",
            value = "^MT-",
            label = "Mitochondrial genes pattern"
          ),

          textInput(
            inputId = "numGenes",
            value = "5",
            label = "Number of highly significant genes to test per cluster, per subset"
          ),


          sliderInput(
            inputId = "resolution",
            label = "Clustering resolution",
            min = 0.1,
            max = 1.0,
            value = 0.5,
            width = 300
          ),
          HTML(
            "*Higher resolutions yield more clusters and require longer processing time"
          ),

          br(),
          hr(),
          actionButton(
            inputId = 'submitButton',
            label = 'Submit',
            class = 'btn btn-primary btn-sm',
            style = "color: white;
                           background-color: grey"
          ),
          actionButton(
            inputId = 'resetButton',
            label = 'Reset',
            icon = icon("refresh"),
            class = 'btn btn-primary btn-sm',
            style = "color: white;
                           background-color: grey"
          )
        )
      ))
    )
  )
)

