tabItem(
  # Tab for calibration input/output
  tabName = "tabAnalysis",
  fluidRow(
    box(
      title = "Analysis Input", status = "primary",solidHeader = TRUE, collapsible = TRUE, width = 12,
      column(width = 4,
             shinyDirButton(id = 'dirIn', label = 'Select Images Folder', title = 'Please select a folder', FALSE, class = "btn-info"),
             conditionalPanel(
               condition = "input.dirIn", br(),
               verbatimTextOutput("dirInAnalysis", placeholder = TRUE)
             ),
             br(),
             fileInput("datafile", "Choose Rdata File",
                       accept = c('.RData'), multiple=F
             ),
             verbatimTextOutput("RDataText", placeholder = TRUE),

             # conditionalPanel(
             #   condition = "input.datafile", br(),
             #   verbatimTextOutput("RDataText", placeholder = FALSE)
             # ),
             conditionalPanel(
               condition = "input.dirIn && input.datafile", br(),
               actionButton("runButtonAnalysis", "run Analysis!")
             )
             
      ),
      column(width = 4,
             h4("Leaf parameters:"),
             textInput("leaf_min_size", "Leaf min size:", value = 1000, width = "150px"),
             textInput("leaf_border_size", "Leaf border size:", value = 3, width = "150px")
      ),
      column(width = 4,
             h4("Lesion parameters:"),
             textInput("lesion_min_size", "Lesion min size:", value = 10, width = "150px"),
             textInput("leaf_border_size", "Lesion border size:", value = 3, width = "150px"),
             selectInput("lesion_color", label = p("Lesion color output"), 
                         choices = list("Black" = 1, "White" = 2), 
                         selected = 1, width = "150px")
      )
    ),
    conditionalPanel(
      condition = "output.codeAnalysis == 0",
      box(
        title = "ERROR", status = "danger",solidHeader = TRUE,
        uiOutput("errAnalysis")
      )
    )
  ),
  fluidRow(
    # conditionalPanel(
    #   condition = "output.codeAnalysis == 1",
      box(
        title = "Analysis output", status = "success",solidHeader = TRUE, width = 12,

        fluidRow(
          column(4,
                 verbatimTextOutput("messAnalysis",placeholder = FALSE)
          )
        # )
      )
    )
  )
)
