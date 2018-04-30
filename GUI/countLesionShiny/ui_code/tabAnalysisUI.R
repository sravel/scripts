tabItem(
  # Tab for calibration input/output
  tabName = "tabAnalysis",
  fluidRow(
    box(
      title = "Analysis Input", status = "primary",solidHeader = TRUE, collapsible = TRUE, width = 12,
      column(width = 4,
             fluidRow(
               shinyDirButton(id = 'dirIn', label = 'Select images folder', title = 'Please select a folder', FALSE, class = "btn-info"),
               conditionalPanel(
                 condition = "input.dirIn", br(),
                 verbatimTextOutput("dirInAnalysis", placeholder = TRUE)
               )
             ),             
             fluidRow(
               shinyDirButton(id = 'dirOut', label = 'Select output results folder', title = 'Please select a folder', FALSE, class = "btn-info"),
               conditionalPanel(
                 condition = "input.dirOut", br(),
                 verbatimTextOutput("dirOutAnalysis", placeholder = TRUE)
               )
             ),
             fluidRow(
               shinyFilesButton('files', label='Load Rdata build in Calibration', title='Please select Rdata file', multiple=T, class = "btn-info"),
               verbatimTextOutput('fileRdata', placeholder = TRUE),
               
               conditionalPanel(
                 condition = "input.dirIn && output.fileRdata && output.codeValidationInt == 1", br(),
                 actionButton("runButtonAnalysis", "run Analysis!")
               )
             )
      ),
      column(width = 4,
             h4("Leaf parameters:"),
             numericInput("leaf_min_size", "Leaf min size:", value = 1000, width = "150px"),
             numericInput("leaf_border_size", "Leaf border size:", value = 3, width = "150px")
             # verbatimTextOutput("value")
             
      ),
      column(width = 4,
             h4("Lesion parameters:"),
             numericInput("lesion_min_size", "Lesion min size:", value = 10, width = "150px"),
             numericInput("lesion_border_size", "Lesion border size:", value = 3, width = "150px"),
             selectInput("lesion_color", label = p("Lesion color output"), 
                         choices = list("Black" = 1, "White" = 2), 
                         selected = 1, width = "150px")
      )
    )
  ),
  fluidRow(
    conditionalPanel(
      condition = 'output.codeValidationInt==0',
      box(
        title = "Warning", status = "warning",solidHeader = TRUE,
        uiOutput("warning")
      )
    )
  ),
  fluidRow(
    conditionalPanel(
      condition = "output.analysisFinish==1",
      box(
        title = "Analysis output", status = "success",solidHeader = TRUE, width = 12,
        verbatimTextOutput("analysisFinish",placeholder = FALSE),
        div(style="display:inline-block",uiOutput("infoButton")),
        DT::dataTableOutput("table2")
      )
    )
  )
  
)

