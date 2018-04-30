tabItem(
  # Tab for calibration input/output
  tabName = "tabCalibration",
  fluidRow(
    box(
      title = "Calibration Input", status = "primary",solidHeader = TRUE, collapsible = TRUE, width = 12,
      column(4,
             tags$div(class = "calibrationTXT", "Calibration input directory must include sub-folders:",  tags$br(),
                      tags$ul(
                        tags$li("limbe"),
                        tags$li("background"),
                        tags$li("lesion")
                      )
             )
      ),
      column(8,
             shinyDirButton(id = 'dir', label = 'Select Data Folder', title = 'Please select a folder', FALSE, class = "btn-info"),
             bsPopover(id = "dir", "Select Input folder", "the input folder must have sub-folders", trigger="hover", options = NULL),
             conditionalPanel(
               condition = "input.dir", br(),
               verbatimTextOutput("dir", placeholder = FALSE),
               actionButton("runButton", "run calibration!")
             )
      )
    ),
    conditionalPanel(
      condition = "output.code == 0",
      box(
        title = "ERROR", status = "danger",solidHeader = TRUE,
        uiOutput("err")
      )
    )
  ),
  fluidRow(
    conditionalPanel(
      condition = "output.code == 1",
      box(
        title = "Calibration output", status = "success",solidHeader = TRUE, width = 12,
        
        fluidRow(
          column(4,
                 verbatimTextOutput("mess",placeholder = FALSE)
          )
        ),
        fluidRow(
          column(4,
                 plotOutput("img")
          ),
          column(8,
                 tableOutput("table")
          )
        )
      )
    )
  )
)
