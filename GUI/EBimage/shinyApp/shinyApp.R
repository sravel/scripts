library(shiny)
library(shinyFiles)
server <- function(input, output){

  # dir
  shinyDirChoose(input, 'dir', roots = c(home = '~'), filetypes = c('', 'txt'))
  dir <- reactive(input$dir)
  output$dir <- renderPrint(dir())

  # path
  path <- reactive({
    home <- normalizePath("~")
    file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
  })

  # files
  output$files <- renderPrint(list.files(path()))

}





ui <- fluidPage(
  HTML("<h1>EBImage analysis tools</h1>"),
  # input frame

  tags$div(
		tags$b("Counts Directory Input:")
		),
	fluidRow(
		column(1,
			shinyFiles::shinyDirButton(id = 'data_folder', label = 'Select Data Folder', title = 'Please select a folder', FALSE, class = "btn-info")
			),
		column(2, verbatimTextOutput("data_folder_show"))
			),

  #output Frame
  fluidRow(
    verbatimTextOutput("fileCalibration"),
    verbatimTextOutput("run")
  )
)

















shinyApp(ui = ui, server = server)
