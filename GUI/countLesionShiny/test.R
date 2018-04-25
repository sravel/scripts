# library(shiny)
# library(shinyFiles)
# 
# # Define UI for application that draws a histogram
# ui <- fluidPage( # Application title
#   mainPanel(
#     shinyDirButton("dir", "Input directory", "Upload"),
#     verbatimTextOutput("dir", placeholder = TRUE)  # added a placeholder
#   ))
# 
# server <- function(input, output) {
#   shinyDirChoose(
#     input,
#     'dir',
#     roots = c(home = '~'),
#     filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
#   )
#   
#   dir <- reactive(input$dir)
#   output$dir <- renderText({  # use renderText instead of renderPrint
#     parseDirPath(c(home = '~'), dir())
#   })
#   
#   observeEvent(ignoreNULL = TRUE,
#                eventExpr = {
#                  input$dir
#                },
#                handlerExpr = {
#                  home <- normalizePath("~")
#                  datapath <<-
#                    file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
#                })
# }
# 
# # Run the application
# shinyApp(ui = ui, server = server)
# 
# library(shiny)
# library(shinyFiles)
# 
# ui <-   fluidPage(
#   
#   shinyFilesButton("blue", "blue band" ,
#                    title = "Please select a folder:",
#                    buttonType = "default", class = NULL, multiple = F)
# )
# 
# 
# server <- function(input,output,session){
#   
#   
#   volumes = getVolumes() 
#   observe({
#     
#     shinyFileChoose(input, "blue", roots = c(home = "~"), session = session)
#     if(!is.null(input$blue)){
#       myOutput1 <<- parseFilePaths(c(home = '~'),input$blue)
#       print(myOutput1)
#       print(myOutput1$datapath)
#       # myblue <- path.expand(myOutput1) #myblue isthen my file path that I can use in my function
#     }
#   })
# }
# 
# shinyApp(ui = ui, server = server)
# 
# 
# 



ui <- shinyUI(
  fileInput("inFile", label="Choose a file", multiple=F)
)

server <- shinyServer(function(input, output, session) {
  values <- reactiveValues()
  
  dat <- reactive({
    if (is.null(inFile$datapath)) {
      dat <- read.file("Samples.RData")
      values$file_name = "Samples.RData"
      values$file_type = "RData"
      values$file_size = file.size("Samples.RData")
      values$file_path = "Samples.RData"
    } else {
      print(inFile$datapath)
      dat <- read.file(inFile$datapath)
      values$file_name = inFile$name
      values$file_size = inFile$size
      values$file_type = inFile$type
      values$file_path = inFile$datapath
    }
  })
})

shinyApp(ui=ui, server=server)
