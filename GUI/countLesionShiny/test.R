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



library(shiny)
library(shinyFiles)


ui <- fluidPage( 
  shinyFilesButton('files', label='File select', title='Please select a file', multiple=T) ,
  verbatimTextOutput('rawInputValue'),
  verbatimTextOutput('filepaths') ,
  downloadButton("downloadFiles", "Download Files")
)

server <- function(input, output) {
  
  roots =  c(wd = '~')
  
  shinyFileChoose(input, 'files', 
                  roots =  roots, 
                  filetypes=c('', 'txt' , 'gz' , 'md5' , 'pdf' , 'fasta' , 'fastq' , 'aln'))
  
  output$rawInputValue <- renderPrint({str(input$files)})
  
  output$filepaths <- renderPrint({parseFilePaths(roots, input$files)})
  
  output$downloadFiles <- downloadHandler(
    filename = function() {
      as.character(parseFilePaths(roots, input$files)$name)
    },
    content = function(file) {
      fullName <- as.character(parseFilePaths(roots, input$files)$datapath)
      file.copy(fullName, file)
    }
  )
}

shinyApp(ui = ui , server = server)