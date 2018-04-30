# library(shiny)
# library(shinyBS)
# shinyApp(
#   ui =
#     fluidPage(
#       sidebarLayout(
#         sidebarPanel(
#           sliderInput("bins",
#                       "Number of bins:",
#                       min = 1,
#                       max = 50,
#                       value = 30),
#           bsTooltip("bins", "The wait times will be broken into this many equally spaced bins",
#                     "right", options = list(container = "body"))
#         ),
#         mainPanel(
#           plotOutput("distPlot")
#         )
#       )
#     ),
#   server =
#     function(input, output, session) {
#       output$distPlot <- renderPlot({
#         
#         # generate bins based on input$bins from ui.R
#         x    <- faithful[, 2]
#         bins <- seq(min(x), max(x), length.out = input$bins + 1)
#         
#         # draw the histogram with the specified number of bins
#         hist(x, breaks = bins, col = 'darkgray', border = 'white')
#         
#       })
#       addPopover(session, "distPlot", "Data", content = paste0("
# 
# Waiting time between ",
#                                                                "eruptions and the duration of the eruption for the Old Faithful geyser ",
#                                                                "in Yellowstone National Park, Wyoming, USA.
# 
# Azzalini, A. and ",
#                                                                "Bowman, A. W. (1990). A look at some data on the Old Faithful geyser. ",
#                                                                "Applied Statistics 39, 357-365.
# "), trigger = 'click')
#     }
# )
# 
# library(shiny)
# library(shinyFiles)
# 
# 
# ui <- fluidPage( 
#   shinyFilesButton('files', label='File select', title='Please select a file', multiple=T) ,
#   verbatimTextOutput('rawInputValue'),
#   verbatimTextOutput('filepaths') ,
#   downloadButton("downloadFiles", "Download Files")
# )
# 
# server <- function(input, output) {
#   
#   roots =  c(wd = '~')
#   
#   shinyFileChoose(input, 'files', 
#                   roots =  roots, 
#                   filetypes=c('', 'txt' , 'gz' , 'md5' , 'pdf' , 'fasta' , 'fastq' , 'aln'))
#   
#   output$rawInputValue <- renderPrint({str(input$files)})
#   
#   output$filepaths <- renderPrint({parseFilePaths(roots, input$files)})
#   
#   output$downloadFiles <- downloadHandler(
#     filename = function() {
#       as.character(parseFilePaths(roots, input$files)$name)
#     },
#     content = function(file) {
#       fullName <- as.character(parseFilePaths(roots, input$files)$datapath)
#       file.copy(fullName, file)
#     }
#   )
# }
# 
# shinyApp(ui = ui , server = server)
# 
# 
# 
library(shiny)
library(DT)

pathImages <- "/home/sebastien/Bayer/AnalyseImagesV4/Images/"

ui <- shinyUI(
  fluidPage(
    actionButton("run", "RUN"),
    DT::dataTableOutput("table")
  )
)

server <- shinyServer(function(input, output) {
  
  files <- eventReactive(input$run, {
    images <- list.files(pathImages, full.names=TRUE)
    return(images)
  })
  
  output$table <- DT::renderDataTable({
    
    dat <-  data.frame(image =sprintf('<img src=%s height="150px"></img>', files()))
    datatable(dat, escape=FALSE)
  })
  
})
shinyApp(ui = ui , server = server)
