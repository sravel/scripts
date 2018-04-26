
if (exists("fileRData")){
  AnalysisDir = datapath
}else{
  AnalysisDir = '~'
}

shinyDirChoose(
  input,'dirIn',
  roots = c(home = homeDir),
  filetypes = c('', "png", "PNG","jpg","JPG","jpeg","JPEG")
)

dirIn <- reactive(input$dirIn)

output$dirInAnalysis <- renderText({
  updateDirAnalysis()
  parseDirPath(c(home = homeDir), dirIn())
})

updateDirAnalysis <- eventReactive(input$dirIn,{
  home <- normalizePath(homeDir)
  datapath <<- file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
  exitStatusAnalysis <<- list(codeAnalysis=-1, messAnalysis = "NULL", errAnalysis = "NULL")
})

#### For RData file
RDataIn <- reactive({
  infile <- input$datafile
  if (is.null(infile)) {
    # User has not uploaded a file yet
    return(NULL)
  }else{
    parseDirPath(c(home = homeDir), infile)
    print(infile)
   }
  })
  

output$RDataText <- renderPrint({
  RDataIn()
})




