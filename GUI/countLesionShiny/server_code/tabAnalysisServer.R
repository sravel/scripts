
if (exists("fileRData")){
  AnalysisDir = datapath
}else{
  AnalysisDir = '~'
}

shinyDirChoose(
  input,'dirIn',
  roots = c(home = homeDir),
  filetypes = c('', "png", "*")
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
shinyFileChoose(
  input, 'files', 
  roots=c(home = AnalysisDir), filetypes=c('.RData')
)

RDataIn <- reactive({
  infile <- input$datafile
  if (is.null(infile)) {
    # User has not uploaded a file yet
    return(NULL)
  }
  objectsLoaded <<- load(input$datafile$name) 
  # the above returns a char vector with names of objects loaded
  df <- eval(parse(text=objectsLoaded[1])) 
  # the above finds the first object and returns it
  return(df)}
  
)

output$RDataText <- renderText({
  RDataIn()
})




