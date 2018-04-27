if (exists("fileRData")){
  # AnalysisDir = '~'
  AnalysisDir = '~/Bayer/AnalyseImagesV4'
  lda1 <<- load(file = fileRData, envir = .GlobalEnv)
}else{
  # AnalysisDir = '~'
  AnalysisDir = '~/Bayer/AnalyseImagesV4'
  lda1 <<- NULL
}


shinyDirChoose(
  input,'dirIn',
  roots = c(home = AnalysisDir),
  filetypes = c('', "png", "PNG","jpg","JPG","jpeg","JPEG")
)

dirIn <- reactive(input$dirIn)

output$dirInAnalysis <- renderText({
  updateDirAnalysis()
  parseDirPath(c(home = AnalysisDir), dirIn())
})

updateDirAnalysis <- eventReactive(input$dirIn,{
  home <- normalizePath(AnalysisDir)
  datapathAnalysis <<- file.path(home, paste(unlist(dirIn()$path[-1]), collapse = .Platform$file.sep))
  exitStatusAnalysis <<- list(codeAnalysis=-1, messAnalysis = "NULL", errAnalysis = "NULL")
  if (exists("fileRData")){
    # AnalysisDir = '~'
    AnalysisDir = '~/Bayer/AnalyseImagesV4'
    lda1 <<- load(file = fileRData, envir = .GlobalEnv)
  }else{
    # AnalysisDir = '~'
    AnalysisDir = '~/Bayer/AnalyseImagesV4'
    lda1 <<- NULL
  }
})

#### For RData file

roots =  c(wd = AnalysisDir)

shinyFileChoose(input, 'files', 
                roots=roots, 
                filetypes=c('', 'rdata' , 'RData'))

output$fileRdata <- renderText({
  RDataIn()
})


RDataIn <- eventReactive(input$files,{
  infile <<- as.character(parseFilePaths(roots, input$files)$datapath)
  if (is.null(infile)) {
    # User has not uploaded a file yet
    return(NULL)
  }else{
    load(file = infile, envir = .GlobalEnv)
    return(infile)
  }
})

##### Rdata already load (same session)
output$RDataLoad <- renderText({
  updateDirAnalysis()
  result()
  if (exists("lda1")){
    1
  } else{
    0
  }
  
})



validate_INT <- function(input,name) {
  if (!is.numeric(input)) {
    paste("Please input a number >= 0 for ",name," !")
  } else {
    NULL
  }
}

values <- reactive({
  if(!is.numeric(input$leaf_min_size) || (input$leaf_min_size <= 0)){
    return(list(codeAnalysis=0, warning=HTML("Please input a number >= 0 for <b>leaf_min_size</b> !")))
  }  
  if(!is.numeric(input$leaf_border_size) || (input$leaf_border_size <= 0)){
    return(list(codeAnalysis=0, warning=HTML("Please input a number >= 0 for <b>leaf_border_size</b> !")))
  }  
  if(!is.numeric(input$lesion_min_size) || (input$lesion_min_size <= 0)){
    return(list(codeAnalysis=0, warning=HTML("Please input a number >= 0 for <b>lesion_min_size</b> !")))
  }  
  if(!is.numeric(input$lesion_border_size) || (input$lesion_border_size <= 0)){
    return(list(codeAnalysis=0, warning=HTML("Please input a number >= 0 for <b>lesion_border_size</b> !")))
  }
  
  return(list(codeAnalysis=1,
              leaf_min_size=as.numeric(input$leaf_min_size),
              leaf_border_size=as.numeric(input$leaf_border_size),
              lesion_min_size=as.numeric(input$lesion_min_size),
              lesion_border_size=as.numeric(input$lesion_border_size)
  ))
})
output$codeAnalysis <- renderText({
  values = values() 
  values$codeAnalysis
})
output$warning <- renderUI({
  values = values() 
  values$warning
})

output$value <- renderPrint({ 
  values = values() 
  # if (values$codeAnalysis == 1){
  values()
  # }
})

outputOptions(output, 'codeAnalysis', suspendWhenHidden = FALSE)
outputOptions(output, 'warning', suspendWhenHidden = FALSE)

########### run analysis

resultAnalysis <- eventReactive(input$runButtonAnalysis,{
  color <- as.numeric(input$lesion_color)
  values = values() 
  pathResult <- "/media/sebastien/Bayer/AnalyseImagesV4/Exemple2/Result" 

  withProgress(message = 'Making Analysis, please wait\n', value = 0, {
  # Increment the progress bar, and update the detail text.
  incProgress(1/6, detail = "analysis start")
    toto <<- analyseFiles(fileRdata=infile,
                           pathResult=pathResult,
                           pathImages=datapathAnalysis,
                           onefileImage=NA, ## analyse du répertoire complet
                           leafMinSize=values$leaf_min_size,
                           leafBorderSize=values$leaf_border_size,
                           lesionBorderSize=values$lesion_border_size,
                           lesionMinSize=values$lesion_min_size,
                           colorLesion=color)
  })
})
  

output$analysisFinish <- renderPrint({
  resultAnalysis()
  toto
})
  
  # analyse.image(path.sample=path.sample,
  #               path.result=path.result,
  #               path.image=path.image,
  #               file.image=NA, ## analyse du répertoire complet
  #               surface.feuille.mini=surface.feuille.mini,
  #               bordure.feuille=bordure.feuille,
  #               bordure.lesion=bordure.lesion,
  #               surface.lesion.mini=surface.lesion.mini,
  #               couleur.lesion=couleur.lesion)
  