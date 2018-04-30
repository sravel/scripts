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
shinyDirChoose(
  input,'dirOut',
  roots = c(home = AnalysisDir),
  filetypes = c('')
)

dirIn <- reactive(input$dirIn)

output$dirInAnalysis <- renderText({
  updateDirAnalysis()
  parseDirPath(c(home = AnalysisDir), dirIn())
})

updateDirAnalysis <- eventReactive(input$dirIn,{
  home <- normalizePath(AnalysisDir)
  datapathAnalysis <<- file.path(home, paste(unlist(dirIn()$path[-1]), collapse = .Platform$file.sep))
  exitStatusAnalysis <<- list(codeValidationInt=-1, codeAnalysis=-1, messAnalysis = "NULL", errAnalysis = "NULL")
})

dirOut <- reactive(input$dirOut)

output$dirOutAnalysis <- renderText({
  updateDirOutAnalysis()
  parseDirPath(c(home = AnalysisDir), dirOut())
})

updateDirOutAnalysis <- eventReactive(input$dirOut,{
  home <- normalizePath(AnalysisDir)
  datapathOutAnalysis <<- file.path(home, paste(unlist(dirOut()$path[-1]), collapse = .Platform$file.sep))
  exitStatusAnalysis <<- list(codeValidationInt=-1, codeAnalysis=-1, messAnalysis = "NULL", errAnalysis = "NULL")
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
    return(list(codeValidationInt=0, warning=HTML("Please input a number >= 0 for <b>leaf_min_size</b> !")))
  }
  if(!is.numeric(input$leaf_border_size) || (input$leaf_border_size <= 0)){
    return(list(codeValidationInt=0, warning=HTML("Please input a number >= 0 for <b>leaf_border_size</b> !")))
  }
  if(!is.numeric(input$lesion_min_size) || (input$lesion_min_size <= 0)){
    return(list(codeValidationInt=0, warning=HTML("Please input a number >= 0 for <b>lesion_min_size</b> !")))
  }
  if(!is.numeric(input$lesion_border_size) || (input$lesion_border_size <= 0)){
    return(list(codeValidationInt=0, warning=HTML("Please input a number >= 0 for <b>lesion_border_size</b> !")))
  }

  return(list(codeValidationInt=1,
              leaf_min_size=as.numeric(input$leaf_min_size),
              leaf_border_size=as.numeric(input$leaf_border_size),
              lesion_min_size=as.numeric(input$lesion_min_size),
              lesion_border_size=as.numeric(input$lesion_border_size)
  ))
})
output$codeValidationInt <- renderText({
  values = values()
  values$codeValidationInt
})
output$warning <- renderUI({
  values = values()
  values$warning
})

output$value <- renderPrint({
  values = values()
  values()
})

########### run analysis

resultAnalysis <- eventReactive(input$runButtonAnalysis,{
  color <- as.numeric(input$lesion_color)
  values = values()

  withProgress(message = 'Making Analysis, please wait\n', value = 0, {
  # Increment the progress bar, and update the detail text.
  incProgress(1/6, detail = "analysis start")
  code <- analyseFiles(fileRdata=infile,
                           pathResult=datapathOutAnalysis,
                           pathImages=datapathAnalysis,
                           onefileImage=NA, ## analyse du répertoire complet
                           leafMinSize=values$leaf_min_size,
                           leafBorderSize=values$leaf_border_size,
                           lesionBorderSize=values$lesion_border_size,
                           lesionMinSize=values$lesion_min_size,
                           colorLesion=color)
  })

  if (code == 0){
    errorMessAnalysis <-tags$div("Error !!!!:")
    exitStatusAnalysis <<-list(codeAnalysis=0, err=errorMess)
  }
  if (code == 1){
    errorMessAnalysis <-tags$div("Error !!!!:")
    exitStatusAnalysis <<-list(codeAnalysis=1, mess="good")
  }
})

#########################################
####   OUTPUT
#########################################

output$analysisFinish <- renderText({
  resultAnalysis()
  exitStatusAnalysis$codeAnalysis
})



output$table2<-DT::renderDataTable({
  resultAnalysis()
  # datapathAnalysisTEST <- "~/Bayer/AnalyseImagesV4/Exemple1/Images/"
  # datapathOutAnalysisTEST <- "~/Bayer/AnalyseImagesV4/Exemple1/"
  addResourcePath("Original",datapathAnalysis) # Images are located outside shiny App
  addResourcePath("LesionColor",datapathOutAnalysis) # Images are located outside shiny App

  addResourcePath("Original",datapathAnalysis) # Images are located outside shiny App
  addResourcePath("LesionColor",datapathOutAnalysis) # Images are located outside shiny App

  LeafNames <- list.files(datapathAnalysis, full.names=FALSE)
  LeafNames2 <- list.files(datapathOutAnalysis, full.names=FALSE, pattern = "*.png")
  LeafTable <<- data.frame(LeafNames,LeafNames2)
  LeafTable <<- within(LeafTable, thumbnail <- paste0("<img src='","Original/",LeafTable$LeafNames,"' height='60'></img>"))
  LeafTable <<- within(LeafTable, thumbnail2 <- paste0("<img src='","LesionColor/",LeafTable$LeafNames2 ,"' height='60'></img>"))

  responseDataFilter2 <- LeafTable[,c(1,2,3,4)]

  displayableData<-DT::datatable(data = as.data.frame(responseDataFilter2, stringAsFactors = FALSE, row.names = NULL),

                                 escape=FALSE,selection="single",rownames=FALSE,colnames=c("Original","LesionColor","Name","toto"),

                                 callback = JS("table.on('dblclick.dt', 'td', function() {
                                                 var row=table.cell(this).index().row;
                                                 Shiny.onInputChange('rows_home',[row, Math.random()])});
                                                 table.on('click.dt', 'td', function() {
                                                 var k=table.cell(this).index().row;
                                                 if(table.rows('.selected').indexes().toArray()!= '' && table.rows('.selected').indexes().toArray() == k){k=-1;}
                                                 Shiny.onInputChange('rows_up_home',[k, Math.random()]);
                                                 Shiny.onInputChange('row_path', table.rows(this).data().toArray());
                                                 });"),

                                 options = list(

                                   paging=TRUE,searching = TRUE,ordering=TRUE,scrollY = 750,scrollCollapse=TRUE,server = TRUE

                                 ))

})

output$infoButton = renderUI({
  s = input$table2_rows_selected # Row number of selected row
  if (length(s)!= 0) {
    tagList(
      actionButton("info", "",icon("info-circle"),style="color:rgb(57,156,8);border-color:rgb(255,255,255)"),

      # Information Dialog Box
      bsModal("ObjectInfo", LeafTable[s,c(2)], "info", size = "large", # Enables Pop up Screen

              img(src= paste0("Original/",LeafTable$LeafNames[s]),width='44%',height='44%'), img(src= paste0("LesionColor/",LeafTable$LeafNames2[s]),width='44%',height='44%')
      )
    )
  }
})

outputOptions(output, 'codeValidationInt', suspendWhenHidden = FALSE)
outputOptions(output, 'analysisFinish', suspendWhenHidden = FALSE)
outputOptions(output, 'warning', suspendWhenHidden = FALSE)
