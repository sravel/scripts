if (exists("fileRData")){
  lda1 <<- load(file = fileRData, envir = .GlobalEnv)
}else{
  lda1 <<- NULL
}


shinyDirChoose(
  input,'dirIn',
  filetypes = c('', "png", "PNG","jpg","JPG","jpeg","JPEG"),
  roots = allVolumesAvail,
  session = session,
  restrictions = system.file(package = 'base')
)
shinyDirChoose(
  input,'dirOut',
  filetypes = c('', "png", "PNG","jpg","JPG","jpeg","JPEG"),
  roots = allVolumesAvail,
  session = session,
  restrictions = system.file(package = 'base')
)

dirIn <- reactive(input$dirIn)

output$dirInAnalysis <- renderText({
  updateDirAnalysis()
})

updateDirAnalysis <- eventReactive(input$dirIn,{
  home <- normalizePath(allVolumesAvail[input$dirIn$root])
  datapathAnalysis <<- file.path(home, paste(unlist(dirIn()$path[-1]), collapse = .Platform$file.sep))
  exitStatusAnalysis <<- list(codeValidationInt=-1, codeAnalysis=-1, messAnalysis = "NULL", errAnalysis = "NULL")
  return(datapathAnalysis)
})

dirOut <- reactive(input$dirOut)

output$dirOutAnalysis <- renderText({
  updateDirOutAnalysis()
})

updateDirOutAnalysis <- eventReactive(input$dirOut,{
  home <- normalizePath(allVolumesAvail[input$dirOut$root])
  datapathOutAnalysis <<- file.path(home, paste(unlist(dirOut()$path[-1]), collapse = .Platform$file.sep))
  exitStatusAnalysis <<- list(codeValidationInt=-1, codeAnalysis=-1, messAnalysis = "NULL", errAnalysis = "NULL")
  return(datapathOutAnalysis)
})

#### For RData file

shinyFileChoose(input, 'files',
                roots=allVolumesAvail,
                filetypes=c('', 'rdata' , 'RData'))

output$fileRdata <- renderText({
  RDataIn()
})


RDataIn <- eventReactive(input$files,{
  infile <<-  normalizePath(as.character(parseFilePaths(roots=allVolumesAvail, input$files)$datapath))
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
  # resultAnalysis()
  datapathAnalysis <- "~/Bayer/AnalyseImagesV4/Exemple1/Images/"
  datapathOutAnalysis <- "~/Bayer/AnalyseImagesV4/Exemple1/"
  addResourcePath("Original",datapathAnalysis) # Images are located outside shiny App
  addResourcePath("LesionColor",datapathOutAnalysis) # Images are located outside shiny App

  addResourcePath("Original",datapathAnalysis) # Images are located outside shiny App
  addResourcePath("LesionColor",datapathOutAnalysis) # Images are located outside shiny App

  LeafNames <- list.files(datapathAnalysis, full.names=FALSE)
  LeafNames2 <- list.files(datapathOutAnalysis, full.names=FALSE, pattern = "*.png")
  LeafTable <<- data.frame(LeafNames,LeafNames2)
  LeafTable <<- within(LeafTable, thumbnail <- paste0("<img src='","Original/",LeafTable$LeafNames,"' height='60'></img>"))
  LeafTable <<- within(LeafTable, thumbnail2 <- paste0("<img src='","LesionColor/",LeafTable$LeafNames2 ,"' height='60'></img>"))

  responseDataFilter2 <- LeafTable[,c(1,3,4)]

  displayableData<-DT::datatable(data = as.data.frame(responseDataFilter2, stringAsFactors = FALSE, row.names = NULL),

                                 escape=FALSE,selection="single",rownames=FALSE,colnames=c("FileName","Original","LesionColor"),

                                 callback = JS("table.on('dblclick.dt', 'td', function() {
                                                 var row=table.cell(this).index().row;
                                                 Shiny.onInputChange('rows_home',[row, Math.random()])});
                                                 table.on('click.dt', 'td', function() {
                                                 var k=table.cell(this).index().row;
                                                 if(table.rows('.selected').indexes().toArray()!= '' && table.rows('.selected').indexes().toArray() == k){k=-1;}
                                                 Shiny.onInputChange('rows_up_home',[k, Math.random()]);
                                                 Shiny.onInputChange('row_path', table.rows(this).data().toArray());
                                                 });"),
                                 style = "bootstrap",
                                 
                                 options = list(

                                   paging=TRUE,searching = TRUE,ordering=TRUE,scrollY = 750,scrollCollapse=TRUE,server = FALSE

                                 ))

})

observeEvent(input$table2_rows_selected,
{
  showModal(      # Information Dialog Box
    modalDialog(
      title = paste(LeafTable[input$table2_rows_selected,c(2)], "info", sep = " "),
      size = "l",
      easyClose = TRUE,
      img(src= paste0("Original/",LeafTable$LeafNames[input$table2_rows_selected]),width='44%',height='44%'), img(src= paste0("LesionColor/",LeafTable$LeafNames2[input$table2_rows_selected]),width='44%',height='44%')
    )
  )
})

outputOptions(output, 'codeValidationInt', suspendWhenHidden = FALSE)
outputOptions(output, 'analysisFinish', suspendWhenHidden = FALSE)
outputOptions(output, 'warning', suspendWhenHidden = FALSE)
