# dir

shinyDirChoose(
  input,
  'dir',
  filetypes = c('', 'txt', 'Rdata', "png", "csv", "*"),
  roots = allVolumesAvail,
  session = session,
  restrictions = system.file(package = 'base')
)

# dir
dir <- reactive(input$dir)
output$dir <- renderText(updateDir())

updateDir <- eventReactive(input$dir,{
  home <- normalizePath(allVolumesAvail[input$dir$root])
  datapath <<- file.path(home,paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
  exitStatus <<- list(code=-1, mess = "NULL", err = "NULL")
  return(datapath)
})

result <- eventReactive(input$runButton,{
  dirCalibration <- existDirCalibration(datapath)
  if(dirCalibration$dirLimbe == TRUE && dirCalibration$dirBackground == TRUE && dirCalibration$dirLesion == TRUE){
    withProgress(message = 'Making calibration, please wait\n', value = 0, {
      
      # Increment the progress bar, and update the detail text.
      incProgress(1/3, detail = "Doing calibration")
      exitStatus <<- apprentissage(datapath,"background","limbe","lesion")
      incProgress(3/3, detail = "End of calibration")
    })
    return(fileRData)
  }
  else{
    # print(paste("else inputdir",datapath))
    
    
    errorMess <-tags$div("Error not find all sub-directories !!!!:",  tags$br(),
                         tags$ul(
                           tags$li(paste("limbe: ", dirCalibration$dirLimbe)),
                           tags$li(paste("background: ", dirCalibration$dirBackground)),
                           tags$li(paste("lesion: ", dirCalibration$dirLesion))
                         )
    )
    exitStatus <<-list(code=0, err=errorMess)
  }
})

exitStatus <- observe(result())

output$code <- renderText({
  result()
  updateDir()
  exitStatus$code
  
})

output$mess <- renderText({
  result()
  fileRData
})
output$err <- renderPrint({
  result()
  exitStatus$err
})

output$img <- renderImage({
  result()
  # Return a list containing the filename
  list(src = plotFileCalibration,
       width = 400,
       height = 400,
       alt = "plot img")
  
}, deleteFile = FALSE)

output$table <- renderTable({
  result()
  outCalibrationTable
  
},striped = TRUE, bordered = TRUE,
align = 'c',
rownames = TRUE)

# observe({ print(output$code) })
# observe({ print(input$dir) })
outputOptions(output, 'code', suspendWhenHidden = FALSE)
outputOptions(output, 'mess', suspendWhenHidden = FALSE)
outputOptions(output, 'err', suspendWhenHidden = FALSE)
