#####################################################################################################
#
# Copyright 2018 CIRAD-INRA
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to CIRAD and South Green developpement plateform
# Version 0.1.0 written by Sebastien RAVEL, François BONNOT, Sajid ALI, FOURNIER Elisabeth
#####################################################################################################




############################################
## Calibration Directory path
############################################
# option to load directory path bottom for calibration folder
shinyDirChoose(
  input,
  'dirCalibration',
  filetypes = c('', 'txt', 'Rdata', "png", "csv", "*"),
  roots = allVolumesAvail,
  session = session,
  restrictions = system.file(package = 'base')
)

# return to UI path selected for calibration
output$dirCalibration <- renderText(updateDir())

# when click to bottom update path
updateDir <- eventReactive(input$dirCalibration,{
  home <- normalizePath(allVolumesAvail[input$dirCalibration$root])
  datapath <<- file.path(home,paste(unlist(input$dirCalibration$path[-1]), collapse = .Platform$file.sep))
  exitStatus <<- list(code=-1, mess = "NULL", err = "NULL")
  return(datapath)
})

############################################
## Run action
############################################

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


############################################
## Output when run click
############################################
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

outputOptions(output, 'code', suspendWhenHidden = FALSE)
outputOptions(output, 'mess', suspendWhenHidden = FALSE)
outputOptions(output, 'err', suspendWhenHidden = FALSE)
