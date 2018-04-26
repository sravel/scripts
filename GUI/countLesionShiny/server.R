#SERVER.R
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinyFiles)


source("fonctions_apprentissage.r")
source("fonctions_analyse.r")

existDirCalibration <- function(datapath){
  list(
    dirLimbe = file.exists(paste(datapath,"/limbe", sep = .Platform$file.sep)),
    dirBackground = file.exists(paste(datapath,"/background", sep = .Platform$file.sep)),
    dirLesion = file.exists(paste(datapath,"/lesion", sep = .Platform$file.sep))
  )
}

#writing server function
shinyServer(function(input, output, session) {

  # Load functions for tab calibration
  source(file.path("server_code", "tabCalibrationServer.R"), local = TRUE)$value
  
  # Load functions for tab analysis
  source(file.path("server_code", "tabAnalysisServer.R"), local = TRUE)$value
  
  # session$onSessionEnded(stopApp)
})
