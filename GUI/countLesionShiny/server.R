#SERVER.R
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinyFiles)


source("fonctions_apprentissage.r")
source("fonctions_analyse.r")

getOwnVolume <- function (exclude=NULL) 
{
  osSystem <- Sys.info()["sysname"]
  if (osSystem == "Darwin") {
    volumes <- list.files("/Volumes/", full.names = T)
    names(volumes) <- basename(volumes)
  }
  else if (osSystem == "Linux") {
    volumes <- c(root = "/")
    home <- c(home = "~")
    media <- list.files("/media", full.names = T)
    names(media) <- media
    volumes <- c(home, media, volumes)
  }
  else if (osSystem == "Windows") {
    volumes <- system("wmic logicaldisk get Caption", intern = T)
    volumes <- sub(" *\\r$", "", volumes)
    keep <- !tolower(volumes) %in% c("caption", "")
    volumes <- volumes[keep]
    volNames <- system("wmic logicaldisk get VolumeName", 
                       intern = T)
    volNames <- sub(" *\\r$", "", volNames)
    volNames <- volNames[keep]
    volNames <- paste0(volNames, ifelse(volNames == "", 
                                        "", " "))
    volNames <- paste0(volNames, "(", volumes, ")")
    names(volumes) <- volNames
  }
  else {
    stop("unsupported OS")
  }
  if (!is.null(exclude)) {
    volumes <- volumes[!names(volumes) %in% exclude]
  }
  volumes
}

allVolumesAvail = getOwnVolume()


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
  output$debug <- renderPrint({
    sessionInfo()
  })
  session$onSessionEnded(stopApp)
})
