#UI.R
#loading shiny library
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinyFiles)

source("fonctions_apprentissage.r")

header <- dashboardHeader(title = "Lesion Count Tools")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Calibration", tabName = "tabCalibration", icon = icon("balance-scale")),
    menuItem("Analysis", tabName = "tabAnalysis", icon = icon("pagelines")),
    menuItem("Theme", tabName = "tabTheme", icon = icon("dashboard"))
  )
)

body <- dashboardBody(
  includeCSS('www/styles.css'),
  tabItems(
    # add tab for calibration
    source(file.path("ui_code", "tabCalibrationUI.R"), local = TRUE)$value,
      
    # add tab for analysis
    source(file.path("ui_code", "tabAnalysisUI.R"), local = TRUE)$value,

    # other tab
    tabItem(
      tabName = "tabTheme",
      h1("Theme"),
      shinythemes::themeSelector()  # <--- Add this somewhere in the UI
    )
  )
)






shinyUI(
  dashboardPage(skin = "green", header, sidebar, body)
)


