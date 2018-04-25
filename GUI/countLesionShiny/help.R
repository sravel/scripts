library(shiny)
library(ggvis)

source("fonctions_apprentissage.r")

path <- "/media/sebastien/Bayer/AnalyseImagesV4/Exemple2/Samples"

apprentissage(path,"background","limbe","lesion")


# runApp(list(
#   ui = bootstrapPage(
#     ggvisOutput("p"),
#     uiOutput("p_ui")
#   ),
#   server = function(..., session) {
#     
#     df4 %>%
#       ggvis(~LD2, ~LD1) %>%
#       layer_points(fill = ~group) %>%
#       bind_shiny("p", "p_ui")
# 
#   }
# ))