library(shinydashboard)
library(shinyjs)

source("ui_import.R",local = T)
source("ui_genome.R", local = T)
source("ui_musica.R",local = T)
source("ui_discover.R", local = T)
source("ui_tables.R", local = T)
source("ui_03_resultvisualization.R", local = T)
source("ui_predict.R", local = T)
source("ui_annotations.R", local = T)
source("ui_compare.R", local = T)

ui <- fluidPage(
  shinyalert::useShinyalert(),
  useShinyjs(),
  dashboardPage(
    dashboardHeader(title = "Musicatk"),
    dashboardSidebar(sidebarMenu(
      menuItem("Import", tabName = "import", icon = icon("th")),
      menuItem("Genome", tabName = "genome", icon = icon("th")),
      menuItem("Musica", tabName = "musica", icon = icon("th")),
      menuItem("Tables", tabName = "tables", icon = icon("th")),
      menuItem("Annotations", tabName = "annotations", icon = icon("th")),
      menuItem("Signatures and Exposures", tabName = "signatures", icon = icon("th"),
               menuSubItem("Discover Signatures and Exposures", "discover"),
               menuSubItem("Predict Signature Exposures", "predict"),
               menuSubItem("Compare Signatures", "compare")),
      menuItem("Data Visualization", tabName = "visualization", icon = icon("th")),
      menuItem("Help", tabName = "widgets", icon = icon("th")))),
    
    dashboardBody(
        tabItems(
          tabItem(tabName = "import", h2("Import Data", shinyPanelImport)),
          tabItem(tabName = "genome", shinyPanelGenome),
          tabItem(tabName = "musica", shinyPanelMusica),
        
###################### Nathan's Code ##########################################
          tabItem(tabName = "tables", h2("Create Tables"), shinyPanelTables),
          tabItem(tabName = "annotations", h2("Add Sample Annotations"), 
                  shinyPanelAnnotations),
          tabItem(tabName = "discover", h2("Discover Signatures and Exposures"), 
                shinyPanelDiscover),
          tabItem(tabName = "predict", h2("Predict Known Signatures"),
                  shinyPanelPredict),
          tabItem(tabName = "compare", h2("Compare Signatures"), 
                  shinyPanelCompare),
###############################################################################
		  tabItem(tabName = "visualization",
                  musicaresultvisualization)
      )
    )
  ))
