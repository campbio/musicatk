library(shinydashboard)
library(shinyjs)
library(shinyBS)

source("ui_import.R",local = T)
source("ui_import_musica.R",local = T)
#source("ui_genome.R", local = T)
source("ui_musica.R",local = T)
#source("ui_test.R",local = T)
source("ui_discover.R", local = T)
source("ui_tables.R", local = T)
source("ui_resultvisualization.R", local = T)
source("ui_predict.R", local = T)
source("ui_annotations.R", local = T)
source("ui_compare.R", local = T)
source("ui_help.R", local = T)
source("ui_heatmap.R",local = T)
source("ui_cluster.R",local = T)
source("ui_differentialanalysis.R", local = T)
ui <- fluidPage(
  shinyalert::useShinyalert(),
  useShinyjs(),
  dashboardPage(
    dashboardHeader(title = "musicatk"),
    dashboardSidebar(sidebarMenu(
      menuItem("Import", tabName = "import",
               menuSubItem("Import Files", "import"),
               menuSubItem("Import Musica Result Object", "musica_result"),
               menuSubItem("Import Annotations", "annotations")),
      #menuItem("Genome", tabName = "genome"),
      menuItem("Create Musica Object", tabName = "musica"),
      menuItem("Build Tables", tabName = "tables"),
      menuItem("Signatures and Exposures", tabName = "signatures",
               menuSubItem("Discover", "discover"),
               menuSubItem("Predict", "predict")),
      menuItem("Data Visualization", tabName = "visualization"),
      menuItem("Additional Analysis", tabName = "downstream",
               menuSubItem("Compare Signatures", tabName = "compare"),
               menuSubItem("Exposure Differential Analysis", 
                           tabName = "differentialanalysis"),
               menuSubItem("Clustering", tabName = "cluster"),
               menuSubItem("Heatmap", tabName = "heatmap")),
      #menuItem("Test", tabName = "test", icon = icon("th")),
      menuItem("Help", tabName = "widgets"))),
    
    dashboardBody(
        tabItems(
          tabItem(tabName = "import", h2("Import Data", shinyPanelImport)),
          tabItem(tabName = "musica_result",h2("Upload Musica"),shinyPanelResult),
          #tabItem(tabName = "genome",  shinyPanelGenome),
          tabItem(tabName = "musica", h2("Create Musica Object"),shinyPanelMusica),
          #tabItem(tabName = "test", shinyPaneltest),
        
###################### Nathan's Code ##########################################
          tabItem(tabName = "tables", h2("Create Tables"), shinyPanelTables),
          tabItem(tabName = "annotations", h2("Add Sample Annotations"), 
                  shinyPanelAnnotations),
          tabItem(tabName = "discover", h2("Discover Signatures and Exposures"),
                shinyPanelDiscover),
          tabItem(tabName = "predict", h2("Predict Signature Exposures"),
                  shinyPanelPredict),
          tabItem(tabName = "compare", h2("Compare Signatures"),
                  shinyPanelCompare),
          tabItem(tabName = "differentialanalysis", h2("Differential Analysis"),
                  shinyPanelDifferentialAnalysis),
###############################################################################
		      tabItem(tabName = "visualization",
                  musicaresultvisualization),
		 tabItem(tabName = "heatmap",h2("Plot heatmap"),
                  shinyPanelheatmap),
          tabItem(tabName = "cluster", cluster_analysis)
      )
    )
  )
)
