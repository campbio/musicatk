library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(shinyalert)
library(shinybusy)

source("ui_import_tcga.R",local = T)
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
source("ui_download.R",local = T)

ui <- fluidPage(
  shinyalert::useShinyalert(),
  shinybusy::add_busy_spinner(),
  shinyjs::useShinyjs(),
  dashboardPage(
    dashboardHeader(title = "musicatk"),
    dashboardSidebar(sidebarMenu(
      id = "menu",
      tags$head(tags$style(".inactiveLink {
                           position: relative;
                           cursor: not-allowed;
                           }")),
      menuItem("Import", tabName = "import",
               menuSubItem("Import Files", "import"),
               menuSubItem("Import TCGA datasets", "import_tcga"),
               menuSubItem("Import Musica Result Object", "musica_result")),
      #menuItem("Genome", tabName = "genome"),
      menuItem("Create Musica Object", tabName = "musica"),
      menuItem("Import Annotations", tabName = "annotations"),
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
      menuItem("Download", tabName = "download"),
      menuItem("Help", tabName = "widgets"))),
    
    dashboardBody(
        tabItems(
          tabItem(tabName = "import_tcga", h2("Import TCGA Datasets", shinyPanelTCGA)),
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
          tabItem(tabName = "cluster", cluster_analysis),
          tabItem(tabName = "download",h2("Download musica result objects"), shinyPanelDownload)
      )
    )
  ),
# tags$script(charset="utf-8", HTML("var x = document.getElementById(\"sidebarItemExpanded\").querySelectorAll(\"li\");
#                var old_html = x[4].innerHTML;
#                var new_html = '<span class=\"musica_wrapper\">' + old_html + '</span>';
#                x[4].innerHTML = new_html;"))
)
