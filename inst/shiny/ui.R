library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(shinyalert)
library(shinybusy)
library(TCGAbiolinks)
source("ui_import_tcga.R", local = TRUE)
source("ui_import.R", local = TRUE)
source("ui_import_musica.R", local = TRUE)
source("ui_musica.R", local = TRUE)
source("ui_discover.R", local = TRUE)
source("ui_tables.R", local = TRUE)
source("ui_resultvisualization.R", local = TRUE)
source("ui_predict.R", local = TRUE)
source("ui_annotations.R", local = TRUE)
source("ui_compare.R", local = TRUE)
source("ui_help.R", local = TRUE)
source("ui_heatmap.R", local = TRUE)
source("ui_cluster.R", local = TRUE)
source("ui_differentialanalysis.R", local = TRUE)
source("ui_download.R", local = TRUE)

ui <- fluidPage(
  #shinyalert::useShinyalert(),
  shinybusy::add_busy_spinner(),
  shinyjs::useShinyjs(),
  tags$head(
    tags$script(HTML("
      $(document).on('shiny:sessioninitialized', function(event) {
        // Code to expand all menu items
        $('.treeview-menu').css('display', 'block');
      });
    "))
  ),
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
               menuSubItem("Import Musica", "musica_result")),
      #menuItem("Genome", tabName = "genome"),
      menuItem("Create Musica Object", tabName = "musica"),
      menuItem("Annotations", tabName = "annotations"),
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
      menuItem(tags$a(href="https://camplab.net/musicatk/current/articles/articles/tutorial_tcga_ui.html", "Help")))),
    dashboardBody(
        tabItems(
          tabItem(tabName = "import_tcga",
                  h2("Import TCGA Datasets"), shiny_panel_tcga),
          tabItem(tabName = "import",
                  h2("Import Data"), shiny_panel_import),
          tabItem(tabName = "musica_result", h2("Upload Musica"),
                  shiny_panel_result),
          tabItem(tabName = "musica", h2("Create Musica Object"),
                  shiny_panel_musica),

###################### Nathan's Code ##########################################
          tabItem(tabName = "tables", h2("Create Tables"), shiny_panel_tables),
          tabItem(tabName = "annotations", h2("Add Sample Annotations"),
                  shiny_panel_annotations),
          tabItem(tabName = "discover", h2("Discover Signatures and Exposures"),
                shiny_panel_discover),
          tabItem(tabName = "predict", h2("Predict Signature Exposures"),
                  shiny_panel_predict),
          tabItem(tabName = "compare", h2("Compare Signatures"),
                  shiny_panel_compare),
          tabItem(tabName = "differentialanalysis", h2("Differential Analysis"),
                  shiny_panel_diffanal),
###############################################################################
          tabItem(tabName = "visualization",
                  musicaresultvisualization),
          tabItem(tabName = "heatmap", h2("Plot heatmap"),
                  shiny_panel_heatmap),
          tabItem(tabName = "cluster", cluster_analysis),
          tabItem(tabName = "download", h2("Download musica objects"),
                  shiny_panel_download)
      )
    )
  )
)
