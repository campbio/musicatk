library(shinydashboard)
source("ui_import.R",local = T)
#source("ui_genome.R", local = T)
source("ui_musica.R",local = T)
source("ui_discover.R", local = T)
source("ui_tables.R", local = T)
source("ui_test.R",local = T)

ui <- fluidPage(
  shinyalert::useShinyalert(),
  
  dashboardPage(
    dashboardHeader(title = "musicatk"),
    dashboardSidebar(sidebarMenu(
      menuItem("Import", tabName = "import", icon = icon("th")),
      #menuItem("Genome", tabName = "genome", icon = icon("th")),
      menuItem("Musica", tabName = "musica", icon = icon("th")),
      menuItem("Tables", tabName = "tables", icon = icon("th")),
      menuItem("Signatures", tabName = "signatures", icon = icon("th")),
      menuItem("Data Visualization", tabName = "visualization", icon = icon("th")),
      menuItem("Help", tabName = "widgets", icon = icon("th")))),
    
    dashboardBody(
        tabItems(
          tabItem(tabName = "import", h2("Import Data", shinyPanelImport)),
          #tabItem(tabName = "genome",  shinyPanelGenome),
          tabItem(tabName = "musica", shinyPanelMusica),
        
###################### Nathan's Code ##########################################
          tabItem(tabName = "tables", h2("Create Tables"), shinyPanelTables),
          tabItem(tabName = "signatures", h2("Signatures and Exposures"), 
                shinyPanelDiscover)
      )
###############################################################################


    )
  ))
