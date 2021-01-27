library(shinydashboard)

source("ui_discover.R", local = T)
source("ui_tables.R", local = T)
source("ui_03_resultvisualization.R", local = TRUE)

ui <- fluidPage(
  shinyalert::useShinyalert(),
  
  dashboardPage(
    dashboardHeader(title = "Musicatk"),
    dashboardSidebar(sidebarMenu(
      menuItem("Import", tabName = "import", icon = icon("th")),
      menuItem("Tables", tabName = "tables", icon = icon("th")),
      menuItem("Signatures", tabName = "signatures", icon = icon("th")),
      menuItem("Data Visualization", tabName = "visualization", icon = icon("th")),
      menuItem("Help", tabName = "widgets", icon = icon("th")))),
    
    dashboardBody(
      tabItems(
        tabItem(tabName = "import"),
        
###################### Nathan's Code ##########################################
        tabItem(tabName = "tables", h2("Create Tables"), shinyPanelTables),
        tabItem(tabName = "signatures", h2("Signatures and Exposures"), 
                shinyPanelDiscover),
###############################################################################
        
        tabItem(tabName = "visualization",
                musicaresultvisualization)
      )
    )
  ))
