library(shinydashboard)
library(shinyjs)

source("ui_discover.R", local = T)
source("ui_tables.R", local = T)
source("ui_predict.R", local = T)

ui <- fluidPage(
  shinyalert::useShinyalert(),
  useShinyjs(),
  # extendShinyjs(text = jscode),
  dashboardPage(
    dashboardHeader(title = "Musicatk"),
    dashboardSidebar(sidebarMenu(
      menuItem("Import", tabName = "import", icon = icon("th")),
      menuItem("Tables", tabName = "tables", icon = icon("th")),
      menuItem("Annotations", tabName = "annotations", icon = icon("th")),
      menuItem("Signatures", tabName = "signatures", icon = icon("th"),
               menuSubItem("Discover Signatures", "discover"),
               menuSubItem("Predict Signatures", "predict")),
      menuItem("Data Visualization", tabName = "visualization", icon = icon("th")),
      menuItem("Help", tabName = "widgets", icon = icon("th")))),
    
    dashboardBody(
        tabItems(
          tabItem(tabName = "import"),
        
###################### Nathan's Code ##########################################
          tabItem(tabName = "tables", h2("Create Tables"), shinyPanelTables),
          tabItem(tabName = "discover", h2("Discover Signatures and Exposures"), 
                shinyPanelDiscover),
          tabItem(tabName = "predict", h2("Predict Known Signatures"),
                  shinyPanelPredict)
###############################################################################

      )
    )
  ))
