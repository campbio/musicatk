library(shinydashboard)

source("ui_discover.R", local = T)
source("ui_tables.R", local = T)

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
        tabItem(
          tabName = "import",
          h2("Import Data "),
          box(
            radioButtons("data", "Select Data:",
                         c("MAF file" = "maf",
                           "VCF file" = "vcf",
                           "Example Data" = "exp")),
            actionButton("get_musica", "Get Musica Data"),
            actionButton("MusicaResults","Get Musica Result Object"),
            fileInput("file_vcf", "VCF File",
                      multiple = TRUE,
                      accept = ".vcf"),
            fileInput("file_maf", "MAF File",
                      multiple = TRUE,
                      accept = ".maf"),
            mainPanel(
              
              # Output: Data file ----
              tableOutput("contents")
              
            )
          
          
        )),
        
###################### Nathan's Code ##########################################
        tabItem(tabName = "tables", h2("Create Tables"), shinyPanelTables),
        tabItem(tabName = "signatures", h2("Signatures and Exposures"), 
                shinyPanelDiscover)
      )
###############################################################################


    )
  ))
