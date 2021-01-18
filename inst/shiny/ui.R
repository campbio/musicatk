library(shinydashboard)

ui <- fluidPage(
  
  dashboardPage(
    dashboardHeader(title = "Musicatk"),
    dashboardSidebar(sidebarMenu(
      menuItem("Import", tabName = "import", icon = icon("th")),
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
            actionButton("MusicaResults","Get Musica Result Object"))
          
          
        ),
        tabItem(tabName = "signatures", h2("Signatures and Exposures"),
            box(
              selectInput("SelectTable", h3("Select Count Table"),
                          choices = list("SBS96" = 1, "SBS192" = 2, 
                                         "DBS" = 3, "Indel" = 4), 
                                         selected = 1),
              textInput("NumberOfSignatures", h3("Number of signatures")),
              
              radioButtons("Method", h3("Method"),
                           choices = list("LDA" = 1, "NMF" = 2)),
              textInput("Seed", h3("Seed")),
              textInput("nStart", h3("Number of random starts")),
              actionButton("MusicaResults", h3("Discover Signatures"))
            )
        )
        
      )
    )
  ))
