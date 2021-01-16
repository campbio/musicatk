ui <- fluidPage(
  
  dashboardPage(
    dashboardHeader(title = "Musicatk"),
    dashboardSidebar(sidebarMenu(
      menuItem("Import", tabName = "import", icon = icon("th")),
      menuItem("Quality Control", tabName = "widgets", icon = icon("th")),
      menuItem("Exploratory data analysis", tabName = "widgets", icon = icon("th")),
      menuItem("Data visualization", tabName = "widgets", icon = icon("th")),
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
            actionButton("get_musica", "Get Musica Data"))
          
          
        )
        
      )
    )
  ))
