library(plotly)
library(shinyBS)
cluster_analysis <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
    tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
  ),
  h2("Clustering Exposures"),
  fluidRow(
    box(
      width = 12,
      h3("Specify Model"),
      uiOutput(outputId = "select_res3"),
      uiOutput(outputId = "select_modality3"),
      uiOutput(outputId = "select_model3")
    )
  ),
  fluidRow(
    div(style = "position: relative;",box(
      width = 12,
      h3("Explore Number of Clusters"),
      selectInput(
        inputId = "metric",
        label = "Metric",
        choices = c("Within Cluster Sum of Squares" = "wss",
                    "Averaged Silhouette Coefficient" = 'silhouette',
                    "Gap Statisitc" = "gap_stat")
      ),
      selectInput(
        inputId = "algorithm1",
        label = "Algorithm",
        choices = c("k-means" = "kmeans", "Hierarchical" = "hclust", 
                    "Hierarchical k-means" = "hkmeans", "k-medoids" = "pam", "CLARA" = "clara"),
        selected = "hclust"
      ),
      add_busy_spinner(spin = "fading-circle"),
      uiOutput(outputId = "no_cluster1"),
      checkboxInput(inputId = "proportional2", label = "Proportional", value = TRUE),
      actionButton(inputId = "explore", label = "Explore"),
      tags$div(id = "insert_explore_plot"),
      tags$a(href = "#", 
             tags$i(class = "fas fa-question-circle"),
             title = "Need help?", 
             `data-toggle` = "popover", 
             `data-trigger` = "focus", 
             `data-content` = "The “Clustering” subtab provides several algorithms from the factoextra and 
             cluster packages to cluster samples based on exposure levels. After selecting the musica result 
             object, it is recommended to use the “Explore Number of Clusters” box to help find the potential 
             number of clusters in your data.",
             `data-html` = "true",
             `data-placement` = "left",
             style = "position: absolute; top: 5px; right: 5px; cursor: pointer;") 
    ))
  ),
  fluidRow(
    div(style = "position: relative;",
    box(
      width = 12,
      h3("Clustering"),
      uiOutput(outputId = "no_cluster2"),
      checkboxInput(inputId = "proportional3", label = "Proportional", value = TRUE),
      selectInput(
        inputId = "algorithm2",
        label = "Clustering Algorithm",
        choices = c("k-means" = "kmeans", "Hierarchical" = "hclust", 
                    "Hierarchical k-means" = "hkmeans", "k-medoids" = "pam", "CLARA" = "clara"),
        selected = "hclust"
      ),
      uiOutput(outputId = "diss"),
      tags$div(id = "hclust"),
      tags$div(id = "clara"),
      tags$div(id = "iter"),
      actionButton(inputId = "cluster_calc", label = "Clustering"),
      tags$div(id = "insert_cluster_table"),
      tags$a(href = "#", 
             tags$i(class = "fas fa-question-circle"),
             title = "Need help?", 
             `data-toggle` = "popover", 
             `data-trigger` = "focus", 
             `data-content` = "The “Clustering” box is where you perform the clustering analysis. 
             In addition to clustering algorithm, several methods for calculating dissimilarity matrix, 
             imported from the philentropy package, are also provided.",
             `data-html` = "true",
             `data-placement` = "left",
             style = "position: absolute; top: 5px; right: 5px; cursor: pointer;")       
    ))
  ),
  fluidRow(
    div(style = "position: relative;",
    box(
      width = 12,
      h3("Visualization"),
      radioButtons(
        inputId = "group2",
        label = "Group By",
        choices = list("None" = "none",
                       "Signature" = "signature",
                       "Annotation" = "annotation"),
        inline = TRUE,
        selected = "none"
      ),
      tags$div(id = "insert_annot2"),
      checkboxInput(inputId = "plotly3", label = "Plotly", value = TRUE),
      actionButton(inputId = "cluster_vis", label = "Visualize"),
      tags$div(id = "cluster_plot_div"),
      tags$a(href = "#", 
             tags$i(class = "fas fa-question-circle"),
             title = "Need help?", 
             `data-toggle` = "popover", 
             `data-trigger` = "focus", 
             `data-content` = "In the “Visualization” box, users can make scatter plots to visualize the 
             clustering results on a UMAP panel. Three types of plots are provided.If “None” is selected, 
             a single scatter plot will be made with points colored by cluster label.If “Signature” is selected, 
             a subplot is made for each combination of cluster label and signature. Points are colored by the 
             level of the specific signature. If “Annotation” is selected, an additional select box will show up 
             and let you choose one type of user-supplied annotation of interest. A subplot is generated for each
             combination of signature and category in the annotation. Points are colored by cluster label.",
             `data-html` = "true",
             `data-placement` = "left",
             style = "position: absolute; top: 5px; right: 5px; cursor: pointer;") 
      
    )
    )
  )
)
