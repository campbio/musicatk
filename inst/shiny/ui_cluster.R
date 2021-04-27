library(plotly)
library(shinyBS)
cluster_analysis <- fluidPage(
  h2("Clustering Exposures"),
  uiOutput(outputId = "select_res3"),
  fluidRow(
    box(
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
                    "Hierarchical k-means" = "hkmeans", "k-medoids" = "pam", "CLARA" = "clara")
      ),
      uiOutput(outputId = "no_cluster1"),
      checkboxInput(inputId = "proportional2", label = "Proportional", value = TRUE),
      actionButton(inputId = "explore", label = "Explore"),
      plotOutput(outputId = "explore_plot")
    )
  ),
  fluidRow(
    box(
      h3("Clustering"),
      uiOutput(outputId = "no_cluster2"),
      checkboxInput(inputId = "proportional3", label = "Proportional", value = TRUE),
      selectInput(
        inputId = "algorithm2",
        label = "Clustering Algorithm",
        choices = c("k-means" = "kmeans", "Hierarchical" = "hclust", 
                    "Hierarchical k-means" = "hkmeans", "k-medoids" = "pam", "CLARA" = "clara")
      ),
      uiOutput(outputId = "diss"),
      tags$div(id = "hclust"),
      tags$div(id = "clara"),
      tags$div(id = "iter"),
      actionButton(inputId = "cluster_calc", label = "Clustering"),
      tags$div(id = "insert_cluster_table")
    )
  ),
  fluidRow(
    box(
      h3("Visualization"),
      radioButtons(
        inputId = "group2",
        label = "Group By",
        choices = list("Signature" = "signature",
                       "Annotation" = "annotation",
                       "None" = "none"),
        inline = TRUE,
        selected = "signature"
      ),
      tags$div(id = "insertannot2"),
      actionButton(inputId = "cluster_vis", label = "Visualize"),
      plotOutput("cluster_plot")
    )
  )
)