library(plotly)
library(shinyBS)
cluster_analysis <- fluidPage(
  h2("Clustering Exposures"),
  fluidRow(
    box(
      width = 12,
      uiOutput(outputId = "select_res3")
    )
  ),
  fluidRow(
    box(
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
      tags$div(id = "insert_explore_plot")
    )
  ),
  fluidRow(
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
      tags$div(id = "insert_cluster_table")
    )
  ),
  fluidRow(
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
      tags$div(id = "cluster_plot_div")
    )
  )
)
