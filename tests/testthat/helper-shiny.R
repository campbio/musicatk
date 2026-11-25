# Packages the Shiny app attaches with library(), in app.R, ui.R and
# server.R, apart from musicatk and TCGAbiolinks, plus ggplot2, which plotly
# attaches. Later ones mask earlier ones, as on the app's search path.
shiny_app_packages <- c("shiny", "shinyjs", "shinydashboard", "shinybusy",
                        "ggplot2", "plotly", "sortable", "shinyBS",
                        "shinyalert", "shinyjqui")

# Environment holding the app's server(), sourced from inst/shiny/server.R,
# for shiny::testServer(). The packages above are reached through a parent
# environment holding their exports, and library() is a no-op, so sourcing
# the app doesn't attach anything for the other tests. TCGAbiolinks is left
# out; only the TCGA import tab uses it. Mock functions before calling this,
# because their exports are copied when it runs.
shiny_server_env <- function() {
  for (pkg in shiny_app_packages) {
    testthat::skip_if_not_installed(pkg)
  }
  pkgs <- new.env(parent = as.environment("package:musicatk"))
  for (pkg in shiny_app_packages) {
    for (name in getNamespaceExports(pkg)) {
      value <- tryCatch(getExportedValue(pkg, name),
                        error = function(e) NULL)
      if (!is.null(value)) {
        assign(name, value, envir = pkgs)
      }
    }
  }
  pkgs$library <- function(...) invisible(NULL)
  env <- new.env(parent = pkgs)
  sys.source(system.file("shiny", "server.R", package = "musicatk"),
             envir = env)
  env
}

# Replaces shinyalert::shinyalert() for the calling test and returns a
# function that lists the alerts shown so far, each as list(title, text,
# type).
record_shinyalerts <- function(env = parent.frame()) {
  alerts <- list()
  testthat::local_mocked_bindings(
    shinyalert = function(title = "", text = "", type = "", ...) {
      alerts[[length(alerts) + 1]] <<- list(title = title, text = text,
                                            type = type)
      invisible(NULL)
    },
    .package = "shinyalert",
    .env = env
  )
  function() alerts
}

# Replaces shiny::showNotification() for the calling test and returns a
# function that gives the notification texts shown so far.
record_notifications <- function(env = parent.frame()) {
  notes <- character()
  testthat::local_mocked_bindings(
    showNotification = function(ui, ...) {
      notes <<- c(notes, as.character(ui))
      invisible(NULL)
    },
    .package = "shiny",
    .env = env
  )
  function() notes
}

# Replaces shiny::showModal() for the calling test and returns a function
# that gives the number of modals shown so far.
record_modals <- function(env = parent.frame()) {
  n <- 0
  testthat::local_mocked_bindings(
    showModal = function(...) {
      n <<- n + 1
      invisible(NULL)
    },
    .package = "shiny",
    .env = env
  )
  function() n
}
