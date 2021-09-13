#' Starts the musicatk interactive Shiny app
#'
#' The musicatk Shiny app allows users
#' to perform mutational signature analysis using an interative graphical user
#' interface (GUI)
#'
#' @param include_version Include the version number in the header.
#' Default \code{TRUE}.
#' @param theme The theme to use for the GUI. Default \code{"yeti"}.
#'
#' @return The shiny app will open. No data will be returned.
#' @export
#' @examples
#' \dontrun{
#' # Start the app
#' musicatk()
#' }
musicatk <- function(include_version=TRUE, theme='yeti') {
  appDir <- system.file("shiny", package = "musicatk")
  shiny::shinyOptions(include_version = include_version)
  shiny::shinyOptions(theme = theme)
  shiny::runApp(appDir, display.mode = "normal")
}
