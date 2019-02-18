
library("shiny")
library("ggplot2")
library("scales")
library("shinythemes")
library("visNetwork")
###################################################
###############network summary#####################
#' @export
pagnet.mra.interface <- function() {
  appDir <- system.file("shiny", package = "PAGnet")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
