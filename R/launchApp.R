launchApp <- function() {
  appDir <- system.file("shiny", "HACSim-RShiny-App-master", package = "HACSim")
  runApp(appDir, display.mode = "normal")
}