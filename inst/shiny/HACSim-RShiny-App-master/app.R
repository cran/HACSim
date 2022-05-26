#
# This is the HACSim R Shiny Web Application. You can run the application by clicking
# the 'Run App' button above.
#

# # install.packages("shiny")
library(shiny)
# install.packages("HACSim")
library(HACSim)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("shinymeta")
library(shinymeta)
# install.packages("shinydashboard")
library(shinydashboard)
# install.packages("shinydashboardPlus")
library(shinydashboardPlus)
# install.packages("shinyWidgets")
library(shinyWidgets)
# install.packages("shinycssloaders")
library(shinycssloaders)
# install.packages("shinyjs")
library(shinyjs)
# install.packages("stringr")
library(stringr)

# load ui elements
source("ui.R")
# load server function
source("server.R")

# Run the application 
shinyApp(ui = ui, server = server)
