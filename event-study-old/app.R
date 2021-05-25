library(shinythemes)
library(shinydashboard)
library(shiny)
library(shinyjs)

ui <- fluidPage(
  
  titlePanel("Econometrics sandbox"),
  
  sidebarLayout(
    
    sidebarPanel(
      p("Oops, this page has moved!"),
      p("Click",
        HTML('<a href="https://datanalytics.worldbank.org/connect/#/apps/674/access" target = "_blank">here</a>'),
        "to be redirected"
      )
    ),
    
    mainPanel(
    )
  )
)

server <- function(input, output) {}

shinyApp(ui, server)

