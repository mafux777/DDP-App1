library(shiny)
library(leaflet)

shinyUI(pageWithSidebar(
        headerPanel("Great Circle Distance"),
        sidebarPanel(
                numericInput(inputId='lat1', label='Latitude (rad)', value = .898451866, step=0.01),
                numericInput(inputId='lon1', label='Longitude (rad)', value = .008052755, step=0.01), # London
                numericInput(inputId='lat2', label='Latitude (rad)', value = .710218106, step=0.01),
                numericInput(inputId='lon2', label='Longitude (rad)', value = -1.294486466, step=0.01), # New York
                tags$a(href = "https://mafux777.github.io/DDP-App1/index.html", "Slidified Doc.")
        ),
        mainPanel(
                h3('Results'),
                leafletOutput('gcMap'),
                textOutput('title')
        )
))
