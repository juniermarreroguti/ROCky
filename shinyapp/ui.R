library(shiny)

# Define UI for miles per gallon application

navbarPage(

    "ROCky",
    tabPanel("Introdução" , uiOutput('page1')),
    tabPanel("Setting", uiOutput('page2')),
    tabPanel("Resultados" , uiOutput('page3'))

)
