library(shiny)
library(ggplot2)

shinyUI(pageWithSidebar(
  
  headerPanel("Hill analysis"),
  
  sidebarPanel(
    withMathJax(),
    sliderInput('n', 'n', min = 3, max = 35, value = 3, step = 2),
    sliderInput('a', 'Alpha', min = 1, max = 10000, value = 18, step = 0.01),
    sliderInput('b', 'Beta', min = 0.05, max = 100, value = 1, step = 0.01),
    sliderInput('g', 'Gamma', min = 2, max = 100, value = 3, step = 0.1),
    sliderInput('time', 'Time', min = 10, max = 2000, value = 100, step = 1),
    sliderInput('size', 'Trajectories', min = 1, max = 20, value = 3, step = 1),
    checkboxInput('drop.start', 'Drop trajectory start', FALSE),
    textInput('proj1', 'proj1', value = "e2"),
    textInput('proj2', 'proj2', value = "e3"),
    textInput('proj3', 'proj3', value = "e1")
  ),
  
  mainPanel(
    verbatimTextOutput("overview"),
    plotOutput('eigenvalues', width = 500, height = 500),
    plotOutput('plot', width = 500, height = 500)
  )
)) 