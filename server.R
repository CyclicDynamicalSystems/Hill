library(shiny)
library(ggplot2)
library(Cairo)
library(tidyr)
library(dplyr)

source("hill-core.R")

shinyServer(function(input, output) {
  output$overview <- renderText({ 
    h <- hill(n = input$n, a = input$a, b = input$b, g = input$g)
    h$overview()
  })
  output$plot <- renderPlot({
    h <- hill(n = input$n, a = input$a, b = input$b, g = input$g)
    sim <- h$simulate.multi(input$size, times = seq(0, input$time, by = 0.1), drop.start = input$drop.start)
    sim$plot.2d(input$proj1, input$proj2)
  }, width = 500, height = 500)
  output$eigenvalues <- renderPlot({
    h <- hill(n = input$n, a = input$a, b = input$b, g = input$g)
    h$eigen$plot()
  }, width = 500, height = 500)
}) 