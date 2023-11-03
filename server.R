#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)


# tmp data test

file <- "https://github.com/rstudio-education/shiny-course/raw/main/movies.RData"
destfile <- "movies.RData"

download.file(file, destfile)

# Load data

load("movies.RData")


# Define server logic required to draw a histogram
function(input, output, session) {
  observeEvent(input$diff_zscore, {
    if(input$diff_zscore == "increase"){
      disable("neg_zscore")
      enable("pos_zscore")
    } else if(input$diff_zscore == "decrease"){
      disable("pos_zscore")
      enable("neg_zscore")
    } else {
      enable("pos_zscore")
      enable("neg_zscore")
    }
  })
  
  rv <- reactiveValues(count = 1)
  
  observeEvent(input$minus,{
    rv$count <- rv$count - 1
    if(rv$count == 1){
      disable("minus")
    } else {
      enable("plus")
    }
  })
  
  observeEvent(input$plus,{
    rv$count <- rv$count + 1
    if(rv$count == 3){
      disable("plus")
    } else {
      enable("minus")
    }
  })
  
  output$count <- renderText(rv$count)
  
  test_data <- eventReactive(input$GO_button, ignoreInit = TRUE, {
    plot1 <- ggplot(data = movies, aes_string(x = "thtr_rel_month", y = "thtr_rel_day")) +
      geom_point()
    plot2 <- ggplot(data = movies, aes_string(x = "thtr_rel_month", y = "thtr_rel_day")) +
      geom_violin()
    plot3 <- ggplot(data = movies, aes_string(x = "thtr_rel_month", y = "thtr_rel_day")) +
      geom_boxplot()
    plots <- list(plot1, plot2, plot3)
    return(plots)
  })
  
  observeEvent(input$GO_button, ignoreInit = TRUE, {
    
    disabled(actionButton("minus", "-1"))
    actionButton("plus", "+1")
    
    output$minus <- renderUI({
      disabled(actionButton("minus", "-1"))
    })
    
    output$plus <- renderUI({
      actionButton("plus", "+1")
    })
  })
  
  output$test <- renderPlot({
    if (rv$count < 1 || rv$count > 3){
      return()
    }
    test_data()[rv$count]
  })

}
