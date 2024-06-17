###### DIMSpect GUI server part ######
# Graphical User Interface for DIMS database v0.2

# load functions
source("utils.R")

options(shiny.maxRequestSize = 300*1024^2)
library(ggplot2) 
library(DT)
library(shinydashboard)
library(shinyFiles)
library(ssh)
library(dplyr)
library(httr2)
library(shinyjs)
library(tidyr)
library(stringr)

server <- function(session, input, output) {
  api <- reactiveVal(request("http://127.0.0.1:8000"))

  selected_metabolites_violinplot <- reactive({
    if (is.null(input$table_violin_rows_selected)) {
      1
    } else {
      input$table_violin_rows_selected
    }
  })
  
  selected_metabolites_barplot <- reactive({
    if (is.null(input$patient_data_table_bar_rows_selected)) {
      1
    } else {
      input$patient_data_table_bar_rows_selected
    }
  })
  
  observeEvent(input$view_patient_selection, {
    if (input$view_patient_selection > 0) {
      enable("zscore_high")
      enable("zscore_low")
      enable("identified_only")
    } else {
      disable("zscore_high")
      disable("zscore_low")
      disable("identified_only")
    }
  })
  
  patient_data <- eventReactive(input$run_patient_query, {
    get_patient_data(api(), patient_query_input(), input$zscore_high, input$zscore_low, input$identified_only)
  })
  
  change_patient_data <- reactive({
    list(
      input$view_patient_selection,
      input$update_view_pat_sel
    )
  })
  
  patient_query_input <- eventReactive(input$view_patient_selection, {
    if(is.null(input$view_patient_selection)){
      return(NULL)
    } else {
      get_info_on_patient(api(), input$patient_id, input$sample_id, input$run_name, input$matrix_selected)
    }
  })
  
  
  #### view occurrences of patient in database ####
  output$patient_list <- DT::renderDT({
      validate(
        need(input$patient_id != "", "Give a patient ID")
      )
      patient_query_input()
      })
  
  observeEvent(input$run_patient_query, {
    updateTabsetPanel(session, "patient_query", selected = "Barplot")
  })
  
  output$patient_data_table_bar <- DT::renderDT({
      df <- patient_data()
      datatable(df, options = list(
        deferRender = TRUE,
        scrollY = TRUE,
        scrollX = TRUE,
        scroller = TRUE,
        autoWidth = FALSE,
        columnDefs = list(list(width = '10%', targets = c(3,12,13)))
        )
      )
    })

  output$table_violin <- DT::renderDT({
      df <- patient_data()
      datatable(df, options = list(
        deferRender = TRUE,
        scrollY = TRUE,
        scrollX = TRUE,
        scroller = TRUE,
        autoWidth = FALSE,
        columnDefs = list(list(width = '10%', targets = c(3,12,13)))
        )
      )
    })

  data_bar_plot <- reactive({
    df <- patient_data()
    metab_rows <- selected_metabolites_barplot()
    plot_data_selected <- select_plot_data(df, metab_rows)
    plot_data_long <- reshape2::melt(plot_data_selected, id.vars = "HMDB_name")
    plot_data_long$color = "green"
    plot_data_long$color[grep("P", plot_data_long$variable)] <- "blue"
    plot_data_long$color[grep(input$sample_id, plot_data_long$variable)] <- "purple"
    plot_data_long
  })
  
  data_violin_plot <- reactive({
    df <- patient_data()
    metab_rows <- selected_metabolites_violinplot()
    data_violin <- df[metab_rows, ]
    data_violin
  })
  
  output$bar_plot <- renderPlot({
    label_split <- function(label) (str_replace_all(label, paste0("(.{15})"), "\\1\n"))
    
    ggplot(data_bar_plot(), aes(variable, value, color)) +
      ggtitle("Metabolite") +
      geom_bar(aes(fill=color), stat="identity") + 
      labs(x="", y="Z-score") +
      facet_grid(HMDB_name ~ ., scales="free", switch = "y", labeller = as_labeller(label_split)) + 
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
            strip.background = element_rect(colour = "black", fill = "white", linewidth = 1, linetype = "solid"),
            strip.text = element_text(size = 12))
  })

  output$violin_plot <- renderPlot({
    all_samples <- unique(patient_query_input()$ID)
    violin_plot <- create_violin_plots(data_violin_plot(), all_samples)
    violin_plot
  })

}