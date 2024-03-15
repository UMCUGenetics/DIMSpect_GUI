###### DIMSpect GUI server part ######
# Graphical User Interface for DIMS database v0.2

# load functions
source("/Users/mraves2/Development/DIMSpect_GUI/utils.R")

# temporary: set directory containing test data
setwd("/Users/mraves2/Development/DIMSpect_GUI")
input_dir <<- "/Users/mraves2/Metabolomics/DIMSdb/Test_data_GUI" 
output_dir <<- "/Users/mraves2/Metabolomics/DIMSdb/Test_data_GUI" 
load("~/Metabolomics/DIMSdb/Test_data_GUI/test_database.RData")
sample_info <- as.data.frame(sample_info)
# load("~/Development/DIMSpect_GUI/test_input.RData")

options(shiny.maxRequestSize = 300*1024^2)
# options(bitmapType='cairo')
library(ggplot2) 
library(DT)
library("shiny")
library("shinydashboard")
library("shinyFiles")
library("ssh")
library("dplyr")

# write sessionInfo to file
writeLines(capture.output(sessionInfo()), paste(output_dir, "sessionInfo.txt", sep="/"))

server <- function(input, output) {
  #### collect input ####
  # NB: alternative patient_id <- reactive({ input$patient_id })
  # doesn't work: Warning: Error in !=: comparison (2) is possible only for atomic and list types
  observe({
      patient_id <<- input$patient_id
      print(patient_id)
  })
  
  observe({
    sample_id <<- input$sample_id
    print(sample_id)
  })
  
  observe({
    run_name <<- input$run_name
    print(run_name)
  })
  
  observe({
    matrix_selected <<- input$matrix_selected
    print(matrix_selected)
  })
  
  observe({
    zscore_high <<- input$zscore_high
    print(zscore_high)
  })
  
  observe({
    zscore_low <<- input$zscore_low
    print(zscore_low)
  })
  
  observe({
    identified_only <<- input$identified_only
    print(identified_only)
  })
  
  observe({
    metabolite_id <<- input$metabolite_id
    print(metabolite_id)
  })
  
  observe({
    input_mass <<- input$input_mass
    print(input_mass)
  })

  observe({
    input_scan_mode <<- input$input_scan_mode
    print(input_scan_mode)
  })
  
  observe({
    selected_metabolites <<- input$patient_data_table_rows_selected
    # take the first metabolite as default
    if (is.null(selected_metabolites)) {
      selected_metabolites <<- 1
    }
    print(selected_metabolites)
  })

  #### view occurrences of patient in database ####
  observeEvent(input$view_patient_selection, {
    if(is.null(input$view_patient_selection)){
      return(NULL)
    } else {
      patient_query_input <<- get_info_on_patient(patient_id, sample_id, run_name, matrix_selected)
      output$patient_list <- shiny::renderTable({
        patient_query_input
      })
    }
    
    # use checkboxes to make a selection from the table
    # output$select_for_query <- renderUI({
    #   checkboxGroupInput(inputId = "run_names", label = "Keep", choices = levels(as.factor(patient_info$run_name)), selected = FALSE)
    # })
    
  }) # end observeEvent view_patient_selection
  
  # display table
  # output$patient_data_table <- DT::renderDataTable({
  #   patient_data
  # })
  
  observeEvent(input$run_patient_query, {
    patient_data <<- get_patient_data(patient_query_input, zscore_high, zscore_low, identified_only)
    
    # display table
    output$patient_data_table <- DT::renderDT({
      datatable(patient_data, options = list(
        deferRender = TRUE,
        scrollY = TRUE,
        scrollX = TRUE,
        scroller = TRUE,
        autoWidth = FALSE,
        columnDefs = list(list(width = '10%', targets = c(3,12,13)))
      ))  # %>% formatStyle(columns = c(2,3), width='20px')
    })
    
    # this gives error Error in [: invalid subscript type 'closure' in line 
    # plot_data_selected <- plot_data[selected_metabolites, ]
    # selected_metabolites <- reactive({ input$patient_data_table_rows_selected })
    # if (is.null(selected_metabolites)) { 
    #   selected_metabolites <<- 1
    # }
    # print(selected_metabolites)
    
    # display bar plots
    # TODO: make barplots reactive to selected_metabolites
    output$bar_plot <- renderPlot({
      if (!is.null(patient_data)) {
        info_columns <- grep("HMDB_name|m_z", colnames(patient_data))
        zscore_columns <- grep("_Zscore", colnames(patient_data))
        plot_names <- gsub("_Zscore", "", colnames(patient_data)[zscore_columns])
        plot_colors <- rep("green", length(plot_names))
        plot_colors[grep("P", plot_names)] <- "blue"
        plot_colors[grep(sample_id, plot_names)] <- "purple"
        # selection of metabolites from table (https://yihui.shinyapps.io/DT-rows)
        plot_data <- patient_data[1, c(info_columns[2], zscore_columns)]
        # reduce plot_data to only for selected metabolites
        plot_data_selected <- plot_data[selected_metabolites, ]
        # put data in long format for ggplot
        # empty_rows <- which(is.na(patient_data$HMDB_code))
        # plot_data_mincols <- patient_data %>% select(-c("HMDB_code", "m_z" , "scanmode", "full_HMDB_code", "full_HMDB_name"))
        # plot_data <- plot_data_mincols[-empty_rows, ]
        plot_data_long <- reshape2::melt(plot_data_selected, id.vars = "HMDB_name")
        plot_data_long$color = "green"
        plot_data_long$color[grep("P", plot_data_long$variable)] <- "blue"
        plot_data_long$color[grep(sample_id, plot_data_long$variable)] <- "purple"
        ggplot(plot_data_long, aes(variable, value, color)) +
            ggtitle("Metabolite") +
            geom_bar(aes(fill=color), stat="identity") + 
            labs(x="", y="Z-score") +
            facet_grid(HMDB_name ~ ., scales="free") + 
            theme(legend.position = "none")
      }
    })
    
  }) # end observeEvent run_patient_query
  
}