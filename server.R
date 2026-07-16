###### DIMSpect GUI server part ######
# Graphical User Interface for DIMS database v0.2

# load functions
source("utils.R")

options(shiny.maxRequestSize = 300 * 1024^2)
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
library(shinyWidgets)
library(bslib)
library(shinycssloaders)
library(shinybusy)
library(thematic)
library(plotly)

thematic::thematic_shiny()

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

  selected_metabolites_sing_pat_plot <- reactive({
    if (is.null(input$single_pat_table_rows_selected)) {
      1
    } else {
      input$single_pat_table_rows_selected
    }
  })

  observeEvent(input$view_patient_selection, {
    if (input$view_patient_selection > 0) {
      enable("pat_zscore_high")
      enable("pat_zscore_low")
      enable("identified_only")
      enable("run_patient_query")
    } else {
      disable("pat_zscore_high")
      disable("pat_zscore_low")
      disable("identified_only")
      disable("run_patient_query")
    }
  })

  observeEvent(input$run_patient_query, {
    showTab(inputId = "patient_query", target = "Barplot", select = FALSE)
    showTab(inputId = "patient_query", target = "Violinplot", select = FALSE)
    showTab(inputId = "patient_query", target = "Single Patient", select = FALSE)
  })


  patient_data <- reactiveValues()

  observeEvent(input$run_patient_query, {
    show_modal_spinner(text = "\n Loading...")
    patient_data$df <- get_patient_data(api(), patient_query_input$df, input$pat_zscore_high,
                                        input$pat_zscore_low, input$identified_only)
    remove_modal_spinner()
  })

  patient_query_input <- reactiveValues()

  observeEvent(input$view_patient_selection, {
    if (is.null(input$view_patient_selection)) {
      return(NULL)
    } else {
      updateTabsetPanel(session, "patient_query", selected = "Patient info")
      show_modal_spinner(text = "\n Loading...")
      patient_query_input$df <- get_info_on_patient(api(), input$patient_id, input$sample_id,
                                                    input$run_name, input$matrix_selected)
      remove_modal_spinner()
    }
  })

  observeEvent(input$update_view_pat_sel, {
    patient_query_input$df <- patient_query_input$df[input$patient_list_rows_selected, ]
  })

  #### view occurrences of patient in database ####
  output$patient_list <- DT::renderDT({
    validate(
      need(input$patient_id != "", "Give a patient ID")
    )
    DT::datatable(patient_query_input$df, rownames = FALSE, filter = "top")
  })

  observeEvent(input$run_patient_query, {
    updateTabsetPanel(session, "patient_query", selected = "Barplot")
  })

  output$patient_data_table_bar <- DT::renderDT({
    validate(
      need(!is.null(patient_data$df), "")
    )
    df <- patient_data$df
    if (nrow(df) > 0) {
      colnames_new <- change_colnames_patient_table(colnames(df))
    }
    datatable(df,
      options = list(
        deferRender = TRUE,
        scrollY = TRUE,
        scrollX = TRUE,
        scroller = TRUE,
        autoWidth = FALSE,
        columnDefs = list(list(width = "10%", targets = c(3, 12, 13)))
      ),
      rownames = FALSE,
      colnames = colnames_new,
      filter = "top",
      escape = FALSE
    ) %>% formatRound(columns = sapply(df, is.numeric), digits = 4)
  })

  output$table_violin <- DT::renderDT({
    validate(
      need(!is.null(patient_data$df), "")
    )
    df <- patient_data$df
    if (nrow(df) > 0) {
      colnames_new <- change_colnames_patient_table(colnames(df))
    }
    datatable(df,
      options = list(
        deferRender = TRUE,
        scrollY = TRUE,
        scrollX = TRUE,
        scroller = TRUE,
        autoWidth = FALSE,
        columnDefs = list(list(width = "10%", targets = c(3, 12, 13)))
      ),
      rownames = FALSE,
      colnames = colnames_new,
      filter = "top",
      escape = FALSE
    ) %>% formatRound(columns = sapply(df, is.numeric), digits = 4)
  })

  data_bar_plot <- reactive({
    df <- patient_data$df
    metab_rows <- selected_metabolites_barplot()
    plot_data_selected <- select_plot_data(df, metab_rows)
    plot_data_long <- reshape2::melt(plot_data_selected, id.vars = c("HMDB_name", "run_name"))
    plot_data_long <- plot_data_long %>% rename(Sample = variable, Zscore = value)
    plot_data_long$Sample <- gsub("_Zscore", "", plot_data_long$Sample)
    plot_data_long <- na.omit(plot_data_long)
    plot_data_long
  })

  data_violin_plot <- reactive({
    df <- patient_data$df
    metab_selected <- df %>%
      slice(selected_metabolites_violinplot()) %>%
      pull(HMDB_name)
    data_violin <- df %>% filter(HMDB_name %in% metab_selected)
    data_violin <- format_data_violin_plot(data_violin, unique(patient_query_input$df$ID))
    data_violin
  })

  output$bar_plot <- renderPlotly({
    validate(
      need(patient_query_input$df$ID != "", "")
    )
    df <- data_bar_plot()
    ggplotly(create_box_plots(data_bar_plot(), unique(patient_query_input$df$ID)), tooltip = "text") %>%
      layout(
        hoverlabel = list(
          bgcolor = "white",
          font = list(family = "Arial", size = 12)
        )
      )
  })

  output$violin_plot <- renderPlotly({
    validate(
      need(patient_query_input$df$ID != "", "")
    )
    all_samples <- unique(patient_query_input$df$ID)
    ggplotly(create_violin_plots(data_violin_plot()), tooltip = "text") %>%
      layout(
        hoverlabel = list(
          bgcolor = "white",
          font = list(family = "Arial", size = 12)
        )
      )
  })

  output$single_pat_plot <- renderPlotly({
    validate(
      need(!is.null(patient_data$df), "")
    )
    metab_plot_data <- format_data_single_pat_plot(patient_data$df, selected_metabolites_sing_pat_plot(),
                                                   unique(patient_query_input$df$ID))
    ggplotly(create_single_pat_plot(metab_plot_data), tooltip = "text") %>%
      layout(
        hoverlabel = list(
          bgcolor = "white",
          font = list(family = "Arial", size = 12)
        )
      )
  })

  output$single_pat_table <- renderDT({
    validate(
      need(!is.null(patient_query_input$df), "")
    )
    patient_table <- patient_query_input$df
    colnames(patient_table) <- c("sample_id", "pat_id", "run_name", "matrix")
    single_pat_table <<- format_data_single_pat_table(patient_data$df, unique(patient_table$sample_id))
    colnames(single_pat_table) <- change_colnames_single_pat_table(colnames(single_pat_table))
    sketch <- make_layout_single_pat_tab(colnames(single_pat_table))
    datatable(single_pat_table,
      container = sketch,
      rownames = FALSE, escape = FALSE,
      options = list(
        filter = "top",
        deferRender = TRUE,
        scrollY = TRUE,
        scrollX = TRUE,
        scroller = TRUE,
        autoWidth = FALSE
      )
    ) %>% formatRound(columns = sapply(single_pat_table, is.numeric), digits = 4)
  })

  #### METABOLITE QUERY ####
  metab_query_input <- reactiveValues()

  observeEvent(input$view_metab_selection, {
    if (input$view_metab_selection > 0) {
      enable("metab_zscore_high")
      enable("metab_zscore_low")
      enable("run_metab_query")
      enable("metab_adducts")
      enable("metab_adducts_pos")
      enable("metab_adducts_neg")
      enable("dropdown_metab")
      enable("metab_plot_ctrls")
    } else {
      disable("metab_zscore_high")
      disable("metab_zscore_low")
      disable("run_metab_query")
      disable("metab_adducts")
      disable("metab_adducts_pos")
      disable("metab_adducts_neg")
      disable("dropdown_metab")
      disable("metab_plot_ctrls")
    }
  })

  observeEvent(input$view_metab_selection, {
    if (is.null(input$view_metab_selection)) {
      return(NULL)
    } else {
      if (input$metabolite_id != "" & !is_valid_hmdb_id(input$metabolite_id)) {
        show_alert(
          title = "Not a valid HMDB ID",
          text = "Please provide a valid HMDB ID starting with HMDB and ending with 5 or 7 numbers.",
          type = "error"
        )
      } else {
        metab_query_input$df <- get_metab_info(api(), input$metabolite_id, input$metabolite_name)
      }
    }
  })

  observeEvent(input$update_view_metab_sel, {
    metab_query_input$df <- metab_query_input$df[input$metab_list_rows_selected, ]
  })

  output$metab_list <- DT::renderDT({
    validate(
      need((input$metabolite_id != "" || input$metabolite_name != ""), "Give a HMDB ID or HMDB name")
    )
    if (!is.null(metab_query_input$df)) {
      df <- metab_query_input$df %>% select(-c(uuid))
      DT::datatable(df,
        rownames = FALSE,
        colnames = c(
          "HMDB ID", "Name", "Molecular formula", "Theoretical mass",
          "Secondary HMDB IDs", "Description"
        ),
        filter = "top",
        escape = FALSE,
        options = list(
          autoWidth = FALSE
        )
      ) %>% formatRound(columns = sapply(df, is.numeric), digits = 4)
    }
  })

  observeEvent(input$run_metab_query, {
    updateTabsetPanel(session, "metab_query", selected = "Metabolite plot")
  })

  metab_zscores <- reactiveValues()

  observeEvent(input$run_metab_query, {
    show_modal_spinner(text = "\n Loading...")
    metab_zscores$df <- get_metab_zscores(api(), metab_query_input$df, input$metab_zscore_high, input$metab_zscore_low)
    remove_modal_spinner()
  })

  output$metab_zscores_table <- DT::renderDT({
    validate(
      need(!is.null(metab_zscores$df), "")
    )
    df <- metab_zscores$df
    DT::datatable(df,
      rownames = FALSE,
      colnames = c(
        "HMDB ID", "Name", "Adduct", "Sample ID", "Z-score",
        "Run name", "Scanmode", "m/z value"
      ),
      filter = "top",
      escape = FALSE,
    ) %>% formatRound(columns = sapply(df, is.numeric), digits = 4)
  })

  output$metab_zscores_plot <- renderPlotly({
    validate(
      need(!is.null(metab_zscores$df), "")
    )
    df <- metab_zscores$df
    if (!input$metab_plot_ctrls) {
      df <- df %>% filter(!str_detect(sample_id, "^C"))
    }
    if (!is.null(input$metab_adducts_pos) | !is.null(input$metab_adducts_neg)) {
      adducts <- c(input$metab_adducts_pos, input$metab_adducts_neg)
      df <- df %>% filter(adduct %in% adducts)
    }

    ggplotly(get_metab_plot(df), tooltip = "text") %>%
      layout(
        hoverlabel = list(
          bgcolor = "white",
          font = list(family = "Arial", size = 12)
        )
      )
  })
}
