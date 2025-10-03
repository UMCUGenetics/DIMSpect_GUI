###### DIMSpect GUI ui part ######

library(shinyjs)
library(bslib)
library(shinyWidgets)
library(plotly)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  theme = bs_theme(preset = "flatly"),
  navbarPage(
    title = div(img(
      src = "dimspect_logo.svg",
      height = 50,
      style = "padding-right: 10px;"
    )),
    windowTitle = "DIMSpect",

    ############ UI: Patient query #######################
    tabPanel(
      "Patient query",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          selectInput(
            inputId = "matrix_selected", label = "Matrix",
            choices = c("plasma", "DBS", "urine", "CSF"), selected = 1
          ),
          textInput(inputId = "patient_id", label = "Patient ID"),
          textInput(inputId = "sample_id", label = "Sample ID"),
          textInput(inputId = "run_name", label = "Run name"),
          actionButton(inputId = "view_patient_selection", label = "Get all patient info"),
          br(), br(),
          actionButton(inputId = "update_view_pat_sel", label = "Include entries for search"),
          hr(style = "border-top: 1px solid #303c54;"),
          disabled(
            numericInput(inputId = "pat_zscore_low", label = "Cut-off for decreased Z-score", value = -1.5, step = 0.1)
          ),
          disabled(
            numericInput(inputId = "pat_zscore_high", label = "Cut-off for elevated Z-score", value = 2, step = 0.1)
          ),
          disabled(
            selectInput(
              inputId = "identified_only", label = "Identified or unidentified peak groups",
              choices = c(
                "Identified only" = "iden", "Unidentified only" = "not_iden",
                "both identified and unidentified" = "all"
              ), selected = 1
            )
          ),
          disabled(
            actionButton(inputId = "run_patient_query", label = "Run patient query")
          )
        ), # end sidebarPanel
        mainPanel(
          width = 9,
          tabsetPanel(
            type = "tabs", id = "patient_query",
            tabPanel("Patient info", DT::dataTableOutput("patient_list")),
            tabPanel(
              "Barplot",
              plotlyOutput("bar_plot"),
              hr(style = "border-top: 1px solid #303c54;"),
              DT::dataTableOutput("patient_data_table_bar")
            ),
            tabPanel(
              "Violinplot",
              plotlyOutput("violin_plot", height = "600px"),
              hr(style = "border-top: 1px solid #303c54;"),
              DT::dataTableOutput("table_violin")
            ),
            tabPanel(
              "Single Patient",
              plotlyOutput("single_pat_plot", height = "600px"),
              hr(style = "border-top: 1px solid #303c54;"),
              DT::DTOutput("single_pat_table")
            )
          )
        ) # end of mainPanel
      ) # end of sidebarLayout (includes sidebarPanel and mainPanel)
    ), # end of tabPanel "Patient query"

    #### UI: Metabolite query   #####
    tabPanel(
      "Metabolite query",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          textInput(inputId = "metabolite_id", label = "HMDB ID"),
          textInput(inputId = "metabolite_name", label = "HMDB name"),
          actionButton(inputId = "view_metab_selection", label = "Get all metabolites"),
          br(), br(),
          actionButton(inputId = "update_view_metab_sel", label = "Include entries for search"),
          hr(style = "border-top: 1px solid #303c54;"),
          disabled(
            numericInput(inputId = "metab_zscore_low", label = "Cut-off for decreased Z-score", value = -1.5, step = 0.1)
          ),
          disabled(
            numericInput(inputId = "metab_zscore_high", label = "Cut-off for elevated Z-score", value = 2, step = 0.1)
          ),
          br(),
          disabled(
            actionButton(inputId = "run_metab_query", label = "Run metabolite query")
          )
        ), # end sidebarPanel
        mainPanel(
          width = 9,
          tabsetPanel(
            type = "tabs", id = "metab_query",
            tabPanel(
              "Metabolite info",
              DT::dataTableOutput("metab_list")
            ),
            tabPanel(
              "Metabolite plot",
              disabled(
                dropdownButton(tags$h4("Change plot"),
                  circle = TRUE,
                  status = "primary",
                  icon = icon("gear"), width = "300px",
                  tooltip = tooltipOptions(title = "Change plot inputs"),
                  inputId = "dropdown_metab",
                  materialSwitch("metab_plot_ctrls", "Show Controls",
                    status = "primary", value = TRUE
                  ),
                  fluidRow(
                    column(
                      width = 5,
                      checkboxGroupInput(
                        label = "Positive adducts",
                        inputId = "metab_adducts_pos",
                        choices = c(
                          "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NaCl]+", "[M+NH4]+", "[M+2Na-H]+",
                          "[M+CH3OH]+", "[M+KCl]+", "[M+NaK-H]+"
                        )
                      )
                    ),
                    column(
                      width = 5,
                      checkboxGroupInput(
                        label = "Negative adducts",
                        inputId = "metab_adducts_neg",
                        choices = c(
                          "[M-H]-", "[M+Cl]-", "[M+For]-", "[M+NaCl]-", "[M+KCl-", "[M+H2PO4]-",
                          "[M+HSO4]-", "[M+Na-H]-", "[M+K-H]-", "[M-H2O]-", "[M-2H]-",
                          "[M+I]-", "[M+Ac]-"
                        )
                      )
                    )
                  )
                )
              ),
              plotlyOutput("metab_zscores_plot"),
              DT::dataTableOutput("metab_zscores_table")
            )
          )
        ) # end of mainPanel
      ) # end of sidebarLayout (includes sidebarPanel and mainPanel)
    ), # end of main tab Metabolite query

    #### UI: m/z query   #####
    tabPanel(
      "Mass over charge query",
      sidebarLayout(
        sidebarPanel(
          textInput(inputId = "input_mass", label = "mass"),
          selectInput(
            inputId = "input_scan_mode", label = "scan mode",
            choices = c("positive", "negative", "neutral"), selected = 1
          )
        ), # end sidebarPanel

        mainPanel(
          uiOutput("sometable"),
        ) # end of mainPanel
      ) # end of sidebarLayout (includes sidebarPanel and mainPanel)
    ), # end of main tab Metabolite query

    #### UI: Help tab  #####
    tabPanel(
      "Help",
      h4("What is this?"),
      p("This is a user interface for accessing and inspecting DI-HRMS (DIMS) metabolomics data."),
      h4("How does it work?"),
      p("To enter a query, use either the Patient tab and enter a patient ID or use the Metabolite tab to search for a metabolite."),
    ) # end of Help tab
  ) # end navbarPage(title = "DIMSpect"
) # end fluidPage
