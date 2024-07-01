###### DIMSpect GUI ui part ######

library(shinyjs)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  # titlePanel(title = span(img(src = "DIMSpect_logo.png", height = 30), " DIMS results database")),
  # add a spinner which is activated when a process takes more than 500 ms to replace the progress bar
  # add_busy_spinner(spin = "folding-cube", position = "full-page", timeout = 500), 
  navbarPage(title = div(img(src="dimspect_logo.svg", 
                             height=50, 
                             style="margin-top: -14px;
                                           padding-right: 10px;
                                           padding-bottom: 10px"),"DIMS results database"), 
             
    ############ UI: Patient query #######################
    tabPanel("Patient query",
             sidebarLayout(
                sidebarPanel(width = 3, 
                  selectInput(inputId = "matrix_selected", label = "Matrix", 
                              choices = c("plasma","DBS","urine", "CSF"), selected = 1),
                  textInput(inputId = "patient_id", label = "Patient ID"),
                  textInput(inputId = "sample_id", label = "Sample ID"),
                  textInput(inputId = "run_name", label = "Run name"),
                  actionButton(inputId = "view_patient_selection", label = "Get all patient info"),
                  br(), br(),
                  actionButton(inputId = "update_view_pat_sel", label = "Exclude entries from query"),
                  hr(style = "border-top: 1px solid #303c54;"),
                  disabled(
                    numericInput(inputId = "zscore_low", label = "Cut-off for decreased Z-score", value = -1.5, step = 0.1)
                    ),
                  disabled(
                    numericInput(inputId = "zscore_high", label = "Cut-off for elevated Z-score", value = 2, step = 0.1)
                    ),
                  disabled(
                    selectInput(inputId = "identified_only", label = "Identified or unidentified peak groups", 
                              choices = c("Identified only" = "iden", "Unidentified only" = "not_iden", 
                                          "both identified and unidentified" = "all"), selected = 1)
                    ),
                  disabled(
                    actionButton(inputId = "run_patient_query", label = "Run patient query")
                  )
                ), # end sidebarPanel
                        
                mainPanel(width = 9,
                  tabsetPanel(type = "tabs", id = "patient_query",
                    tabPanel("Patient info",DT::dataTableOutput('patient_list')),
                    
                    tabPanel("Barplot", 
                             plotOutput('bar_plot', height = "600px"),
                             DT::dataTableOutput("patient_data_table_bar")),
                    tabPanel("Violinplot", 
                             plotOutput("violin_plot", height = "600px"),
                             DT::dataTableOutput("table_violin"))
                  )
                ) #end of mainPanel
                        
              ) #end of sidebarLayout (includes sidebarPanel and mainPanel)
                      
            ), #end of tabPanel "Patient query"
             
             
    #### UI: Metabolite query   #####       
    tabPanel("Metabolite query",
             sidebarLayout(
               sidebarPanel(
                 textInput(inputId = "metabolite_id", label = "HMDB ID")
               ), # end sidebarPanel
               
               mainPanel(
                 uiOutput('someplots'),
               ) #end of mainPanel
             ) #end of sidebarLayout (includes sidebarPanel and mainPanel)
    ), #end of main tab Metabolite query

    #### UI: m/z query   #####       
    tabPanel("Mass over charge query",
             sidebarLayout(
               sidebarPanel(
                 textInput(inputId = "input_mass", label = "mass"),
                 selectInput(inputId = "input_scan_mode", label = "scan mode", 
                             choices = c("positive", "negative", "neutral"), selected = 1)
               ), # end sidebarPanel
               
               mainPanel(
                 uiOutput('sometable'),
               ) #end of mainPanel
             ) #end of sidebarLayout (includes sidebarPanel and mainPanel)
    ), #end of main tab Metabolite query
    
    #### UI: Help tab  ##### 
    tabPanel("Help",
             h2("What is this?"),
             h4("This is a user interface for accessing and inspecting DI-HRMS (DIMS) metabolomics data."),
             h2("How does it work?"),
             h4("To enter a query, use either the Patient tab and enter a patient ID or use the Metabolite tab to search for a metabolite."),
    ) #end of Help tab
             
  ) # end navbarPage(title = "DIMSpect"
  
) # end fluidPage