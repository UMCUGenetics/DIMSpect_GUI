###### DIMSpect GUI ui part ######

ui <- fluidPage(
  titlePanel(title = span(img(src = "DIMSpect_logo.png", height = 30), " DIMS results database")),
  # add a spinner which is activated when a process takes more than 500 ms to replace the progress bar
  # add_busy_spinner(spin = "folding-cube", position = "full-page", timeout = 500), 
  navbarPage(title = "DIMSpect", 
             
    ############ UI: Patient query #######################
    tabPanel("Patient query",
             sidebarLayout(
                sidebarPanel(
                  selectInput(inputId = "matrix_selected", label = "Matrix", 
                              choices = c("plasma","DBS","urine", "CSF"), selected = 1),
                  textInput(inputId = "patient_id", label = "Patient ID"),
                  textInput(inputId = "sample_id", label = "Sample ID"),
                  textInput(inputId = "run_name", label = "Run name"),
                  actionButton(inputId = "view_patient_selection", label = "View all occurrences of patient ID"),
                  uiOutput("sel_for_query"),
                  actionButton(inputId = "update_view_pat_sel", label = "Werkt nog niet: Exclude entries from query"),
                  uiOutput("sample_selection"),
                  textInput(inputId = "zscore_high", label = "Cut-off for elevated Z-score", value = 2),
                  textInput(inputId = "zscore_low", label = "Cut-off for decreased Z-score", value = -1.5),
                  selectInput(inputId = "identified_only", label = "Identified or unidentified peak groups", 
                              choices = c("Identified only", "Unidentified only", "both identified and unidentified"), selected = 1),
                  actionButton(inputId = "run_patient_query", label = "Run patient query")
                ), # end sidebarPanel
                        
                mainPanel(
                  tableOutput(outputId = 'patient_list'),
                  # plotOutput('violin_plot'),
                  plotOutput('bar_plot'),
                  # tableOutput(outputId = 'patient_data_table')
                  DT::dataTableOutput("patient_data_table")
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