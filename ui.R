library(shiny)
library(shinycssloaders)
library(shinyWidgets)

# Define UI for application that draws a histogram
fluidPage(
  useShinyjs(),
  
  navbarPage("DIMSdb",
    tabPanel("Patients",
      sidebarLayout(
        sidebarPanel(
          shinyjs::useShinyjs(),
          textInput("pat_samp_id", label = "Patient or sample id"),
          selectInput("matrix_type", label = "Select material", choices = list("Dried blood spot"="DBS","Urine"="U","Plasma"="P"), selected = "BS"),
          radioButtons("sample_set", label = "All samples?", choices = list("Yes"="Y","No"="N"), selected = "Y"),
          radioButtons("diff_zscore",label = "Select ", choices = list("Increased"="increase","Decreased"="decrease","All"="all")),
          numericInput("pos_zscore", label = "Select positive Z-score cut-off", value = 2, min = 0, step = 0.5),
          numericInput("neg_zscore", label = "Select negative Z-score cut-off", value = -2.5, max = 0, step = 0.5),
          actionButton("GO_button", label = "GO"),
          ),
        mainPanel(
          plotOutput("test"),
          uiOutput("minus"),
          uiOutput("plus"),
          br(),
          textOutput("value"),
        )
      )
    )
  )
)
