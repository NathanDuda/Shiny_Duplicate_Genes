library(shiny)
library(plotly)

# Define the UI
ui <- fluidPage(
  titlePanel("Duplicate Gene Analyses"),
  sidebarLayout(
    sidebarPanel(
      textInput("inputPath", "Path to OrthoFinder Output Directory:", ""),
      textInput("expressionPath", "Expression File Path:", ""),
      textInput("annotationPath", "Annotation gtf/gff File Path:", ""),
      textInput("outputPath", "Output File Path:", ""),
      fileInput("file", "Choose a file"),
      checkboxGroupInput("scripts", "Select analyses:",
                         choices = c("Functionalization", "Duplication Mechanism", 'Tissue Specificity', "Selection")
      ),
      actionButton("runButton", "Run Selected Analyses")
    ),
    mainPanel(
      verbatimTextOutput("statusOutput"),
      textOutput("output_text_1"),
      textOutput("output_text_2"),
      textOutput("output_text_3"),
      textOutput("output_text_4"),
      textOutput("output_text_5"),
      #plotOutput("output_plot_func"), #1
      plotlyOutput('output_plot_func'),
      plotlyOutput("output_plot_mech")
    )
  )
)

# Define the server
server <- function(input, output, session) {
  # Variables to store paths and results
  input_path <- NULL
  output_path <- NULL
  result <- NULL
  
  scripts_root <- 'C:/Users/17735/Downloads/Shiny_Duplicate_Genes/scripts/'
  
  # Event handler for the "Run Selected Analyses" button
  observeEvent(input$runButton, {
    req(!is.null(input$inputPath) && !is.null(input$outputPath), "Please submit input and output paths first.")
    
    # Get input and output file paths
    input_path <- input$inputPath
    expression_path <- input$expressionPath
    output_path <- input$outputPath
    annotation_path <- input$annotationPath
    
    source(paste0(scripts_root,'startup.R'))
    source(paste0(scripts_root,'OrthoFinder_Output_Cleanup.R'))
    n_dups <- clean_orthofinder(input_path, output_path)
    result <- run_selected_scripts(input_path, expression_path, output_path, annotation_path, input$scripts)
    
    
    output$output_text_1 <- renderText({paste("A total of ", n_dups, " duplicate pairs were found:")})
    
    

    
    
  })
  
  # Function to run selected scripts
  run_selected_scripts <- function(input_path, expression_path, output_path, annotation_path, selected_scripts) {

    #func_counts <- NULL # allows these to not be chosen 
    #mech_counts <- NULL
    

    # Run selected scripts
    if ("Functionalization" %in% selected_scripts) {
      source(paste0(scripts_root,'Ancestral_Copy.R'))
      get_ancestral_copy(input_path, output_path, expression_path)
      source(paste0(scripts_root,'Functionalization.R'))
      func_counts <- functionalization(output_path, expression_path)
      output$output_text_2 <- renderText({paste(func_counts['Conserved'], " conserved,")})
      output$output_text_3 <- renderText({paste(func_counts['Neofunctionalized'], " neofunctionalized,")})
      output$output_text_4 <- renderText({paste(func_counts['Specialized'], " specialized, and")})
      output$output_text_5 <- renderText({paste(func_counts['Subfunctionalized'], " subfunctionalized.")})
      output$output_plot_func <- renderPlotly({visualize_func(func_counts)})
    }
    if ("Duplication Mechanism" %in% selected_scripts) {
      source(paste0(scripts_root,'Duplication_Mechanism.R'))
      mech_counts <- duplication_mechanism(output_path, annotation_path)
      output$output_plot_mech <- renderPlotly({visualize_mech(mech_counts)})
    }
    if ("Duplication Mechanism" %in% selected_scripts) {
      source(paste0(scripts_root,'Tau_Calculation.R'))
      calculate_tau(expression_path, output_path)
      
    }
    
    if ("Selection" %in% selected_scripts) {
      #run_hyphy_command <- paste("bash ./scripts/HyPhy.sh", shQuote(input_path), shQuote(output_path), collapse = " ")
      #system(run_hyphy_command, intern = TRUE)
    }
    

    return(0)
  }
  
}

# Run the application
shinyApp(ui, server)
