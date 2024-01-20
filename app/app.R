library(shiny)

# Define the UI
ui <- fluidPage(
  titlePanel("Script Runner App"),
  sidebarLayout(
    sidebarPanel(
      textInput("inputPath", "Path to OrthoFinder Output Directory:", ""),
      textInput("expressionPath", "Expression File Path:", ""),
      textInput("annotationPath", "Annotation gtf/gff File Path:", ""),
      textInput("outputPath", "Output File Path:", ""),
      fileInput("file", "Choose a file"),
      checkboxGroupInput("scripts", "Select analyses:",
                         choices = c("Functionalization", "Duplication Mechanism", "Selection")
      ),
      actionButton("runButton", "Run Selected Scripts")
    ),
    mainPanel(
      verbatimTextOutput("statusOutput"),
      textOutput("output_text_1"),
      textOutput("output_text_2"),
      textOutput("output_text_3"),
      textOutput("output_text_4"),
      textOutput("output_text_5"),
      plotOutput("output_plot_1"),
      plotOutput("output_plot_2")
    )
  )
)

# Define the server
server <- function(input, output, session) {
  # Variables to store paths and results
  input_path <- NULL
  output_path <- NULL
  result <- NULL
  
  # Event handler for the "Run Selected Scripts" button
  observeEvent(input$runButton, {
    req(!is.null(input$inputPath) && !is.null(input$outputPath), "Please submit input and output paths first.")
    
    # Get input and output file paths
    input_path <- input$inputPath
    expression_path <- input$expressionPath
    output_path <- input$outputPath
    annotation_path <- input$annotationPath
    
    source('C:/Users/17735/Downloads/Shiny_Duplicate_Genes/scripts/startup.R')
    source('C:/Users/17735/Downloads/Shiny_Duplicate_Genes/scripts/OrthoFinder_Output_Cleanup.R')
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
      source('C:/Users/17735/Downloads/Shiny_Duplicate_Genes/scripts/Ancestral_Copy.R')
      get_ancestral_copy(input_path, output_path, expression_path)
      source('C:/Users/17735/Downloads/Shiny_Duplicate_Genes/scripts/Functionalization.R')
      func_counts <- functionalization(output_path, expression_path)
      output$output_text_2 <- renderText({paste(func_counts['Conserved'], " conserved,")})
      output$output_text_3 <- renderText({paste(func_counts['Neofunctionalized'], " neofunctionalized,")})
      output$output_text_4 <- renderText({paste(func_counts['Specialized'], " specialized, and")})
      output$output_text_5 <- renderText({paste(func_counts['Subfunctionalized'], " subfunctionalized.")})
      output$output_plot_1 <- renderPlot({visualize_func(func_counts)})
    }
    if ("Duplication Mechanism" %in% selected_scripts) {
      source('C:/Users/17735/Downloads/Shiny_Duplicate_Genes/scripts/Duplication_Mechanism.R')
      mech_counts <- duplication_mechanism(output_path, annotation_path)
      output$output_plot_2 <- renderPlot({visualize_mech(mech_counts)})
    }
    if ("Selection" %in% selected_scripts) {
      #run_hyphy_command <- paste("bash ./scripts/HyPhy.sh", shQuote(input_path), shQuote(output_path), collapse = " ")
      #system(run_hyphy_command, intern = TRUE)
    }
    
    # func_counts wont exist if Functionalization not chosen
    
    return(0)
  }
  
}

# Run the application
shinyApp(ui, server)
