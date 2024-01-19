library(shiny)

# Define the UI
ui <- fluidPage(
  titlePanel("Script Runner App"),
  sidebarLayout(
    sidebarPanel(
      textInput("inputPath", "Input File Path:", ""),
      textInput("expressionPath", "Expression File Path:", ""),
      textInput("outputPath", "Output File Path:", ""),
      fileInput("file", "Choose a file"),
      checkboxGroupInput("scripts", "Select analyses:",
                         choices = c("Functionalization", "Duplication Mechanism", "Selection")
      ),
      actionButton("runButton", "Run Selected Scripts")
    ),
    mainPanel(
      #verbatimTextOutput("statusOutput"),
      textOutput("output_text_1"),
      textOutput("output_text_blank"),
      textOutput("output_text_2"),
      textOutput("output_text_3"),
      textOutput("output_text_4"),
      textOutput("output_text_5"),
      textOutput("output_text_blank"),
      plotOutput("outputPlot"),
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
    req(!is.null(input$inputPath) && !is.null(input$outputPath), "Please submit both paths first.")
    
    # Get input and output file paths
    input_path <- input$inputPath
    expression_path <- input$expressionPath
    output_path <- input$outputPath
    
    source('C:/Users/17735/Downloads/Shiny_Duplicate_Genes/scripts/startup.R')
    source('C:/Users/17735/Downloads/Shiny_Duplicate_Genes/scripts/OrthoFinder_Output_Cleanup.R')
    n_dups <- clean_orthofinder(input_path, output_path)
    func_counts <- run_selected_scripts(input_path, expression_path, output_path, input$scripts)
    
    
    output$output_text_1 <- renderText({paste("A total of ", n_dups, " duplicate pairs were found:")})
    output$output_text_blank <-  renderText({paste('')})
    output$output_text_2 <- renderText({paste(func_counts['Conserved'], " conserved,")})
    output$output_text_3 <- renderText({paste(func_counts['Neofunctionalized'], " neofunctionalized,")})
    output$output_text_4 <- renderText({paste(func_counts['Specialized'], " specialized, and")})
    output$output_text_5 <- renderText({paste(func_counts['Subfunctionalized.'], " subfunctionalized.")})
    
    output$outputPlot <- renderPlot({
      func_counts <- as.data.frame(func_counts)
      ggplot(func_counts, aes(x="", y=Freq, group=unique(Var1), fill=Var1)) +
        geom_bar(width = 1, stat = "identity", position = position_fill()) +
        geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
        theme_void() +
        theme(legend.title = element_blank(), panel.background = element_rect(fill = "white", color = "white"), legend.background = element_rect(fill = "white", color = 'white')) +
        coord_polar("y") +
        scale_fill_manual(values=c("#3B7BBD", "#E23F51", "#6CBC4D",'#F18244'))
      
    })
    
  })
  
  # Function to run selected scripts
  run_selected_scripts <- function(input_path, expression_path, output_path, selected_scripts) {

    
    
    # Run selected scripts
    if ("Functionalization" %in% selected_scripts) {
      source('C:/Users/17735/Downloads/Shiny_Duplicate_Genes/scripts/Ancestral_Copy.R')
      get_ancestral_copy(output_path, expression_path)
      source('C:/Users/17735/Downloads/Shiny_Duplicate_Genes/scripts/Functionalization.R')
      func_counts <- functionalization(output_path, expression_path)
      
    }
    
    if ("Duplication Mechanism" %in% selected_scripts) {
      # Run Script 2
      # Example: results$script2 <- another_function(data)
    }
    
    if ("Selection" %in% selected_scripts) {
      # Run Script 3
      # Example: results$script3 <- yet_another_function(data)
    }
    
    # func_counts wont exist if Functionalization not chosen
    
    return(func_counts)
  }
  
  # Save the result to a file when the app is closed
}

# Run the application
shinyApp(ui, server)
