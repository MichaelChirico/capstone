library(shiny)

ui <- shinyUI(fluidPage(
  titlePanel("Predicted Weekly Fire Locations in Seattle"),
  h4("Predicted vs. Actual Intensity for the Whole City"),
  fluidRow(column(6, plotOutput("pred_full", height = "500px")),
           column(6, plotOutput("actu_full", height = "500px"))),
  hr(),
  h4("Predicted vs. Actual Intensity for Downtown Seattle"),
  fluidRow(column(6, plotOutput("pred_down", height = "500px")),
           column(6, plotOutput("actu_down", height = "500px"))),
  hr(),
  fluidRow(column(12, selectInput(
    "wk", "Week:", width = "100%",
    choices = list(
      'March 1-7, 2017' = '20170301',
      'March 8-14, 2017' = '20170308',
      'March 15-21, 2017' = '20170315',
      'March 22-28, 2017' = '20170322'
    )
  )))
))

# server simply runs traj_yr
server <- shinyServer(function(input, output) {
  output$trajectory <- renderPlot(traj_yr(input$yr))
  output$pred_full <- renderImage(
    list(src = paste0('./images/pred_all_', input$wk, '.png'),
         alt = "Predictions for All of Seattle")
  )
  output$actu_full <- renderImage(
    list(src = paste0('./images/actu_all_', input$wk, '.png'),
         alt = "Actual Hotspots for All of Seattle")
  )
  output$pred_down <- renderImage(
    list(src = paste0('./images/pred_dwn_', input$wk, '.png'),
         alt = "Predictions for Downtown Seattle")
  )
  output$actu_down <- renderImage(
    list(src = paste0('./images/actu_dwn_', input$wk, '.png'),
         alt = "Actual Hotspots for Downtown Seattle")
  )
})

# Run the application 
shinyApp(ui = ui, server = server)
