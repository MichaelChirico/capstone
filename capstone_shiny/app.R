library(shiny)

wks = c('20170301', '20170308', '20170315', '20170322')
url_wk = setNames(rep('http', 4L), wks)
URLs = list()
URLs$pred_full = url_wk
url_wk[wks] = ''
URLs$actu_full = url_wk
url_wk[wks] = ''
URLs$pred_down = url_wk
url_wk[wks] = ''
URLs$actu_down = url_wk

ui <- shinyUI(fluidPage(
  titlePanel("Predicted Weekly Fire Locations in Seattle"),
  fluidRow(column(6, plotOutput("pred_full", height = "600px")),
           column(6, plotOutput("actu_full", height = "600px"))),
  fluidRow(column(6, plotOutput("pred_down", height = "600px")),
           column(6, plotOutput("actu_down", height = "600px"))),
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
    list(src = URLs[["pred_full"]][input$wk],
         alt = "Predictions for All of Seattle")
  )
  output$actu_full <- renderImage(
    list(src = URLs[["actu_full"]][input$wk],
         alt = "Actual Hotspots for All of Seattle")
  )
  output$pred_down <- renderImage(
    list(src = URLs[["pred_down"]][input$wk],
         alt = "Predictions for Downtown Seattle")
  )
  output$actu_down <- renderImage(
    list(src = URLs[["actu_down"]][input$wk],
         alt = "Actual Hotspots for Downtown Seattle")
  )
})

# Run the application 
shinyApp(ui = ui, server = server)