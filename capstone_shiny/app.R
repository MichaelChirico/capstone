library(shiny)

ui <- shinyUI(fluidPage(
  tabsetPanel(
    tabPanel("About",
             h3("About Me"),
             p(paste(
               "I am Michael Chirico, a PhD Economist and",
               "Data Incubator Fellow with interests in",
               "applications of data analysis and machine",
               "learning to public policy and good governance.",
               "You can find the code for this project on my GitHub",
               "(MichaelChirico/capstone)."
             )),
             hr(),
             h3("About the Project"),
             p(paste(
               "This project aims to forecast weekly hotspots of",
               "relatively high concentrations of fire events",
               "in the City of Seattle. Sparse geospatial event",
               "prediction is well-known to be a difficult endeavor.", 
               "To tackle this problem for emergency events, I deploy an",
               "autoregressive model of fire event intensity",
               "including lagged kernel density and random Fourier",
               "features. A given gridding of the city is fed into",
               "a large-scale regularized Poisson regression fit",
               "by Vowpal Wabbit. This model allows flexibility in",
               "forecasted cell shape, orientation, RFF expansion",
               "density, and more; in all, 12 hyperparameters",
               "are discovered using Bayesian Optimization."
             )),
             p(paste(
               "The data used in this process consists of",
               "geospatially tagged and categorized 911 events",
               "provided in a log by the Seattle Fire Department.",
               "This data was filtered to exclude events lying",
               "outside Seattle's city limits and to focus on",
               "fire-related events (the data also includes",
               "medical emergencies). A few events in the log were",
               "also missing geospatial information or had",
               "mal-formed time labels. In total, there are",
               "roughly 90,000 fire events in 7 years of data."
             )),
             p(paste(
               "Note: This technique was debuted by Team",
               "Kernel Glitches (myself, Seth Flaxman, Charles",
               "Loeffler, and Pau Pereira-Batlle) in our",
               "submission to the NIJ Real-Time Crime",
               "Forecasting competition, where we tied for the",
               "best performance in the Large Team category.",
               "You can find the code associated thereto on my",
               "GitHub page at MichaelChirico/portland."
             ))),
    tabPanel("Fires in Seattle",
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
))

# server simply runs traj_yr
server <- shinyServer(function(input, output) {
  output$trajectory <- renderPlot(traj_yr(input$yr))
  output$pred_full <- renderImage(
    list(src = paste0('./images/pred_all_', input$wk, '.png'),
         alt = "Predictions for All of Seattle"),
    deleteFile = FALSE
  )
  output$actu_full <- renderImage(
    list(src = paste0('./images/actu_all_', input$wk, '.png'),
         alt = "Actual Hotspots for All of Seattle"),
    deleteFile = FALSE
  )
  output$pred_down <- renderImage(
    list(src = paste0('./images/pred_dwn_', input$wk, '.png'),
         alt = "Predictions for Downtown Seattle"),
    deleteFile = FALSE
  )
  output$actu_down <- renderImage(
    list(src = paste0('./images/actu_dwn_', input$wk, '.png'),
         alt = "Actual Hotspots for Downtown Seattle"),
    deleteFile = FALSE
  )
})

# Run the application 
shinyApp(ui = ui, server = server)
