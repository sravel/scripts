mymtcars = head(mtcars)
for_pop_up = 1:6

app <- shinyApp(
  ui = fluidPage(
    
    DT::dataTableOutput("mydatatable")
  ),
  
  
  server =  shinyServer(function(input, output, session) {
    
    mycars = head(mtcars)
    output$mydatatable = DT::renderDataTable(mycars, selection = 'single',  
                                             rownames = FALSE, options = list(dom = 't'))
    
    observeEvent(input$mydatatable_rows_selected,
                 {
                   showModal(modalDialog(
                     title = "You have selected a row!",
                     mycars[input$mydatatable_rows_selected,]
                   ))
                 })
    
    
    
  })
)