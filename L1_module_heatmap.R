

heatmapUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h3("Interaction network"),
      selectizeInput(
        ns("rtk"), label = "Select receptor tyrosine kinases (RTKs):", multiple = TRUE, choices = rtk), 
      selectizeInput(
        ns("bp"), label = "Select interactors of receptor tyrosine kinases:", multiple = TRUE, choices = NULL),
      actionButton(
        ns("go"), label = "Search"
      )
    ),
    mainPanel(
      h4("Interacting network between the given proteins and RTKs"),
      fluidRow(
        column(6, plotlyOutput(ns("heatmap"), height = 950)),
        column(6, DTOutput(ns("tab")))
      )
    )
  )
}

heatmapServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      
      observe({
        ss <- unique(tab_f_rtk()$bp.gene.name)
        updateSelectizeInput(session = session, inputId = "bp", choices = ss, server = TRUE)
      })
      
      tab_f_rtk <- eventReactive(input$rtk, {
        dat[dat$rtk %in% input$rtk, ]
      })
      
      tab_f_bp <- eventReactive(input$go, {
        req( d <- tab_f_rtk() )
        req( input$bp )
        d[d$bp.gene.name %in% input$bp, ]
      })
      
      output$heatmap <- renderPlotly({
        aaa <<- tab_f_bp()
        hm( tab_f_bp() )
      })
      
      output$tab <- renderDT({
        tab <- tab_f_bp()[, scc2]
        colnames(tab) <- names(scc2)
        formatDTScrollY( tab, height = "800px" )
      })
    }
  )
}
