pyUI <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    
    sidebarPanel(
      style = "background: white", #; height: 800px",
      tags$head(tags$style(HTML(sprintf("#%s {height: 550px}", ns("sites"))))), # ; font-size: 50px;
      selectizeInput(ns("rtk"), label = "Select a receptor tyrosine kinase (RTK):", multiple = FALSE, choices = rtk), 
      multiInput(ns("sites"), label = "Select pY residues", choices = "", width = "100%") 
    ),
    mainPanel(
      wellPanel(
        style = "background: white",
        fluidRow(
          strong(h4("RTK interacting proteins with SH2 or PTB domain:")),
          column(6, plotOutput(ns("linearSeq"), height = "600px")),
          column(6, DTOutput(ns("domainTab"))),
          strong(h4("All RTK interacting proteins:")),
          column(12, DTOutput(ns("allTab")))
        )
      )
    )
  )
}

pyServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      
      siteOpts <- eventReactive(input$rtk, {
        pySites[[input$rtk]]
      })
      
      observe({
        req(siteOpts())
        updateMultiInput(session, "sites", choices = siteOpts())
      })
      
      at <- reactive({
        req(input$rtk)
        d0 <- dat[dat$rtk %in% input$rtk, ]
        d0 <- d0[order(d0$af.filter.batch, decreasing = TRUE), ]
        if (length(input$sites) > 0)
          d0 <- d0[d0$pY.id2 %in% input$sites, ]
        d0
      })
      
      plotPrep <- reactive({
        req(input$rtk)
        u <- NULL
        if (length(input$sites) > 0)
          u <- input$sites
        drawSeqPrep(pi[[input$rtk]], dat = dat_domain, pY = u)
        
      })
      
      output$linearSeq <- renderPlot({
        drawSeq(plotPrep())
      })
      
      output$domainTab <- renderDT({
        t1 <- plotPrep()$df2[, scc]
        t1 <- t1[order(t1$af.filter.batch, decreasing = TRUE), ]
        colnames(t1) <- names(scc)
        formatDTScrollY(t1)
      })
      
      output$allTab <- renderDT({
        t1 <- at()[, scc2]
        colnames(t1) <- names(scc2)
        formatDTScrollY(t1, height = "450px")
        })
    }
  )
}
