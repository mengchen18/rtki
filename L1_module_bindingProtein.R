

bpUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      style = "background: white;", # height: 800px",
      selectizeInput(
        ns("bp"), label = "Select an interactor of receptor tyrosine kinase (RTK):", multiple = FALSE, choices = NULL)
    ),
    mainPanel(
      
      wellPanel(        
        style = "background: white",
        h4("pY sites interacting with the given protein (Select multiple points to perform motif analysis)."),
        fluidRow(
          column(6, plotlyOutput(ns('dotscore'), height = "600px")),
          column(6, DTOutput(ns("allTab"))),
          uiOutput(ns("motif_ui"))
        )
      )
    )
  )
}


bpServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- NS(id)
      updateSelectizeInput(session = session, inputId = "bp", choices = bp, server = TRUE)
      
      at <- reactive({
        req(input$bp)
        d0 <- dat[dat$bp.gene.name %in% input$bp, ]
        d0[order(d0$af.filter.batch, decreasing = TRUE), ]
      })
      
      selected <- reactiveVal(NULL)
      observeEvent(input$bp, {
        selected(NULL)
      })
      
      observe({
        selected ( event_data("plotly_selected", source = "selected_rtk_pY") )
      })
      
      fg <- reactive({
        req(selected())
        req(at())
        at()$pY.sequence[selected()$x]
      })
      
      tab <- reactive({
        if (is.null(selected()))
          return(at())
        at()[selected()$x, ]
      })
      
      output$dotscore <- renderPlotly({
        req(v <- at())
        v$rank <- seq_len(nrow(v))
        fig <- plot_ly(data = v, source = "selected_rtk_pY")
        fig <- add_trace(
          fig, x = ~ rank, y = ~ af.filter.batch, 
          marker = list(sizemode = 'diameter', opacity = 0.75), type = "scatter", mode = "markers",
          text = ~ pY.id, hoverinfo = 'text', showlegend = FALSE)
        layout(
          fig,
          xaxis = list(title = 'Rank of pY sites', range = list(0, nrow(v)+1)),
          yaxis = list(title = 'Enrichment score')
        )
      })
      
      output$allTab <- renderDT({
        iii <- c(3, 4, 2, 5)
        t1 <- tab()[, scc2[iii]]
        t1 <- t1[order(t1$af.filter.batch, decreasing = TRUE), ]
        colnames(t1) <- names(scc2[iii])
        formatDTScrollY(t1)
        
        })
      
      v7 <- callModule(ptmotif_module, id = "motif", fg.seqs = fg, bg.pfm = motif_bg)
      
      output$motif_ui <- renderUI({
        req(fg())
        tagList(
          h4("Motif analysis using pY sequence selected above."),
          ptmotif_ui(ns("motif"))
        )
      })
    }
  )
}
