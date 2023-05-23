#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#' @importFrom ggplot2 geom_vline
#' 
ptmotif_ui <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    tabPanel(
      "Ratio selected/all", 
      fluidRow(
        column(6, plotOutput(ns("plt"), height = "400px")),
        column(6, omicsViewer:::dataTableDownload_ui(ns("seqtable_rat")))
      )
    ),
    tabPanel(
      "Selected", 
      fluidRow(
        column(6, plotOutput(ns("plt.fg"), height = "400px")),
        column(6, omicsViewer:::dataTableDownload_ui(ns("seqtable_fg")))
      )
    ),
    tabPanel(
      "All", 
      fluidRow(
        column(6, plotOutput(ns("plt.bg"), height = "400px")),
        column(6, omicsViewer:::dataTableDownload_ui(ns("seqtable_bg")))
      )
    )
  )
}

ptmotif_module <- function(
  input, output, session, fg.seqs = reactive(), bg.pfm = NULL
) {
  ns <- session$ns
  
  fg.pfm <- reactive({
    req(fg.seqs())
    omicsViewer:::aaFreq(fg.seqs())
  })
  
  logo <- reactive({    
    req(bg.pfm)    
    req(fg.pfm())
    omicsViewer:::motifRF(fg.pfm = fg.pfm(), bg.pfm = bg.pfm)
  })
  
  ##
  output$plt <- renderPlot({
    req( logo() )
    ggseqlogo::ggseqlogo( data = logo() ) + geom_vline(
      xintercept = (ncol(logo())+1)/2, linetype="dashed", color = "orange", size=1.5
    )
  })
  
  output$plt.fg <- renderPlot({
    req( d <- fg.pfm() )
    ggseqlogo::ggseqlogo( data = d ) + geom_vline(
      xintercept = (ncol(d)+1)/2, linetype="dashed", color = "orange", size=1.5
    )
  })
  
  output$plt.bg <- renderPlot({
    req( d <- bg.pfm )
    ggseqlogo::ggseqlogo( data = d ) + geom_vline(
      xintercept = (ncol(d)+1)/2, linetype="dashed", color = "orange", size=1.5
    )
  })
  
  mat2df <- function(x) {
    data.frame(Name = rownames(x), x, stringsAsFactors = FALSE)
  }
  
  callModule(
    omicsViewer:::dataTableDownload_module,
    id = "seqtable_fg", reactive_table = reactive(mat2df(fg.pfm())), prefix = "seqLogoPFM_foreground", pageLength = 10)
  
  callModule(
    omicsViewer:::dataTableDownload_module,
    id = "seqtable_bg", reactive_table = reactive(mat2df(bg.pfm)), prefix = "seqLogoPFM_background", pageLength = 10)
  
  callModule(
    omicsViewer:::dataTableDownload_module,
    id = "seqtable_rat", reactive_table = reactive(mat2df(logo())), prefix = "seqLogoPFM_ratio", pageLength = 10)
  
}