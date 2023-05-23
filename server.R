server <- function(input, output, session) {
  
  showModal(modalDialog(
    title = tt, size = "l",
    radioGroupButtons(
      inputId = "func",
      label = "What to do?",
      choices = c(
        "Given a receptor tyrosine kinase (RTK), exploring the interacting proteins of each phosphotyrosine in the given RTK" = "s_py", 
        "Given a protein, exploring its interacting receptor tyrosine kinases (RTKs) and phosphotyrosines" = "s_bp", 
        "Given multiple receptor tyrosine kinases (RTKs) and proteins, exploring the interacting network between the given proteins and RTKs" =  "s_heat"),
      direction = "vertical", width = "100%", selected = character(0)
    ),
    footer = NULL
  ))
  
  observe({
    req(input$func)
    val <- switch (
      input$func,
      "s_py" = "Search a receptor tyrosin kinase", 
      "s_bp" = "Search an RTK interacting protein", 
      "s_heat" = "Interaction network"
    )
    updateTabsetPanel(inputId = "main", selected = val)
    removeModal()
  })
  
  pyServer("pymod")
  bpServer("bpmod")
  heatmapServer("heatmap")
}