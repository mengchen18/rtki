ui <- fluidPage(
  h1(tt),
  tabsetPanel(
    id = "main",

    tabPanel(
      "Search a receptor tyrosin kinase", 
      pyUI("pymod")
    ),
    tabPanel(
      "Search an RTK interacting protein",
      bpUI("bpmod")
    ),
    tabPanel(
      "Interaction network",
      heatmapUI("heatmap")
    )
  )
)


