# Standalone Shiny entry point for nanoamp.
suppressPackageStartupMessages(library(nanoamp))
shiny::shinyApp(
  ui = nanoamp:::nanoamp_gui_ui(),
  server = nanoamp:::nanoamp_gui_server
)
