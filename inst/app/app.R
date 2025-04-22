source(system.file("app", "dependencies.R", package = "eRNAkit"))
source(system.file("app", "global.R", package = "eRNAkit"))
source(system.file("app", "ui.R", package = "eRNAkit"))
source(system.file("app", "server.R", package = "eRNAkit"))

shinyApp(ui=ui, server=server)
