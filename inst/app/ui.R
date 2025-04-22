

ui <- dashboardPage(
  skin = "black",
  title = "eRNAkit",

  # HEADER ------------------------------------------------------------------

  dashboardHeader(
    title = "eRNAkit",
    titleWidth = 300,
    dropdownMenu(
      type = "notifications",
      headerText = strong("DB"),
      icon = icon("question"),
      badgeStatus = NULL,
      notificationItem(
        text = (" DB description"),
        icon = icon("spinner")
      ),
      notificationItem(
        text = paste(db$data[1], db$information[1], sep = " : "),
        icon = icon("address-card")
      ),
      notificationItem(
        text = paste(db$data[2], db$information[2], sep = " : "),
        icon = icon("calendar")
      ),
      notificationItem(
        text = paste(db$data[3], db$information[3], sep = " : "),
        icon = icon("calendar")
      ),
      notificationItem(
        text = paste(db$data[4], db$information[4], sep = " : "),
        icon = icon("user-md")
      ),
      notificationItem(
        text = paste(db$data[5], db$information[5], sep = " : "),
        icon = icon("user-md")
      ),

      notificationItem(
        text = strong(paste(db$data[6], "Information on the db's", sep = " : ")),
        icon = icon("exclamation")
      )
    ),
    tags$li(
      a(
        strong("ABOUT eRNAkit"),
        height = 40,
        href = "https://github.com/caanene1/RBPInper/blob/main/README.md",
        title = "",
        target = "_blank"
      ),
      class = "dropdown"
    )
  ),

  # SIDEBAR -----------------------------------------------------------------

  dashboardSidebar(
    width = 300,
    ## Need to work out how to use this for Demo
    introBox(data.step = 3, data.intro = db$information[1],  # intro tour
             div(class = "inlay", style = "height:15px;width:100%;background-color: #ecf0f5;"),
             sidebarMenu(
               introBox(data.step = 1, data.intro = db$information[3],  # intro tour
                        div(id = "sidebar_button",
                            bsButton(inputId = "demo",
                                     label = "MODULES",
                                     icon = icon("play-circle"),
                                     style = "danger")
                        )
               ),

               div(class = "inlay", style = "height:15px;width:100%;background-color: #ecf0f5;"),

               # New Menu Items for search options
               menuItem(
                 "by eRNA",
                 tabName = "searchByERNA",
                 icon = icon("search"),
                 textInput(inputId = "byERNA", label = "Enter eRNA:", placeholder = "e.g., en1000"),
                 actionButton(inputId = "byERNAButton", label = "Search")
               ),

               menuItem(
                 "by Coordinate",
                 tabName = "searchByCoordinate",
                 icon = icon("search"),
                 textInput(inputId = "byCoordinate", label = "Enter Coordinate:", placeholder = "e.g., 1:10550-10835"),
                 actionButton(inputId = "byCoordinateButton", label = "Search")
               ),

               menuItem(
                 "by Gene",
                 tabName = "searchByGenes",
                 icon = icon("search"),
                 textInput(inputId = "byGene", label = "Enter Ensemble Gene ID:", placeholder = "e.g., ENSG00000142611"),
                 actionButton(inputId = "byGeneButton", label = "Search")
               ),

               menuItem(
                 "Filter Status",
                 tabName = "filterStatus",
                 icon = icon("info-circle")
               ),

               tags$div(
                 style = "padding: 10px; font-size: 14px; color: darkred;",
                 textOutput("filterStatusMessage")
               ),

               br(), br(), br(),

               div(id = "download_button",
                   downloadButton("downloadData",
                                  label = "DOWNLOAD",
                                  class = "btn btn-danger",
                                  style = "width: 100%; margin-top: 10px; margin-bottom: 15px;")
               )
             )
    )
  )  ,


  # BODY --------------------------------------------------------------------

  dashboardBody(
    tags$head(
      tags$link(
        rel = "stylesheet",
        type = "text/css",
        href = "radar_style.css")
    ),

    useShinyjs(),
    introjsUI(),

    # MAIN BODY ---------------------------------------------------------------
    fluidRow(
      column(width = 2, bsButton("loc", label = "LOCALISATION", icon = icon("spinner", class = "spinner-box"), style = "success")),
      column(width = 2, bsButton("r2r", label = "TARGET", icon = icon("spinner", class = "spinner-box"), style = "success")),
      column(width = 2, bsButton("ee", label = "EXPRESSION", icon = icon("user"), style = "success")),
      column(width = 2),
      column( width = 2, bsButton("pacy", label = "RIBSOME", icon = icon("flask", class = "flask-box"), style = "success")),
      column(width = 2, bsButton("ko", label = "STRESS eRNA", icon = icon("flask", class = "flask-box"), style = "success"))
    ),


    # This matches the dynamic server side
    fluidRow(
      column(width = 12, uiOutput("mainPanelUI"))
    )


  ))


