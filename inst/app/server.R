
server <- function(input, output, session) {


  # DEFINE SETS -------------------------------------------------

  # Store current selected view (button)
  selected_view <- reactiveVal(NULL)

  #show intro modal
  observeEvent("", {
    showModal(modalDialog(
      includeHTML(system.file("app", "intro_text.html",
                              package = "eRNAkit")),
      easyClose = TRUE,
      footer = tagList(
        actionButton(inputId = "intro", label = "CLOSE")
      )
    ))
  })

  observeEvent(input$intro,{
    removeModal()
  })


  # DEFINE Search process -------------------------------------------------

  # Filtered place holder
  filtered_db <- reactiveVal(NULL)
  filter_status <- reactiveVal("Search to start DB filtering ...")
  ## First status update
  output$filterStatusMessage <- renderText({
    filter_status()
  })

  # Search by eRNA
  observeEvent(input$byERNAButton, {
    req(input$byERNA)
    dRes <- eRNAkit::subDB(emi, "E",
                          eRNAkit::stringS(input$byERNA))

    if (eRNAkit::check_null(dRes)) {
      showNotification("No matching eRNAs found.", type = "warning")
      filtered_db(NULL)
      return()
    }

    dRes[["subloc"]] <- eRNAkit::winLoc(dRes$loc)
    filtered_db(dRes)
    filter_status("DB filtering done! You can now proceed.")
  })

  # Search by Coordinates
  observeEvent(input$byCoordinateButton, {
    req(input$byCoordinate)

    ehits <- eRNAkit::byRange(eRNAkit::getGRange(emi$core),
                     input$byCoordinate)

    if (nrow(ehits) == 0) {
      showNotification("No matching coordinates found.", type = "warning")
      filtered_db(NULL)
      return()
    }

    dRes <- eRNAkit::subDB(emi, "E", unique(ehits$ID))
    dRes[["subloc"]] <- eRNAkit::winLoc(dRes$loc)
    filtered_db(dRes)
    filter_status("DB filtering done! You can now proceed.")
  })

  # Search by Gene
  observeEvent(input$byGeneButton, {
    req(input$byGene)
    eRs <- eRNAkit::subDB(emi, "G", eRNAkit::stringS(input$byGene))

    if (eRNAkit::check_null(eRs)) {
      showNotification("No matching Gene pair found.", type = "warning")
      filtered_db(NULL)
      return()
    }

    eRs <- unique(eRs$R2R$E)

    dRes <- eRNAkit::subDB(emi, "E", eRs)
    dRes$R2R <- dRes$R2R[dRes$R2R$G %in% eRNAkit::stringS(input$byGene), ]
    dRes[["subloc"]] <- eRNAkit::winLoc(dRes$loc)

    filtered_db(dRes)
    filter_status("DB filtering done! You can now proceed.")

  })


  ## Give status on filtering here
  output$filterStatusMessage <- renderText({
    filter_status()
  })


  # Event Logic
  observeEvent(input$loc,  { selected_view("loc") })
  observeEvent(input$r2r,  { selected_view("r2r") })
  observeEvent(input$ee,   { selected_view("ee") })
  observeEvent(input$pacy,   { selected_view("pacy") })
  observeEvent(input$ko,   { selected_view("ko") })


  ## Dynamic layout
  output$mainPanelUI <- renderUI({
    view <- selected_view()
    req(view)

    if (view == "loc") {
      tagList(fluidRow(
          column(width = 8, plotOutput("plot1")),
          column(width = 4, plotOutput("plot2"))),
        fluidRow(column(width = 8, DTOutput("table1")),
                 column(width = 4, plotOutput("plot3"))) )

    } else if (view == "r2r") {
      tagList(fluidRow(
        column(width = 6, plotOutput("plot4")),
        column(width = 6, plotOutput("plot5"))),
        fluidRow(column(width = 12, DTOutput("table2"))))
    } else if (view == "ee") {
      tagList(fluidRow(
          column(width = 6, plotOutput("plot6")),
          column(width = 6, plotOutput("plot7"))),
        fluidRow(column(width = 12, DTOutput("table3"))))

    } else if (view == "pacy") {
      tagList(fluidRow(column(width = 12, DTOutput("table4"))))
    } else if (view == "ko") {
      tagList(fluidRow(column(width = 12, DTOutput("table5"))))
    }
  })

## LOCALISATION View
  output$plot1 <- renderPlot({ req(selected_view() == "loc")
    req(filtered_db()$subloc$p1)
    eRNAkit::plotLoc(filtered_db()$subloc$p1, "p1") }, res = 120)

  output$table1 <- renderDT({ req(selected_view() == "loc")
    req(filtered_db()$subloc$d)
    datatable(filtered_db()$subloc$d) })

  output$plot2 <- renderPlot({ req(selected_view() == "loc")
    req(filtered_db()$subloc$p2)
    eRNAkit::plotLoc(filtered_db()$subloc$p2, "p2") }, res = 120)


  output$plot3 <- renderPlot({ req(selected_view() == "loc")
    req(filtered_db()$subloc$p3)
    eRNAkit::plotLoc(filtered_db()$subloc$p3, "p3") }, res=120)




## TARGET View
  output$plot4 <- renderPlot({ req(selected_view() == "r2r")
    req(filtered_db()$R2R)
    df <- filtered_db()$R2R
    #
    graph <- graph_from_data_frame(df[, c("E", "G")],
                                   directed = FALSE)
    E(graph)$color <- "black"

    l <- layout_with_fr(graph)
    l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

    #
    plot(graph,
         layout = l * 1,
         rescale=F,
         vertex.color = "darkred",
         vertex.label.cex = 0.4,
         vertex.label.dist = 2,
         edge.width = 0.2,
         vertex.label.position = 1,
         margin = 0,
         main = "eRNA to mRNA Targets") }, res=120)

  output$plot5 <- renderPlot({ req(selected_view() == "r2r")
    req(filtered_db()$R2R)
    ggplot(filtered_db()$R2R, aes(x = reorder(pair, -count), y = count)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(hjust = 1, angle = 90,  size = 8)) +
      labs(title = "Observed Interactions Counts", x = "", y = "Counts") }, res=120)

  output$table2 <- renderDT({ req(selected_view() == "r2r")
    df <- filtered_db()$R2R
    req(df)
    datatable(df[, c("pair", "Gs", "biotype", "method", "D2D",
                     "chr.EG", "RSV", "VSV", "cycloheximide", "harringtonine",
                     "DDX3X", "arsenite", "SRSF1")]) })


## EXPRESSION View
  output$plot6 <- renderPlot({ req(selected_view() == "ee")
    df <- filtered_db()$EOC
    df <- df[df$type %in% "organ", ]
    req(df)
    eRNAkit::plotOC(df, "Organ")
    }, res=120)

  output$plot7 <- renderPlot({ req(selected_view() == "ee")
    df <- filtered_db()$EOC
    df <- df[df$type %in% "cell", ]
    req(df)
    eRNAkit::plotOC(df, "Cells")
    }, res=120)

  output$table3 <- renderDT({ req(selected_view() == "ee")
    req(filtered_db()$EOC)

    datatable(filtered_db()$EOC)
    })

  output$table4 <- renderDT({ req(selected_view() == "pacy")
    df <- filtered_db()$PAcy
    df2 <- filtered_db()$core
    req(df)
    req(df2)
    df2 <- df2[c("E", "source")]
    datatable(merge(df, df2, by="E"))
  })


  # DEFINE download logic -------------------------------------------------
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("filtered-results-", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      saveDB <- filtered_db()
      saveDB$subloc <- NULL
      req(saveDB)

      saveDB <- lapply(saveDB, function(x) {
        if (is.null(x)) data.frame() else x
      })

      writexl::write_xlsx(saveDB, path = file)
    }
  )


}




