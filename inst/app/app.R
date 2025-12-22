library(shiny)
library(DT)
library(shinycssloaders)
library(SummarizedExperiment)
library(msdialSE)

ui <- fluidPage(
  titlePanel("MS-DIAL alignment file to SummarizedExperiment (lipidomics)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "MS-DIAL alignment CSV", accept = ".csv"),
      numericInput("ann", "Annotation columns (N)", value = 35, min = 1, step = 1),
      numericInput("hdr", "Header rows", value = 5, min = 1, step = 1),
      numericInput("start", "Data start row", value = 6, min = 1, step = 1),
      hr(),
      # In UI sidebarPanel(...)
      selectInput("norm", "Normalization", choices = c("none","log2","sum","median","pqn","quantile","qc_loess")),
      conditionalPanel(
        condition = "input.norm == 'log2'",
        numericInput("norm_offset", "log2 offset", value = 1e-9, min = 0, step = 1e-9)
      ),
      conditionalPanel(
        condition = "input.norm == 'qc_loess'",
        textInput("qc_label", "QC label in colData$class", value = "QC"),
        textInput("order_col", "Injection order column (colData)", value = "order"),
        numericInput("loess_span", "LOESS span", value = 0.75, min = 0.2, max = 1, step = 0.05)
      ),
      actionButton("build", "Build SE"),
      hr(),
      downloadButton("dl_rds", "Download SE (.rds)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary", verbatimTextOutput("summary") %>% withSpinner()),
        tabPanel("rowData (head)", DTOutput("rd") %>% withSpinner()),
        tabPanel("colData", DTOutput("cd") %>% withSpinner())
      )
    )
  )
)

server <- function(input, output, session) {
  se_react <- eventReactive(input$build, {
    req(input$file)
    se <- msdialSE::load_lipidomics_se(
      csv_path       = input$file$datapath,
      annotation_cols= seq_len(input$ann),
      header_rows    = input$hdr,
      data_start_row = input$start
    )
    # In server inside eventReactive after load_lipidomics_se(...)
    norm_ctl <- switch(input$norm,
                       log2     = list(offset = input$norm_offset),
                       qc_loess = list(qc_label = input$qc_label, order_col = input$order_col, span = input$loess_span),
                       list()
    )
    se <- msdialSE::normalize_se(se, method = input$norm, control = norm_ctl)
    se
  }, ignoreInit = TRUE)
  
  output$summary <- renderPrint({
    se <- se_react()
    if (is.null(se)) return("No SE built yet.")
    list(
      dimensions = dim(SummarizedExperiment::assay(se, "abundance")),
      assay_names = SummarizedExperiment::assayNames(se),
      rowData_cols = colnames(SummarizedExperiment::rowData(se)),
      colData_cols = colnames(SummarizedExperiment::colData(se)),
      subclass_example = head(SummarizedExperiment::rowData(se)$subclass)
    )
  })
  
  output$rd <- renderDT({
    se <- se_react(); req(se)
    as.data.frame(SummarizedExperiment::rowData(se)) |>
      head(50) |>
      datatable(options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$cd <- renderDT({
    se <- se_react(); req(se)
    as.data.frame(SummarizedExperiment::colData(se)) |>
      datatable(options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$dl_rds <- downloadHandler(
    filename = function() sprintf("lipidomics_se_%s.rds", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      se <- se_react(); req(se)
      saveRDS(se, file)
    }
  )
}

shinyApp(ui, server)
