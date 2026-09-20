# ---------------------------------------------------------------------------
# R Shiny graphical user interface
# ---------------------------------------------------------------------------

nanoamp_gui_default_outdir <- function() {
  file.path(path.expand("~"), "nanoamp_results")
}

nanoamp_gui_has_dt <- function() {
  requireNamespace("DT", quietly = TRUE)
}

nanoamp_gui_table_output <- function(id) {
  if (nanoamp_gui_has_dt()) DT::DTOutput(id) else shiny::tableOutput(id)
}

nanoamp_gui_render_table <- function(expr) {
  if (nanoamp_gui_has_dt()) {
    DT::renderDT(expr, options = list(pageLength = 10, scrollX = TRUE))
  } else {
    shiny::renderTable(expr)
  }
}

nanoamp_gui_drop_sequence <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df <- as.data.frame(df)
  df[, setdiff(names(df), "sequence"), drop = FALSE]
}

nanoamp_gui_open_dir <- function(path) {
  if (is.null(path) || !dir.exists(path)) return(invisible(FALSE))
  if (.Platform$OS.type == "windows") {
    try(shell.exec(path), silent = TRUE)
  } else if (identical(Sys.info()[["sysname"]], "Darwin")) {
    try(system2("open", path), silent = TRUE)
  } else {
    try(system2("xdg-open", path), silent = TRUE)
  }
  invisible(TRUE)
}

nanoamp_gui_check <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required for the GUI. Install it with install.packages('shiny').",
         call. = FALSE)
  }
  invisible(TRUE)
}

nanoamp_gui_ui <- function() {
  shiny::fluidPage(
    shiny::titlePanel("nanoamp - Nanopore Amplicon Haplotype Analysis"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 4,
        shiny::fileInput("reads", "FASTQ file", accept = c(".fastq", ".fq", ".gz")),
        shiny::fileInput("reference", "Reference FASTA",
                         accept = c(".fa", ".fasta", ".fna", ".seq")),
        shiny::textInput("outdir", "Output directory",
                         value = nanoamp_gui_default_outdir()),
        shiny::radioButtons(
          "mode", "Analysis mode",
          choices = c(
            "A - reference-guided (recommended)" = "A",
            "B - de novo clustering" = "B",
            "C - exact match (diagnostic)" = "C"
          ),
          selected = "A"
        ),
        shiny::numericInput("top_n", "Top N haplotypes", value = 20, min = 1, max = 1000),
        shiny::hr(),
        shiny::strong("Advanced parameters"),
        shiny::numericInput("min_reads", "min_reads", value = 3, min = 1),
        shiny::numericInput("min_freq", "min_freq", value = 0.02, min = 0, max = 1, step = 0.01),
        shiny::numericInput("min_identity", "min_identity", value = 0.90, min = 0, max = 1, step = 0.01),
        shiny::numericInput("identity_cutoff", "identity_cutoff (Mode B)",
                            value = 0.99, min = 0.5, max = 1, step = 0.001),
        shiny::numericInput("min_cluster_reads", "min_cluster_reads (Mode B)",
                            value = 2, min = 1),
        shiny::selectInput("consensus_method", "consensus_method (Mode B)",
                           choices = c("decipher", "medoid"), selected = "decipher"),
        shiny::selectInput("aligner", "Alignment backend",
                           choices = c("minimap2", "r"), selected = "minimap2"),
        shiny::numericInput("threads", "threads", value = 4, min = 1, max = 64),
        shiny::checkboxInput("keep_intermediates", "Keep BAM and intermediate files", value = TRUE),
        shiny::hr(),
        shiny::actionButton("run", "Run analysis", class = "btn-primary"),
        shiny::br(), shiny::br(),
        shiny::actionButton("open_dir", "Open output folder"),
        shiny::hr(),
        shiny::downloadButton("download_haplotypes", "Download haplotypes.tsv"),
        shiny::downloadButton("download_variants", "Download variants.tsv")
      ),
      shiny::mainPanel(
        width = 8,
        shiny::verbatimTextOutput("status"),
        shiny::tabsetPanel(
          shiny::tabPanel("Log", shiny::verbatimTextOutput("log")),
          shiny::tabPanel("Haplotypes", nanoamp_gui_table_output("haplotypes")),
          shiny::tabPanel("Variants", nanoamp_gui_table_output("variants")),
          shiny::tabPanel("QC", shiny::tableOutput("qc"))
        )
      )
    )
  )
}

nanoamp_gui_server <- function(input, output, session) {
  values <- shiny::reactiveValues(
    result = NULL, log = "Ready.", outdir = NULL, status = "idle"
  )

  shiny::observeEvent(input$run, {
    shiny::req(input$reads, input$reference)
    outdir <- path.expand(input$outdir)
    if (!nzchar(outdir)) {
      values$log <- "Please provide an output directory."
      return()
    }
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    values$outdir <- outdir
    values$status <- "running"
    values$log <- "Running analysis..."

    log_text <- ""
    res <- NULL
    shiny::withProgress(message = "Running nanoamp analysis", value = 0, {
      shiny::incProgress(0.2, detail = "Aligning reads and calling variants")
      res <- tryCatch({
        log_text <- paste(utils::capture.output({
          res <- run_haplotype_analysis(
            reads = input$reads$datapath,
            reference = input$reference$datapath,
            outdir = outdir,
            mode = input$mode,
            top_n = as.integer(input$top_n),
            min_reads = as.integer(input$min_reads),
            min_freq = as.numeric(input$min_freq),
            min_identity = as.numeric(input$min_identity),
            identity_cutoff = as.numeric(input$identity_cutoff),
            min_cluster_reads = as.integer(input$min_cluster_reads),
            consensus_method = input$consensus_method,
            aligner = input$aligner,
            threads = as.integer(input$threads),
            keep_intermediates = isTRUE(input$keep_intermediates)
          )
        }, type = "output"), collapse = "\n")
        res
      }, error = function(e) {
        log_text <<- paste(log_text, conditionMessage(e), sep = "\n")
        NULL
      })
      shiny::incProgress(0.7, detail = "Finalizing outputs")
    })

    values$result <- res
    if (nzchar(log_text)) values$log <- log_text
    values$status <- if (is.null(res)) "error" else "done"
  })

  output$status <- shiny::renderText({
    if (identical(values$status, "done")) {
      paste0("Analysis complete. Output: ", values$outdir)
    } else if (identical(values$status, "error")) {
      "Analysis failed. See the Log tab."
    } else if (identical(values$status, "running")) {
      "Analysis running..."
    } else {
      "Ready."
    }
  })

  output$log <- shiny::renderText(values$log)

  output$haplotypes <- nanoamp_gui_render_table({
    shiny::req(values$result)
    nanoamp_gui_drop_sequence(values$result$haplotypes)
  })

  output$variants <- nanoamp_gui_render_table({
    shiny::req(values$result)
    values$result$variants
  })

  output$qc <- shiny::renderTable({
    shiny::req(values$result)
    qc <- values$result$qc
    data.frame(metric = names(qc), value = unlist(qc, use.names = FALSE))
  })

  output$download_haplotypes <- shiny::downloadHandler(
    filename = function() "haplotypes.tsv",
    content = function(file) {
      shiny::req(values$outdir)
      src <- file.path(values$outdir, "haplotypes.tsv")
      if (file.exists(src)) file.copy(src, file, overwrite = TRUE)
      else writeLines("No results yet.", file)
    }
  )

  output$download_variants <- shiny::downloadHandler(
    filename = function() "variants.tsv",
    content = function(file) {
      shiny::req(values$outdir)
      src <- file.path(values$outdir, "variants.tsv")
      if (file.exists(src)) file.copy(src, file, overwrite = TRUE)
      else writeLines("No results yet.", file)
    }
  )

  shiny::observeEvent(input$open_dir, {
    nanoamp_gui_open_dir(values$outdir)
  })
}

#' Create the nanoamp Shiny app object
#'
#' Build the Shiny application object used by the GUI. This is useful for
#' testing and for embedding the app in other launchers.
#'
#' @return A `shiny.appobj`.
#' @export
nanoamp_gui_app <- function() {
  nanoamp_gui_check()
  shiny::shinyApp(ui = nanoamp_gui_ui(), server = nanoamp_gui_server)
}

#' Launch the nanoamp graphical user interface
#'
#' Start a local Shiny server and open the GUI in a web browser.
#'
#' @param launch.browser Whether to open the app in the default browser.
#' @param port Port to listen on; `NULL` lets Shiny choose a port.
#' @param host Host address to bind.
#'
#' @return Invisibly, the value returned by [shiny::runApp()].
#' @export
nanoamp_gui <- function(launch.browser = TRUE, port = NULL, host = "127.0.0.1") {
  nanoamp_gui_check()
  shiny::runApp(
    nanoamp_gui_app(),
    launch.browser = launch.browser,
    port = port,
    host = host
  )
}
