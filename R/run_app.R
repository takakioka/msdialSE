#' Run the bundled Shiny app
#' @export
run_app <- function() {
  app_dir <- system.file("app", package = "msdialSE")
  if (app_dir == "") stop("App not found. Try re-installing msdialSE.", call. = FALSE)
  shiny::runApp(app_dir, display.mode = "normal")
}
