#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
#'
#' @return Initialization of CITEViz app
#' @examples
#'
#' run_app()
#'
run_app <- function(on_start = NULL,
    options = list(),
    enable_bookmarking = NULL,
    ui_pattern = "/",
    ...) {
    with_golem_options(
        app = shinyApp(
            ui = app_ui,
            server = app_server,
            # on_start = on_start,
            options = options,
            # enable_bookmarking = enable_bookmarking,
            # ui_pattern = ui_pattern
        ),
        golem_opts = list(...)
    )
}
