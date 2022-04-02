#create gating datatable for display in UI
#' Title
#'
#' @param temp_reactive_gating_df 
#'
#' @return
#' @export
#'
#' @examples
create_gating_dt <- function(temp_reactive_gating_df) {
    temp_gating_dt <- DT::datatable(temp_reactive_gating_df, 
                                    rownames = TRUE, 
                                    selection = "single", #only let user click 1 row at a time so only that one gate can be shown in plot. Otherwise, app will crash.
                                    editable = list(target = "cell", disable = list(columns = c(0:4, 6:ncol(temp_reactive_gating_df)))), # only cells in column 5(subset_name) can be edited by user
                                    extensions = c("Buttons", "Scroller"),
                                    options = list(deferRender = TRUE,
                                                scroller = TRUE,
                                                scrollY = 300,
                                                scrollX = TRUE,
                                                dom = "frtipB", 
                                                buttons = c("copy", "print"),
                                                columnDefs = list(list(visible = FALSE, targets = c(11:ncol(temp_reactive_gating_df))))
                                    )) %>%
    formatRound(c("Percent_Subsetted_From_Previous", "Percent_Subsetted_From_Total"), digits = 4)
    
    return(temp_gating_dt)
}