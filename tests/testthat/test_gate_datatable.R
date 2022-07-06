#---------------testing create_gating_dt()-----------------#

# create an empty gating data frame using create_gating_df()
test_df <- create_gating_df()

# create Gate objects for testing
test_gate1 <- Gate(
    counter = as.integer(1),
    assay_name = "test1",
    input_cells = list(c("a", "b", "c")),
    input_coords = data.frame(x = c(1, 2, 3), y = c(4, 5, 6)),
    subset_cells = list(c("a", "b")),
    subset_coords = data.frame(x = c(1, 2), y = c(4, 5)),
    x_axis = "ADT-A",
    y_axis = "ADT-B",
    gate_coords = list(x = c(1, 2, 3, 4), y = c(5, 6, 7, 8)),
    name_subset_cells = "test_cells_A",
    num_input_cells = as.integer(1000),
    num_subset_cells = as.integer(500),
    total_num_cells_in_sample = as.integer(1000),
    pct_subset_from_previous = 50,
    pct_subset_from_total = 50
)

# update_gating_df() requires a list of gate objects to choose from

test_gate2 <- Gate(
    counter = as.integer(2),
    assay_name = "test2",
    input_cells = list(c("a", "b")),
    input_coords = data.frame(x = c(1, 2), y = c(4, 5)),
    subset_cells = list(c("a")),
    subset_coords = data.frame(x = c(1), y = c(4)),
    x_axis = "ADT-C",
    y_axis = "ADT-D",
    gate_coords = list(x = c(10, 20, 30, 40), y = c(50, 60, 70, 80)),
    name_subset_cells = "test_cells_B",
    num_input_cells = as.integer(500),
    num_subset_cells = as.integer(100),
    total_num_cells_in_sample = as.integer(1000),
    pct_subset_from_previous = 20,
    pct_subset_from_total = 10
)

test_gate_list <- list("gate1" = test_gate1, "gate2" = test_gate2)

test_df <- update_gating_df(gate_name_string = "gate1", reactive_gate_list = test_gate_list, temp_gating_df = test_df)
test_df <- update_gating_df(gate_name_string = "gate2", reactive_gate_list = test_gate_list, temp_gating_df = test_df)

test_dt <- create_gating_dt(test_df)

test_that("create_gating_dt() ran correctly", {
    expect_visible(test_dt)
    expect_type(object = test_dt, type = "list")
    expect_silent(test_dt)
    expect_error(create_gating_dt(test_gate1))
})
