#---------------testing create_gating_df()-----------------#

# create an empty gating data frame using create_gating_df()
test_df <- create_gating_df()

test_that("create_gating_df creates properly formatted empty dataframe", {
    expect_equal(dim(test_df), c(0, 18))
    expect_type(test_df$Gate_ID, "character")
    expect_type(test_df$Assay_Name, "character")
    expect_type(test_df$X_Axis, "character")
    expect_type(test_df$Y_Axis, "character")
    expect_type(test_df$Subset_Name, "character")
    expect_type(test_df$Num_Input_Cells, "integer")
    expect_type(test_df$Num_Subset_Cells, "integer")
    expect_type(test_df$Total_Cells_in_Sample, "integer")
    expect_type(test_df$Percent_Subsetted_From_Previous, "double")
    expect_type(test_df$Percent_Subsetted_From_Total, "double")
    expect_type(test_df$Input_Cells, "list")
    expect_type(test_df$Input_X_Coordinates, "list")
    expect_type(test_df$Input_Y_Coordinates, "list")
    expect_type(test_df$Subset_Cells, "list")
    expect_type(test_df$Subset_X_Coordinates, "list")
    expect_type(test_df$Subset_Y_Coordinates, "list")
    expect_type(test_df$Gate_X_Coordinates, "list")
    expect_type(test_df$Gate_Y_Coordinates, "list")
})

#---------------testing Gate()-----------------#

# create Gate objects for testing
test_gate_1 <- Gate(
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


test_that("Gate() create correctly formatted Gate object", {
    expect_type(test_gate_1@assay_name, "character")
    expect_type(test_gate_1@num_input_cells, "integer")
    expect_type(test_gate_1@input_cells, "list")
    expect_equal(test_gate_1@x_axis, "ADT-A")
    expect_equal(test_gate_1@y_axis, "ADT-B")
})

#---------------testing update_gating_df()-----------------#

# update_gating_df() requires a list of gate objects to choose from

test_gate_2 <- Gate(
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

test_gate_list <- list("gate1" = test_gate_1, "gate2" = test_gate_2)

test_df <- update_gating_df(gate_name_string = "gate1", reactive_gate_list = test_gate_list, temp_gating_df = test_df)
test_df <- update_gating_df(gate_name_string = "gate2", reactive_gate_list = test_gate_list, temp_gating_df = test_df)

test_that("update_gating_df() updated gating dataframes correctly", {
    expect_equal(test_df$Assay_Name, c("test1", "test2"))
    expect_equal(test_df$Num_Input_Cells[2], test_gate_2@num_input_cells)
    expect_equal(test_df$Num_Input_Cells[1], test_gate_1@num_input_cells)
    expect_equal(test_df$Input_X_Coordinates[[1]], test_gate_1@input_coords$x)
})

#---------------testing SetSubsetName()-----------------#

test_gate_3 <- SetSubsetName(test_gate_1, "test_cells_C")

test_that("Changing a gate's label", {
    expect_equal(test_gate_3@name_subset_cells, "test_cells_C")
})