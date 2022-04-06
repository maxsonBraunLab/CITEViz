#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' 
#' @import shiny
#' 
#' @noRd
app_server <- function( input, output, session ) {
  
  # ---------- Packages list ---------- 
  # Change the package library list below to fit Bioconductr guidelines!!! 
  # Specify necessary packages in NAMESPACE and DESCRIPTION files instead!
  
  library(bslib)    #needed for custom bootstrap theme
  library(DT)       #needed to create HTML datatable widgets
  library(ggplot2)
  library(plotly)
  library(RColorBrewer)
  library(Seurat)
  library(shiny)
  library(shinyBS)  #bootstrap components for shiny
  library(shinyjs)
  library(stats)
  library(tidyverse)
  library(vembedr)  #needed to embed Youtube videos in landing/welcome tab
  
  
  #set max file upload size to 3gb (default is only 5mb) since rds files can be really big
  #takes in an integer argument for max filesize in megabytes
  #to improve: set max file upload size based on user's hardware limitations?
  options(shiny.maxRequestSize = 3000 * 1024^2)
  
  # Your application server logic 
  # ---------- Class definitions ----------
  
  # define a Gate class for gate objects
  #' Title
  #'
  #' @slot counter integer. 
  #' @slot assay_name character. 
  #' @slot input_cells list. 
  #' @slot input_coords data.frame. 
  #' @slot subset_cells list. 
  #' @slot subset_coords data.frame. 
  #' @slot x_axis character. 
  #' @slot y_axis character. 
  #' @slot gate_coords list. 
  #' @slot name_subset_cells character. 
  #' @slot num_input_cells integer. 
  #' @slot num_subset_cells integer. 
  #' @slot total_num_cells_in_sample integer. 
  #' @slot pct_subset_from_previous numeric. 
  #' @slot pct_subset_from_total numeric. 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  setClass("Gate", slots = list(
    counter = "integer", #to keep track of gate number in the case of multiple gating steps
    assay_name = "character", #to keep track of which assay (RNA, ADT, etc) the gate data is coming from
    
    ## cells and coordinates, seurat data
    input_cells = "list",
    input_coords = "data.frame",
    subset_cells = "list",
    subset_coords = "data.frame",
    
    ## gate information
    x_axis = "character",
    y_axis = "character",
    gate_coords = "list",
    
    ## cell summary data
    name_subset_cells = "character", 
    num_input_cells = "integer",
    num_subset_cells = "integer",
    total_num_cells_in_sample = "integer",
    pct_subset_from_previous = "numeric", #percent of cells in the current gate that were subsetted from previous gate (100 * num_subset_cells/num_input_cells)
    pct_subset_from_total = "numeric" #percent of cells in the current gate that were subsetted from the original total number of cells in the sample
  ))
  
  
  #constructor for Gate objects
  #' Title
  #'
  #' @param counter 
  #' @param assay_name 
  #' @param input_cells 
  #' @param input_coords 
  #' @param subset_cells 
  #' @param subset_coords 
  #' @param x_axis 
  #' @param y_axis 
  #' @param gate_coords 
  #' @param name_subset_cells 
  #' @param num_input_cells 
  #' @param num_subset_cells 
  #' @param total_num_cells_in_sample 
  #' @param pct_subset_from_previous 
  #' @param pct_subset_from_total 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  Gate <- function(counter=NA_integer_, assay_name = NA_character_, input_cells = list(), input_coords = data.frame(), 
                   subset_cells = list(), subset_coords = data.frame(), x_axis = NA_character_, y_axis = NA_character_,
                   gate_coords = list(), name_subset_cells = NA_character_, num_input_cells = NA_integer_, num_subset_cells = NA_integer_,
                   total_num_cells_in_sample = NA_integer_, pct_subset_from_previous = NA_real_, pct_subset_from_total = NA_real_) {
    new("Gate", 
        counter=counter, 
        assay_name = assay_name, 
        input_cells = input_cells, 
        input_coords = input_coords, 
        subset_cells = subset_cells, 
        subset_coords = subset_coords, 
        x_axis = x_axis, 
        y_axis = y_axis,
        gate_coords = gate_coords, 
        name_subset_cells = name_subset_cells, 
        num_input_cells = num_input_cells, 
        num_subset_cells = num_subset_cells,
        total_num_cells_in_sample = total_num_cells_in_sample, 
        pct_subset_from_previous = pct_subset_from_previous, 
        pct_subset_from_total = pct_subset_from_total)
  }
  
  
  # ---------- Function definitions ----------
  
  #create custom color palette with color interpolation to use for QA plots
  #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
  #takes in an integer argument for number of colors to generate in palette
  #returns a character vector of color hex codes for the desired number of colors
  #to improve: add a colorblind-friendly palette option
  #' Title
  #'
  #' @param num_colors 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  get_palette <- function(num_colors = integer(3)) {
    base_palette <- colorRampPalette(colors = c("turquoise3",
                                                "darkgreen",
                                                "red",
                                                "royalblue1",
                                                "orange2",
                                                "lightseagreen",
                                                "yellow",
                                                "royalblue4",
                                                "yellow3",
                                                "darkorchid1",
                                                "lawngreen",
                                                "magenta3"))
    return(base_palette(num_colors))
  }
  
  
#create gate objects from user input selections depending on if the user is forward-gating or back-gating
#' Title
#'
#' @param is_forward_gating 
#' @param assay_count_data 
#' @param gate_counter 
#' @param reactive_gate_list 
#' @param reactive_selected_gate 
#' @param reactive_last_buttons_clicked 
#'
#' @return
#' @export
#'
#' @examples
  create_gate_from_input <- function(is_forward_gating = TRUE, assay_count_data, gate_counter, reactive_gate_list, reactive_selected_gate, reactive_last_buttons_clicked) {
    # "ui_input_suffix" refers to the suffix that is at the end of the UI input elements for forward and backgating
    # So for forward-gating, an example of an input from a UI element would be input$Assay. The corresponding input in the backgating page would be input$Assay_bg. the suffix for the backgating UI elements is "_bg"
    ui_input_suffix <- ""
    
    sel <- NA_character_
    brushed_coords <- NA_character_
    input_cells <- list()
    input_coords <- data.frame()

    # forward-gating logic
    if (is_forward_gating == TRUE) {
      # get plotly event data
      sel <- event_data("plotly_selected", source = "C")
      brushed_coords <- event_data("plotly_brushed", source = "C")

      # Check for Gate object. If one does not exists, create initial input for first gate.
      if (length(reactive_gate_list) == 0 | reactive_last_buttons_clicked$second_to_last == "reset_button" | reactive_last_buttons_clicked$second_to_last == "clear_all_gates_button") {
        input_cells <- list(rownames(assay_count_data))
        input_coords <- data.frame(x = assay_count_data[,1], y = assay_count_data[,2])
      }
      else {
        #set up inputs to gate object if there is already a gate object
        # if user has not selected a previous gate from datatable to re-gate from
        if (is.null(input$gating_pg_table_rows_selected)) {
          input_cells <- reactive_gate_list[[paste0("gate_", gate_counter - 1)]]@subset_cells
          input_coords <- data.frame(x = brushed_coords$x, y = brushed_coords$y)
        }
        # else if user has selected a previous gate from datatable that they want to re-gate from
        else {
          input_cells <- reactive_gate_list[[reactive_selected_gate]]@subset_cells
          input_coords <- reactive_gate_list[[reactive_selected_gate]]@subset_coords
        }
      }
    }
    # back-gating logic
    else {
      ui_input_suffix <- "_bg"

      # get plotly event data and input cell data
      sel <- event_data("plotly_selected", source = "D")
      brushed_coords <- event_data("plotly_brushed", source = "D")
      input_cells <- list(rownames(count_data))
      input_coords <- data.frame(x = count_data[,1], y = count_data[,2])
    }

    #assign values to variables
    assay_name <- input[[paste0("Assay", ui_input_suffix)]]
    subset_cells <- list(sel$customdata)
    subset_coords <- data.frame(x = assay_count_data[sel$customdata,1], y = assay_count_data[sel$customdata,2])
    x_axis <- input[[paste0("x_feature", ui_input_suffix)]]
    y_axis <- input[[paste0("y_feature", ui_input_suffix)]]
    gate_coords <- list(x = c(brushed_coords$x), y = c(brushed_coords$y))
    name_subset_cells <- input[[paste0("cell_subset_name", ui_input_suffix)]]
    num_input_cells <- length(unlist(input_cells))
    num_subset_cells <- length(unlist(subset_cells))
    total_num_cells_in_sample <- length(rownames(assay_count_data))
    pct_subset_from_previous <- 100 * num_subset_cells / num_input_cells
    pct_subset_from_total <- 100 * num_subset_cells / total_num_cells_in_sample

    # construct Gate object
    gate_obj <- Gate(counter = gate_counter,
                    assay_name = assay_name,
                    input_cells = input_cells,
                    input_coords = input_coords,
                    subset_cells = subset_cells,
                    subset_coords = subset_coords,
                    x_axis = x_axis,
                    y_axis = y_axis,
                    gate_coords = gate_coords,
                    name_subset_cells = name_subset_cells,
                    num_input_cells = num_input_cells,
                    num_subset_cells = num_subset_cells,
                    total_num_cells_in_sample = total_num_cells_in_sample,
                    pct_subset_from_previous = pct_subset_from_previous,
                    pct_subset_from_total = pct_subset_from_total)
    return(gate_obj)
  }
  
  
  #' Title
  #'
  #' @param gating_reactiveValues 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  get_reactive_gate_list <- function(gating_reactiveValues) {
    reactive_gate_list <- reactive({
      unordered_list <- reactiveValuesToList(gating_reactiveValues)
      ordered_list <- unordered_list[order(names(unordered_list))]
      ordered_list
    })
    return(reactive_gate_list)
  }
  
  #create gating dataframe
  #' Title
  #'
  #' @return
  #' @export
  #'
  #' @examples
  create_gating_df <- function() {
    data.frame(
      Gate_ID = character(),
      Assay_Name = character(),
      X_Axis = character(),
      Y_Axis = character(),
      Subset_Name = character(),
      Num_Input_Cells = integer(),
      Num_Subset_Cells = integer(),
      Total_Cells_in_Sample = integer(),
      Percent_Subsetted_From_Previous = numeric(),
      Percent_Subsetted_From_Total = numeric(),
      
      Input_Cells = I(list()),
      Input_X_Coordinates = I(list()),
      Input_Y_Coordinates = I(list()),
      Subset_Cells = I(list()),
      Subset_X_Coordinates = I(list()),
      Subset_Y_Coordinates = I(list()),
      Gate_X_Coordinates = I(list()),
      Gate_Y_Coordinates = I(list())
    )
  }
  
  #function to populate/update gating data frame
  #' Title
  #'
  #' @param gate_name_string 
  #' @param reactive_gate_list 
  #' @param temp_gating_df 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  update_gating_df <- function(gate_name_string, reactive_gate_list, temp_gating_df) {
    gate_obj <- reactive_gate_list[[gate_name_string]]
    if (!is.null(gate_obj)) {
      temp_gating_df <- add_row(temp_gating_df,
                                Gate_ID = gate_name_string,
                                Assay_Name = gate_obj@assay_name,
                                X_Axis = gate_obj@x_axis,
                                Y_Axis = gate_obj@y_axis,
                                Subset_Name = gate_obj@name_subset_cells,
                                Num_Input_Cells = gate_obj@num_input_cells,
                                Num_Subset_Cells = gate_obj@num_subset_cells,
                                
                                Total_Cells_in_Sample = gate_obj@total_num_cells_in_sample,
                                Percent_Subsetted_From_Previous = gate_obj@pct_subset_from_previous,
                                Percent_Subsetted_From_Total = gate_obj@pct_subset_from_total,
                                
                                Input_Cells = gate_obj@input_cells,
                                Input_X_Coordinates = list(gate_obj@input_coords$x),
                                Input_Y_Coordinates = list(gate_obj@input_coords$y),
                                Subset_Cells = gate_obj@subset_cells,
                                Subset_X_Coordinates = list(gate_obj@subset_coords$x),
                                Subset_Y_Coordinates = list(gate_obj@subset_coords$y),
                                Gate_X_Coordinates = list(gate_obj@gate_coords$x),
                                Gate_Y_Coordinates = list(gate_obj@gate_coords$y)
      )
    }
    return(temp_gating_df)
  }
  
  
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
  
  
  # reset all values in gate_reactive_values to NULL
  #' Title
  #'
  #' @param gate_name_string 
  #' @param local_gate_reactive_values 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  set_gates_to_null <- function(gate_name_string, local_gate_reactive_values) {
    local_gate_reactive_values[[gate_name_string]] <- NULL
    return(local_gate_reactive_values[[gate_name_string]])
  }
  
  
  # ---------- Filehandling ---------- 
  # Putting everything in an observe function will put everything in the server function into the same environment allowing for
  # a single read of the uploaded seurat object instead of a read everytime myso is called in render* function.
  # This speeds up the code immensely
  # The main function of this initial observe is to allow for a single upload of a Seurat Object over all pages.
  # Filetype validation: 
  # Note that only RDS files can be inputted by the user due to the UI fileInput() argument `accept = ".rds"`. 
  # This works on an actual web browser but not in the RStudio viewer.
  
  observe({
    
    inp_file <- input$rds_input_file
    if (is.null(inp_file)) {
      return(NULL)
    }
    
    #read in RDS file
    rds_obj <- readRDS(inp_file$datapath)
    
    if (typeof(rds_obj) == "list") {
      #check if integrated obj exists before retrieving it from RDS that was read in
      integrated_obj_index <- grep("integrated", names(rds_obj), ignore.case = TRUE)
      
      #if integrated obj is not in list of objs and each seurat obj in the sample list is also wrapped in a list
      if (length(integrated_obj_index) == 0) {
        myso <- rds_obj[[1]]
      }
      else {
        #read in integrated obj
        myso <- rds_obj[[integrated_obj_index]]
      }
    }
    else {
      #if integrated obj is not in list of objs, and the obj from the RDS file is just a single Seurat obj and not a list
      myso <- rds_obj
    }
    
    #check if Seurat object has "reductions" slot
    reduction_validation <- reactive({ 
      validate(
        need(
          length(myso@reductions) > 0,
          message = "Seurat object does not contain reductions. Please check input RDS file."
        )
      )
    })
    
    output$reduction_validation_status <- renderText({ reduction_validation() })
    
    
    # ---------- ***** QA ***** ---------- 
    # QA plots generated from Maxson-Braun lab's CITE-seq data preprocessing pipeline
    
    observe({
      
      updateSelectInput(
        session = session,
        inputId = "color_qa",
        choices = colnames(myso@meta.data[lapply(myso@meta.data, class) %in% c("factor", "character")]),
        selected = colnames(myso@meta.data[lapply(myso@meta.data, class) %in% c("factor", "character")])[1]
      )
      
      # ----- QA distribution plot -----
      # reactive distribution plot
      #reactive function will rerun this expression every time distribution_plot is called, which should be only when QA or color_qa choice is changed
      distribution_plot <- reactive({
        
        req(input$rds_input_file, input$QA, input$color_qa)
        color <- input$color_qa
        
        #This assigns the params variable a list of strings that act a varying parameters depending on input of QA
        params <- switch(input$QA,
                         "RNA Count Per Cell" = c("nCount_RNA", "Distribution of Counts per Cell", "Number of Counts"),
                         "Gene Count Per Cell" = c("nFeature_RNA", "Distribution of Genes Detected per Cell", "Number of Unique Genes"),
                         "Percent Mitochondria" = c("percentMito", "Distribution of Mito GE per Cell", "Mitochondrial Ratio"),
                         "ADT Count Per Cell" = c("nCount_ADT", "ADT Counts per Cell","Number of Counts"),
                         "Unique ADTs Per Cell" = c("nFeature_ADT", "Distribution of CITE-seq Antibodies per Cell", "Number of Unique Antibodies")
        )
        
        #creates list of x-coordinates for quantiles of data. 
        quant <- quantile(x = myso@meta.data[,params[1]],
                          probs = c(0.5,0.75,0.95),
                          na.rm = TRUE)
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(myso@meta.data[[color]])))
        
        #generate base plot template with features that all QA distribution plots will have
        base_distrib_plot <- myso@meta.data %>%
          ggplot2::ggplot(aes(x = eval(parse(text = params[1])), fill = eval(parse(text = color)), color = eval(parse(text = color)))) + #eval(parse(text = x)) necessary to turn string into variable name format
          labs(fill = color, color = color) + 
          theme(plot.title = element_text(hjust=0.5)) +
          ggtitle(params[2]) +
          xlab(params[3]) +
          geom_vline(xintercept = quant, size = 0.5, alpha = 0.5, linetype = "dashed", color = "grey30") +
          scale_color_manual(values = custom_palette) +
          scale_fill_manual(values = custom_palette)
        
        #initialize QA distribution plot before if/else statements below so that the plot object can be accessed outside of the if/else statements
        final_distrib_plot <- base_distrib_plot
        
        #create density/bar plot for selected input. If integrated object is uploaded, then the original identity of the cells will separate into graphs per sample
        if (input$QA %in% "ADT Count Per Cell") {
          final_distrib_plot <- base_distrib_plot + scale_x_log10() + geom_density(alpha = 0.25) + ylab("Density")
        } 
        else if (input$QA %in% "Unique ADTs Per Cell") {
          final_distrib_plot <- base_distrib_plot + geom_bar(alpha = 0.5, position = "dodge") + ylab("Frequency")
        }
        else {
          final_distrib_plot <- base_distrib_plot + geom_density(alpha = 0.25) + ylab("Density")
        }
        
        #show distribution plot
        final_distrib_plot %>% 
          plotly::ggplotly() %>% 
          config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          layout(title = list(font = list(size = 14)), hovermode = FALSE) 
      })
      
      
      # ----- QA box plot -----
      # reactive box plot 
      box_plot <- reactive({
        
        req(input$rds_input_file, input$QA, input$color_qa)
        color <- input$color_qa
        
        #This assigns the params variable a list of strings that act a varying parameters depending on input of QA
        params <- switch(input$QA,
                         "RNA Count Per Cell" = c("nCount_RNA", "Distribution of Counts per Cell", "Number of Counts"),
                         "Gene Count Per Cell" = c("nFeature_RNA", "Distribution of Genes Detected per Cell", "Number of Unique Genes"),
                         "Percent Mitochondria" = c("percentMito", "Distribution of Mito GE per Cell", "Mitochondrial Ratio"),
                         "ADT Count Per Cell" = c("nCount_ADT", "ADT Counts per Cell","Number of Counts"),
                         "Unique ADTs Per Cell" = c("nFeature_ADT", "Distribution of CITE-seq Antibodies per Cell", "Number of Unique Antibodies")
        )
        
        #instead of switch and hardcoded values, try colnames() of the input data so user can select whatever cols are in their data
        
        #creates list of x-coordinates for quantiles of data. 
        quant <- quantile(x = myso@meta.data[,params[1]],
                          probs = c(0.5,0.75,0.95),
                          na.rm = TRUE)
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(myso@meta.data[[color]])))
        
        #generate base plot template with features that all QA boxplots will have
        base_box_plot <- myso@meta.data %>%
          ggplot2::ggplot(aes(x = eval(parse(text = color)), y = eval(parse(text = params[1])), fill = eval(parse(text = color)), color = eval(parse(text = color)))) + #eval(parse(text = x)) necessary to turn string into variable name format
          labs(fill = color, color = color) + 
          geom_boxplot(alpha = 0.5, width=0.5) + 
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
          theme(plot.title = element_text(hjust = 0.5)) +
          ggtitle(params[2]) +
          xlab("Sample") +
          ylab(params[3]) +
          geom_violin(alpha = 0.2) +
          geom_hline(yintercept = quant, size = 0.5, alpha = 0.5, linetype = "dashed", color = "grey30") +
          scale_color_manual(values = custom_palette) +
          scale_fill_manual(values = custom_palette)
        
        #initialize QA box plot before if/else statements below so that the plot object can be accessed outside of the if/else statements
        final_box_plot <- base_box_plot
        
        #create box plot for selected input. If integrated object is uploaded, then the original identity of the cells will separate into boxes per sample
        if (input$QA %in% "ADT Count Per Cell") {
          final_box_plot <- base_box_plot + scale_y_log10()
        } 
        
        #show box plot
        final_box_plot %>% 
          plotly::ggplotly() %>% 
          config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          layout(title = list(font = list(size = 14)), hovermode = FALSE)
      })
      
      # ----- render QA plots -----
      # render the reactive plotly QA plots 
      output$distrib_plot <- plotly::renderPlotly({ distribution_plot() })
      output$box_plot <- plotly::renderPlotly({ box_plot() })
      
    }) #end of QA observe() wrapper
    
    
    # ---------- ***** Clustering ***** ---------- 
    observe({
      
      #changes the selectInput "reduction" dropdown contents to include all reductions in Seurat Object
      updateSelectInput(
        session = session,
        inputId = "reduction",
        choices = sort(names(myso@reductions)),
        selected = dplyr::last(sort(names(myso@reductions)))
      )
      
      #changes selectInput "color" dropdown contents to include all metadata in Seurat Object
      updateSelectInput(
        session = session,
        inputId = "color1",
        choices = colnames(myso@meta.data[lapply(myso@meta.data, class) %in% c("factor", "character")]),
        selected = colnames(myso@meta.data[lapply(myso@meta.data, class) %in% c("factor", "character")])[1]
      )
      
      # ----- Reactive 2D reduction graph -----
      reduc_plot <- reactive({
        req(input$rds_input_file, input$reduction, input$color1)
        
        #create string for reduction to plot
        reduc <- input$reduction
        
        #selected metadata to color clusters by
        color <- input$color1
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(myso@meta.data[[color]])))
        
        #create dataframe from reduction selected
        cell_data <- data.frame(eval(parse(text = paste0("myso@reductions$", reduc, "@cell.embeddings"))))
        
        #create list containing all column names of cell_data
        cell_col <- colnames(cell_data)
        
        #show plot
        plotly::plot_ly(cell_data, 
                        x = ~cell_data[,1], y = ~cell_data[,2],
                        customdata = rownames(cell_data), #customdata is printed in cell selection and used to find metadata
                        color = ~eval(parse(text = paste0("myso@meta.data$", color))), #color by selected metadata in object; need to incorporate 
                        colors = custom_palette,
                        type = "scatter", 
                        mode = "markers",
                        marker = list(size = 3, width = 2),
                        source = "A") %>%
          
          config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          #Layout changes the aesthetic of the plot
          layout(
            title = toupper(reduc),
            xaxis = list(title = cell_col[1]),
            yaxis = list(title = cell_col[2]),
            dragmode = "select") %>% 
          #Determines the mode of drag interactions. "select" and "lasso" apply only to scatter traces with markers or text. "orbit" and "turntable" apply only to 3D scenes.
          
          event_register("plotly_selected")
      })
      
      # ----- Reactive 3D reduction graph -----
      reduc_plot_3d <- reactive({
        req(input$rds_input_file, input$reduction, input$color1)
        
        #create string for reduction to plot
        reduc <- input$reduction
        
        #selected metadata to color clusters by
        color <- input$color1
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(myso@meta.data[[color]])))
        
        #create dataframe from reduction selected
        cell_data <- data.frame(eval(parse(text = paste0("myso@reductions$", reduc, "@cell.embeddings"))))
        
        #create list containing all column names of cell_data
        cell_col <- colnames(cell_data)
        
        #show plot
        plotly::plot_ly(cell_data, 
                        x = ~cell_data[,1], y = ~cell_data[,2], z = ~cell_data[,3],
                        customdata = rownames(cell_data), #customdata is printed in cell selection and used to find metadata
                        color = ~eval(parse(text = paste0("myso@meta.data$", color))), #color by selected metadata in object; need to incorporate 
                        colors =  custom_palette,
                        type = "scatter3d", 
                        mode = "markers",
                        marker = list(size = 2, width = 1)) %>%
          
          config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          #Layout changes the aesthetic of the plot
          layout(
            title = toupper(paste(reduc, "(3D)")),
            scene = list(xaxis = list(title = cell_col[1]),
                         yaxis = list(title = cell_col[2]),
                         zaxis = list(title = cell_col[3])),
            dragmode = "orbit") #Determines the mode of drag interactions. "select" and "lasso" apply only to scatter traces with markers or text. "orbit" and "turntable" apply only to 3D scenes.
      })
      
      # ----- render clustering plots -----
      output$output_2dplot_1 <- plotly::renderPlotly({ reduc_plot() })
      output$output_3dplot_1 <- plotly::renderPlotly({ reduc_plot_3d() })
      
      # ----- datatable of metadata for cells selected in plotly -----
      # `server = FALSE` helps make it so that user can copy entire datatable to clipboard, not just the rows that are currently visible on screen
      output$cluster_pg_selected <- DT::renderDT(server = FALSE, {
        #currently returns every column of metadata dataframe. May want to select specific columns in the future
        selected_metadata_df <- myso@meta.data[event_data("plotly_selected", source = "A")$customdata, ]
        cluster_dt <- DT::datatable(selected_metadata_df, 
                                    rownames = TRUE,
                                    selection = "none", #make it so no rows can be selected (bc we currently have no need to select rows)
                                    extensions = c("Buttons", "Scroller", "FixedColumns"),
                                    options = list(deferRender = TRUE,
                                                   scroller = TRUE,
                                                   scrollY = 400,
                                                   scrollX = TRUE,
                                                   dom = "lfrtipB",
                                                   buttons = c("copy", "print"),
                                                   fixedColumns = list(leftColumns = 1)
                                    )
        ) 
        if (is.null(cluster_dt)) "Brushed points appear here (double-click to clear)" else cluster_dt
      })
      
    })  #end of clustering observed wrapper 
    
    
    # ------- ***** Co-Expression ***** ----------
    
    ### choose assay for x and y axes and then display dropdowns
    
    observe({
      
      # ----- update/render UI elements -----
      
      updateSelectInput(
        #session = session,
        inputId = "reduction_expr",
        choices = sort(names(myso@reductions)),
        selected = dplyr::last(sort(names(myso@reductions)))
      )
      
      output$Assay_x_axis = renderUI({
        selectInput(inputId = "Assay_x_axis",
                    label = "Choose assay for x-axis colorscale:",
                    choices = sort(names(myso@assays)),
                    selected = sort(names(myso@assays))[1]
        )
      })
      
      output$Assay_y_axis = renderUI({
        selectInput(inputId = "Assay_y_axis",
                    label = "Choose assay for y-axis colorscale:",
                    choices = sort(names(myso@assays)),
                    selected = sort(names(myso@assays))[1]
        )
      })
      
      
      output$x_axis_feature = renderUI({
        req(input$Assay_x_axis)
        feature_path <- paste0('myso@assays$', input$Assay_x_axis, '@data')
        selectInput(inputId = "x_axis_feature",
                    label = "Choose feature for x-axis colorscale:",
                    choices = rownames(eval(parse(text=feature_path))),
                    selected = rownames(eval(parse(text=feature_path)))[1]
        )
      })
      
      output$y_axis_feature = renderUI({
        req(input$Assay_y_axis)
        feature_path <- paste0('myso@assays$', input$Assay_y_axis, '@data')
        selectInput(inputId = "y_axis_feature",
                    label = "Choose feature for y-axis colorscale:",
                    choices = rownames(eval(parse(text=feature_path))),
                    selected = rownames(eval(parse(text=feature_path)))[2]
        )
      })
      
      
      # ----- 1D gene/ADT expression reactive reduction graph -----
      expr_reduc_plot_1d <- eventReactive(
        list(
          #list of input events that can trigger reactive 
          input$rds_input_file,
          input$reduction_expr,
          input$x_axis_feature
        ), 
        {
          req(input$rds_input_file, input$reduction_expr, input$x_axis_feature)
          
          #create string for reduction to plot
          reduc <- input$reduction_expr
          
          #selected metadata to color clusters by
          color_x <- input$x_axis_feature
          
          SeuratObject::DefaultAssay(myso) <- input$Assay_x_axis
          count_data <- SeuratObject::FetchData(object = myso, vars = color_x, slot = "data")
          
          # #interpolate the base color palette so that exact number of colors in custom palette is 3 (for no, low, and high expression values).
          # custom_palette <- colorRampPalette(c("#E4E4E4", "darkmagenta"))(length(unique(count_data[,input$x_axis_feature])))
          
          #create dataframe from reduction selected
          cell_data <- data.frame(eval(parse(text = paste0("myso@reductions$", reduc, "@cell.embeddings"))))
          
          #create list containing all column names of cell_data
          cell_col <- colnames(cell_data)
          
          #show plot
          plotly::plot_ly(cell_data, 
                          x = ~cell_data[,1], y = ~cell_data[,2],
                          customdata = rownames(cell_data), #customdata is printed in cell selection and used to find metadata
                          #colors = custom_palette,
                          type = "scatter", 
                          mode = "markers",
                          marker = list(size = 3,
                                        color = ~count_data[, color_x], 
                                        colorbar = list(title = color_x,
                                                        len = 0.5),
                                        colorscale = "Viridis",
                                        reversescale = TRUE)) %>%

            config(toImageButtonOptions = list(format = "png",
                                               scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
            ) %>%
            #Layout changes the aesthetic of the plot
            layout(
              title = toupper(reduc),
              xaxis = list(title = cell_col[1]),
              yaxis = list(title = cell_col[2]),
              dragmode = "select")
            #Determines the mode of drag interactions. "select" and "lasso" apply only to scatter traces with markers or text. "orbit" and "turntable" apply only to 3D scenes.
      })
      
      
      # ----- 2D gene/ADT expression reactive reduction graph -----
      expr_reduc_plot_2d <- eventReactive(
        list(
          #list of input events that can trigger reactive plot
          input$rds_input_file,
          input$reduction_expr,
          input$x_axis_feature,
          input$y_axis_feature
        ), 
        {
          req(input$rds_input_file, input$reduction_expr, input$x_axis_feature, input$y_axis_feature)
          
          #create string for reduction to plot
          reduc <- input$reduction_expr
          
          #selected metadata to color clusters by
          color_x <- input$x_axis_feature
          color_y <- input$y_axis_feature
          
          SeuratObject::DefaultAssay(myso) <- input$Assay_x_axis
          count_data_x <- SeuratObject::FetchData(object = myso, vars = color_x, slot = "data")
          
          SeuratObject::DefaultAssay(myso) <- input$Assay_y_axis
          count_data_y <- SeuratObject::FetchData(object = myso, vars = color_y, slot = "data")
          
          #create dataframe from reduction selected
          cell_data <- data.frame(eval(parse(text = paste0("myso@reductions$", reduc, "@cell.embeddings"))))
          
          #create list containing all column names of cell_data
          cell_col <- colnames(cell_data)
          
          #show plot
          plotly::plot_ly(cell_data,
                          x = ~cell_data[,1], y = ~cell_data[,2],
                          customdata = rownames(cell_data), #customdata is printed in cell selection and used to find metadata
                          type = "scatter",
                          mode = "markers",
                          marker = list(size = 3,
                                        color = ~count_data_x[, color_x],
                                        colorbar = list(title = color_x,
                                                        len = 0.5, 
                                                        yanchor = "bottom"),
                                        colorscale = "Reds"),
                          opacity = 1) %>%
            
            add_markers(x = ~cell_data[,1], y = ~cell_data[,2],
                        marker = list(size=3,
                                      color = ~count_data_y[, color_y],
                                      colorbar = list(title = color_y,
                                                      len = 0.5,
                                                      yanchor = "top"),
                                      colorscale = "Blues",
                                      reversescale = TRUE #essential for how plotly works
                                      ),
                        opacity = 0.5) %>%

          config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
            #Layout changes the aesthetic of the plot
            layout(
              showlegend = FALSE,
              title = toupper(reduc),
              xaxis = list(title = cell_col[1]),
              yaxis = list(title = cell_col[2]),
              dragmode = "select") 
            #Determines the mode of drag interactions. "select" and "lasso" apply only to scatter traces with markers or text. "orbit" and "turntable" apply only to 3D scenes.
      })
      
      
      #render reactive reduction plots
      output$exploration_reduct_1d <- renderPlotly({
        expr_reduc_plot_1d()
      })
      
      output$exploration_reduct_2d <- renderPlotly({
        expr_reduc_plot_2d()
      })
      
    })  # belongs to OBSERVE WRAPPER for co-expression tab
    
    
    # ---------- ***** Gating ***** ----------
    observe({
      
      # ----- update/render UI elements -----
      output$Assay <- renderUI({
        selectInput(inputId = "Assay",
                    label = "Choose assay:",
                    choices = sort(names(myso@assays)),
                    selected = sort(names(myso@assays))[1])
      })
      
      output$x_feature <- renderUI({
        req(input$Assay)
        feature_path <- paste0('myso@assays$', input$Assay, '@data')
        selectInput(
          inputId = "x_feature",
          label = "Choose x-axis feature:",
          choices = rownames(eval(parse(text=feature_path))),
          selected = rownames(eval(parse(text=feature_path)))[1])
      })
      
      output$y_feature <- renderUI({
        req(input$Assay)
        feature_path <- paste0('myso@assays$', input$Assay, '@data')
        selectInput(
          inputId = "y_feature",
          label = "Choose y-axis feature:",
          choices = rownames(eval(parse(text=feature_path))),
          selected = rownames(eval(parse(text=feature_path)))[2])
      })
      
      #changes the selectInput "reduction" dropdown contents to include all reductions in Seurat Object
      updateSelectInput(
        session = session,
        inputId = "reduction_g",
        choices = sort(names(myso@reductions)),
        selected = dplyr::last(sort(names(myso@reductions)))
      )
      updateSelectInput(
        session = session,
        inputId = "color2",
        choices = colnames(myso@meta.data[lapply(myso@meta.data, class) %in% c("factor", "character")]),
        selected = colnames(myso@meta.data[lapply(myso@meta.data, class) %in% c("factor", "character")])[1]
      )
      
      # ----- Last-clicked buttons tracker -----
      # keep track of last 2 buttons clicked (either gate button, reset button, or clear-all-gates button) in gating tab bc this will determine which data to use for reactive ADT scatterplot
      # second-to-last clicked button will determine what data will be used for input gating cells and input cell count
      # default is "NA" (using NULL causes problems when rendering reactive featurescatter)
      # but if reset, gate, or clear-all-gates button is clicked, then the values of last_buttons_clicked will reflect that click
      
      last_buttons_clicked <- reactiveValues(last = "NA", second_to_last = "NA")
      
      observeEvent(input$reset_adt_scatter, {
        req(input$x_feature, input$y_feature) #needed so that nothing happens if button is pressed before initial scatterplot renders upon startup
        last_buttons_clicked$second_to_last <- last_buttons_clicked$last
        last_buttons_clicked$last <- "reset_button"
      })
      observeEvent(input$gate, {
        req(input$x_feature, input$y_feature) #needed so that nothing happens if button is pressed before initial scatterplot renders upon startup
        last_buttons_clicked$second_to_last <- last_buttons_clicked$last
        last_buttons_clicked$last <- "gate_button"
      })
      
      
      # ----- reactive scatterplot of assay features -----
      reactive_featurescatter <- eventReactive(
        list(
          #list of input events that can trigger reactive featurescatter
          input$rds_input_file,
          input$gate,
          input$reset_adt_scatter,
          input$clear_all_gates,
          input$x_feature,
          input$y_feature,
          input$gating_pg_table_rows_selected
        ), 
        { 
          #code to execute when one of the above input events occurs
          req(input$x_feature, input$y_feature)
          
          SeuratObject::DefaultAssay(myso) <- input$Assay
          count_data <- SeuratObject::FetchData(object = myso, vars = c(input$x_feature, input$y_feature), slot = "data")
          
          #generate dataframe for custom colorscale for contour plot, where each hex color code is mapped to a specific z-value between 0 and 1 (inclusive)
          #colorscale needs to be in this format for Plotly's add_histogram2dcontour(colorscale = ...) parameter
          customColorScale <- data.frame(
            z = c(0.0, 0.20, 0.40, 0.60, 0.80, 1.0),
            col = c("#FFFFFF", "#4564FE", "#76EFFF", "#FFF900", "#FFA300", "#FF1818")
          )
          
          #initialize a base_scatterplot variable before if/else statements below so that the plot object can be accessed outside of the if/else statements
          base_scatterplot <- NULL
          
          # generate a base scatterplot based on user input
          if ((last_buttons_clicked$last == "NA" | last_buttons_clicked$last == "reset_button" | last_buttons_clicked$last == "clear_all_gates_button") & is.null(input$gating_pg_table_rows_selected)) {
            base_scatterplot <- plotly::plot_ly(count_data,
                                                x = ~count_data[,input$x_feature], y = ~count_data[,input$y_feature],
                                                customdata = rownames(count_data),
                                                mode = "markers",
                                                source = "C") %>% 
              add_histogram2dcontour(showscale = FALSE, ncontours = 10, colorscale = customColorScale, 
                                     contours = list(coloring='heatmap')) %>%
              add_markers(x = count_data[,input$x_feature],
                          y = count_data[,input$y_feature],
                          marker = list(size=2),
                          color = I("black"),
                          alpha = 0.6) 
          }
          else {
            selected_cell_barcodes <- NULL
            count_data_subset <- NULL
            
            if (last_buttons_clicked$last == "gate_button" & is.null(input$gating_pg_table_rows_selected)) {
              selected_cell_barcodes <- event_data("plotly_selected", source = "C")$customdata
              count_data_subset <- count_data[rownames(count_data) %in% selected_cell_barcodes, ]
            }
            else if (!is.null(input$gating_pg_table_rows_selected)) {
              selected_cell_barcodes <- gate_list()[[selected_gate()]]@subset_cells[[1]]
              count_data_subset <- count_data[rownames(count_data) %in% selected_cell_barcodes, ]
            }
            base_scatterplot <- plotly::plot_ly(count_data_subset,
                                                x = ~count_data_subset[,input$x_feature], y = ~count_data_subset[,input$y_feature],
                                                customdata = rownames(count_data_subset),
                                                mode = "markers",
                                                source = "C") %>% 
              add_histogram2dcontour(showscale=FALSE, ncontours=10, colorscale = customColorScale, 
                                     contours = list(coloring='heatmap')) %>%
              add_markers(x = count_data_subset[,input$x_feature], 
                          y = count_data_subset[,input$y_feature], 
                          marker = list(size=2.5), 
                          color = I("black"), 
                          alpha = 0.6)
          }
          
          # add configuration and layout options to base scatterplot, register selection events
          base_scatterplot %>%
            config(toImageButtonOptions = list(format = "png",
                                               scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
            ) %>%
            #Layout changes the aesthetic of the plot
            layout(
              title = "Normalized Feature Scatter Plot",
              xaxis = list(title = input$x_feature),
              yaxis = list(title = input$y_feature),
              showlegend = FALSE,
              dragmode = "select") %>% #Determines the mode of drag interactions. "select" and "lasso" apply only to scatter traces with markers or text. "orbit" and "turntable" apply only to 3D scenes.
            
            event_register("plotly_selected")
        })
      
      # ----- reactive gating 2D reduction graph -----
      gating_reduc_plot <- reactive({
        req(input$rds_input_file, input$color2, input$reduction_g)
        
        #create string for reduction to plot
        reduc <- input$reduction_g
        
        #selected metadata to color clusters by
        color <- input$color2
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(myso@meta.data[[color]])))
        
        #create dataframe from reduction selected
        cell_data <- data.frame(eval(parse(text = paste0("myso@reductions$", reduc, "@cell.embeddings"))))
        
        #create list containing all column names of cell_data
        cell_col <- colnames(cell_data)
        
        # initialize selected_cells before if/else statements so it can be accessed outside of the statements
        selected_cells <- NULL
        
        if (!is.null(input$gating_pg_table_rows_selected)) {
          selected_cells <- gate_list()[[selected_gate()]]@subset_cells[[1]]
        }
        else {
          selected_cells <- event_data("plotly_selected", source = "C")$customdata
        }
        
        if (is.null(selected_cells)) {
          plotly_color_list <- c(paste0("myso@meta.data$", color), custom_palette)
        } 
        else {
          plotly_color_list <- c("rownames(cell_data) %in% selected_cells", 'c("grey", "black")')
        }
        
        plotly::plot_ly(cell_data, 
                        x = ~cell_data[,1], y = ~cell_data[,2],
                        customdata = rownames(cell_data), #customdata is printed in cell selection and used to find metadata
                        color = ~eval(parse(text = plotly_color_list[1])), #color by selected metadata in object; need to incorporate 
                        colors = ~eval(parse(text = plotly_color_list[2])),
                        type = 'scatter', 
                        mode = 'markers',
                        marker = list(size = 3, width=2)) %>%
          
          config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          #Layout changes the aesthetic of the plot
          layout(
            title = toupper(reduc),
            xaxis = list(title = cell_col[1]),
            yaxis = list(title = cell_col[2]))
      })
      
      
      # ----- render gating plots -----
      #render reactive feature scatterplot
      output$featurescatter_2d <- plotly::renderPlotly({ reactive_featurescatter() })
      
      #render reactive UMAP
      output$gating_reduc_2d <- plotly::renderPlotly({ gating_reduc_plot() })
      
      
      # ----- generate Gate objects -----
      gate_reactive_values <- reactiveValues()
      gate_list <- get_reactive_gate_list(gate_reactive_values)
      counter_reactive <- reactiveVal(as.integer(0))
      
      # holds value of gate that is selected from gating datatable in case user wants to go back and re-gate from selected gate
      selected_gate <- reactiveVal(NULL)
      
      #events triggered by clicking gate button
      observeEvent(input$gate, {
        req(input$x_feature, input$y_feature)
        
        SeuratObject::DefaultAssay(myso) <- input$Assay
        count_data <- SeuratObject::FetchData(object = myso, vars = c(input$x_feature, input$y_feature), slot = "data")
        
        # get plotly event data
        sel <- event_data("plotly_selected", source = "C")
        brushed_coords <- event_data("plotly_brushed", source = "C")

        # increment counter every time gate button is clicked
        counter <- as.integer(counter_reactive() + 1)
        counter_reactive(counter)

        # create gate object based on UI input
        gate_reactive_values[[paste0("gate_", counter)]] <- create_gate_from_input(is_forward_gating = TRUE,
                                                                                   assay_count_data = count_data,
                                                                                   gate_counter = counter,
                                                                                   reactive_gate_list = gate_list(),
                                                                                   reactive_selected_gate = selected_gate(),
                                                                                   reactive_last_buttons_clicked = last_buttons_clicked)
      }) #end of gating logic for events triggered by clicking gate button
      
      
      # ----- generate reactive gating dataframe -----
      reactive_gating_df <- reactive({
        
        #gating dataframe (full version) for download
        full_gating_df <- create_gating_df()
        
        if (!rlang::is_empty(names(gate_list()))) {
          full_gating_list <- lapply(names(gate_list()), update_gating_df, reactive_gate_list = gate_list(), temp_gating_df = full_gating_df)
          full_gating_df <- do.call(rbind, full_gating_list)
        }
        full_gating_df
      })
      
      # returns data table of summary data for each gate object
      output$gating_pg_table <- DT::renderDT({
        # data table is generated here
        gating_dt <- create_gating_dt(reactive_gating_df())
        if (is.null(gating_dt)) "Brushed points appear here (double-click to clear)" else gating_dt
      })
      
      
      # ----- datatable event handlers -----
      
      # Update back-end gating data when user edits name of cell subset in front-end datatable
      observeEvent(input$gating_pg_table_cell_edit, {
        cell_edit_data <- input$gating_pg_table_cell_edit
        gate_id <- reactive_gating_df()$Gate_ID[cell_edit_data$row]
        # [input$gating_pg_table_cell_edit$row, 1]
        gate_reactive_values[[gate_id]]@name_subset_cells <- cell_edit_data$value
      })
      
      # when a datatable row is selected, set the selected_gate reactive value to the gate that was selected in the datatable
      # this ensures that if the user wants to re-gate based on this gate, then the corresponding input gate stats are correct
      observeEvent(input$gating_pg_table_rows_selected, {
        row_index <- input$gating_pg_table_rows_selected
        selected_gate(reactive_gating_df()$Gate_ID[row_index])
        
        #update assay, x and y axis dropdowns to reflect the same axes shown in selected gate, so that if user gates from this selected gate, the right axes are recorded for the new gate
        updateSelectInput(
          session = session,
          inputId = "Assay",
          selected = gate_list()[[selected_gate()]]@assay_name
        )
        updateSelectInput(
          session = session,
          inputId = "x_feature",
          selected = gate_list()[[selected_gate()]]@x_axis
        )
        updateSelectInput(
          session = session,
          inputId = "y_feature",
          selected = gate_list()[[selected_gate()]]@y_axis
        )
      })
      
      # when clear all gates button is clicked or new rds file is uploaded, reset gating info
      observeEvent(
        list(
          #list of input events that can trigger resetting of gating info
          input$clear_all_gates,
          input$rds_input_file
        ), 
        {
          counter_reactive(as.integer(0))
          
          # reset all values in gate_reactive_values to NULL
          if (!rlang::is_empty(names(gate_reactive_values))) {
            gate_reactive_values <- lapply(names(gate_reactive_values), set_gates_to_null, local_gate_reactive_values = gate_reactive_values)
          }
          
          #update last-buttons tracker
          last_buttons_clicked$second_to_last <- last_buttons_clicked$last
          last_buttons_clicked$last <- "clear_all_gates_button"
        })
      
      
      # ----- gating data download handlers -----
      output$download_as_list_rds <- downloadHandler(
        filename = "gate_info_list.rds",
        content = function(file) {
          saveRDS(gate_list(), file = file)
        }
      )
      output$download_as_df_rds <- downloadHandler(
        filename = "gate_info_df.rds",
        content = function(file) {
          saveRDS(reactive_gating_df(), file = file)
        }
      )
      
    })  # belongs to OBSERVE WRAPPER for gating tab
    
    
    # ---------- ***** BackGating ***** ----------
    observe({
      
      # ----- update/render UI elements -----
      output$Assay_bg <- renderUI({
        selectInput(inputId = "Assay_bg",
                    label = "Choose assay:",
                    choices = sort(names(myso@assays)),
                    selected = sort(names(myso@assays))[1])
      })
      
      output$x_feature_bg <- renderUI({
        req(input$Assay_bg)
        feature_path <- paste0('myso@assays$', input$Assay_bg, '@data')
        selectInput(
          inputId = "x_feature_bg",
          label = "Choose x-axis feature:",
          choices = rownames(eval(parse(text=feature_path))),
          selected = rownames(eval(parse(text=feature_path)))[1])
      })
      
      output$y_feature_bg <- renderUI({
        req(input$Assay_bg)
        feature_path <- paste0('myso@assays$', input$Assay_bg, '@data')
        selectInput(
          inputId = "y_feature_bg",
          label = "Choose y-axis feature:",
          choices = rownames(eval(parse(text=feature_path))),
          selected = rownames(eval(parse(text=feature_path)))[2])
      })
      
      #changes the selectInput "reduction" dropdown contents to include all reductions in Seurat Object
      updateSelectInput(
        session = session,
        inputId = "reduction_bg",
        choices = sort(names(myso@reductions)),
        selected = dplyr::last(sort(names(myso@reductions)))
      )
      updateSelectInput(
        session = session,
        inputId = "color2_bg",
        choices = colnames(myso@meta.data[lapply(myso@meta.data, class) %in% c("factor", "character")]),
        selected = colnames(myso@meta.data[lapply(myso@meta.data, class) %in% c("factor", "character")])[1]
      )
      
      # ----- Last-clicked buttons tracker -----
      # keep track of last 2 buttons clicked (either gate button, reset button, or clear-all-gates button) in gating tab bc this will determine which data to use for reactive ADT scatterplot
      # second-to-last clicked button will determine what data will be used for input gating cells and input cell count
      # default is "NA" (using NULL causes problems when rendering reactive featurescatter)
      # but if reset, gate, or clear-all-gates button is clicked, then the values of last_buttons_clicked will reflect that click
      
      last_buttons_clicked_bg <- reactiveValues(last = "NA", second_to_last = "NA")
      
      observeEvent(input$reset_adt_scatter_bg, {
        req(input$x_feature_bg, input$y_feature_bg) #needed so that nothing happens if button is pressed before initial scatterplot renders upon startup
        last_buttons_clicked_bg$second_to_last <- last_buttons_clicked_bg$last
        last_buttons_clicked_bg$last <- "reset_button"
      })
      observeEvent(input$gate_bg, {
        req(input$x_feature_bg, input$y_feature_bg) #needed so that nothing happens if button is pressed before initial scatterplot renders upon startup
        last_buttons_clicked_bg$second_to_last <- last_buttons_clicked_bg$last
        last_buttons_clicked_bg$last <- "gate_button"
      })
      
      # ----- backgate reactive scatterplot of assay features -----
      reactive_featurescatter_bg <- eventReactive(
        list(
          #list of input events that can trigger reactive featurescatter
          input$rds_input_file,
          input$gate_bg,
          input$reset_adt_scatter_bg,
          input$clear_all_gates_bg,
          input$x_feature_bg,
          input$y_feature_bg,
          input$gating_pg_table_bg_rows_selected,
          event_data("plotly_selected", source = "D")
        ), 
        {
          #code to execute when one of the above input events occurs
          req(input$x_feature_bg, input$y_feature_bg)
          
          SeuratObject::DefaultAssay(myso) <- input$Assay_bg
          count_data <- SeuratObject::FetchData(object = myso, vars = c(input$x_feature_bg, input$y_feature_bg), slot = "data")
          selected_cell_barcodes <- NULL
          
          if (is.null(input$gating_pg_table_bg_rows_selected)) {
            selected_cell_barcodes <- event_data("plotly_selected", source = "D")$customdata
          }
          else {
            selected_cell_barcodes <- gate_list_bg()[[selected_gate_bg()]]@subset_cells[[1]]
          }
          
          plotly::plot_ly(count_data,
                          x = ~count_data[,input$x_feature_bg], y = ~count_data[,input$y_feature_bg],
                          customdata = rownames(count_data),
                          mode = "markers",
                          color = rownames(count_data) %in% selected_cell_barcodes, #color cells by whether they're in the selection or not
                          colors = c("grey", "black")) %>% 
            add_histogram2dcontour(showscale = FALSE, ncontours = 10,
                                   colorscale = NULL,
                                   contours = list(coloring='none'),
                                   color = I("magenta3"),
                                   size = I(1.5)) %>%
            add_markers(x = count_data[,input$x_feature_bg],
                        y = count_data[,input$y_feature_bg],
                        marker = list(size=2),
                        alpha = 0.6) %>%
            config(toImageButtonOptions = list(format = "png",
                                               scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
            ) %>%
            #Layout changes the aesthetic of the plot
            layout(
              title = "Normalized Feature Scatter Plot",
              xaxis = list(title = input$x_feature_bg),
              yaxis = list(title = input$y_feature_bg),
              showlegend = FALSE)
        })
      
      # ----- backgate reactive reduction plot of selected cell features -----
      gating_reduc_plot_bg <- reactive({
        req(input$rds_input_file, input$color2_bg, input$reduction_bg)
        
        #creates string for reduction to plot
        reduc <- input$reduction_bg
        
        #selected metadata to color clusters by
        color <- input$color2_bg
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(myso@meta.data[[color]])))
        plotly_color_list <- c(paste0("myso@meta.data$", color), custom_palette)
        
        #creates dataframe from reduction selected
        cell_data <- data.frame(eval(parse(text = paste0("myso@reductions$", reduc, "@cell.embeddings"))))
        
        #creates list containing all column names of cell_data
        cell_col <- colnames(cell_data)
        
        plotly::plot_ly(cell_data, 
                        x = ~cell_data[,1], y = ~cell_data[,2],
                        customdata = rownames(cell_data), #customdata is printed in cell selection and used to find metadata
                        color = ~eval(parse(text = plotly_color_list[1])), #color by selected metadata in object; need to incorporate 
                        colors = ~eval(parse(text = plotly_color_list[2])),
                        type = 'scatter',
                        mode = 'markers',
                        marker = list(size = 3, width=2),
                        source = "D") %>%
          config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          #Layout changes the aesthetic of the plot
          layout(title = toupper(reduc),
                 xaxis = list(title = cell_col[1]),
                 yaxis = list(title = cell_col[2]),
                 dragmode = "select") %>%
          event_register("plotly_selected")
      })
      
      
      # ----- render backgating plots -----
      #render reactive feature scatterplot
      output$featurescatter_2d_bg <- plotly::renderPlotly({reactive_featurescatter_bg()})
      
      #render reactive UMAP
      output$gating_reduc_2d_bg <- plotly::renderPlotly({gating_reduc_plot_bg()})
      
      
      # ----- generate Gate objects -----
      gate_reactive_values_bg <- reactiveValues()
      gate_list_bg <- get_reactive_gate_list(gate_reactive_values_bg)
      
      counter_reactive_bg <- reactiveVal(as.integer(0))
      
      # holds value of gate that is selected from gating datatable in case user wants to go back and re-gate from selected gate
      selected_gate_bg <- reactiveVal(NULL)
      
      #events triggered by clicking gate button
      observeEvent(input$gate_bg, {
        req(input$x_feature_bg, input$y_feature_bg)
        
        SeuratObject::DefaultAssay(myso) <- input$Assay_bg
        count_data <- SeuratObject::FetchData(object = myso, vars = c(input$x_feature_bg, input$y_feature_bg), slot = "data")
        
        # get plotly event data
        sel <- event_data("plotly_selected", source = "D")
        brushed_coords <- event_data("plotly_brushed", source = "D")

        # increment counter every time gate button is clicked
        counter <- as.integer(counter_reactive_bg() + 1)
        counter_reactive_bg(counter)
        
        # create gate object based on UI input
        gate_reactive_values_bg[[paste0("gate_", counter)]] <- create_gate_from_input(is_forward_gating = FALSE,
                                                                                     assay_count_data = count_data,
                                                                                     gate_counter = counter,
                                                                                     reactive_gate_list = gate_list_bg(),
                                                                                     reactive_selected_gate = selected_gate_bg(),
                                                                                     reactive_last_buttons_clicked = last_buttons_clicked_bg)
      }) #end of gating logic for events triggered by clicking gate button
      
      
      # ----- generate reactive backgating dataframe -----
      reactive_gating_df_bg <- reactive({
        
        #gating dataframe (full version) for download
        full_gating_df <- create_gating_df()
        
        if (!rlang::is_empty(names(gate_list_bg()))) {
          full_gating_list <- lapply(names(gate_list_bg()), update_gating_df, reactive_gate_list = gate_list_bg(), temp_gating_df = full_gating_df)
          full_gating_df <- do.call(rbind, full_gating_list)
        }
        full_gating_df
      })
      
      # returns data table of summary data for each gate object
      output$gating_pg_table_bg <- DT::renderDT({
        # data table is generated here
        gating_dt <- create_gating_dt(reactive_gating_df_bg())
        if (is.null(gating_dt)) "Brushed points appear here (double-click to clear)" else gating_dt
      })
      
      
      # ----- datatable event handlers -----
      
      # Update back-end gating data when user edits name of cell subset in front-end datatable
      observeEvent(input$gating_pg_table_bg_cell_edit, {
        cell_edit_data <- input$gating_pg_table_bg_cell_edit
        gate_id <- reactive_gating_df_bg()$Gate_ID[cell_edit_data$row]
        gate_reactive_values_bg[[gate_id]]@name_subset_cells <- cell_edit_data$value
      })
      
      # when a datatable row is selected, set the selected_gate_bg reactive value to the gate that was selected in the datatable
      # this ensures that if the user wants to re-gate based on this gate, then the corresponding input gate stats are correct
      observeEvent(input$gating_pg_table_bg_rows_selected, {
        row_index <- input$gating_pg_table_bg_rows_selected
        selected_gate_bg(reactive_gating_df_bg()$Gate_ID[row_index])
        
        #update assay, x and y axis dropdowns to reflect the same axes shown in selected gate, so that if user gates from this selected gate, the right axes are recorded for the new gate
        updateSelectInput(
          session = session,
          inputId = "Assay_bg",
          selected = gate_list_bg()[[selected_gate_bg()]]@assay_name
        )
        updateSelectInput(
          session = session,
          inputId = "x_feature_bg",
          selected = gate_list_bg()[[selected_gate_bg()]]@x_axis
        )
        updateSelectInput(
          session = session,
          inputId = "y_feature_bg",
          selected = gate_list_bg()[[selected_gate_bg()]]@y_axis
        )
      })
      
      # when clear all gates button is clicked or new rds file is uploaded, reset gating info
      observeEvent(
        list(
          #list of input events that can trigger resetting of gating info
          input$clear_all_gates_bg,
          input$rds_input_file
        ), 
        {
          counter_reactive_bg(as.integer(0))
          
          # reset all values in gate_reactive_values to NULL
          if (!rlang::is_empty(names(gate_reactive_values_bg))) {
            gate_reactive_values_bg <- lapply(names(gate_reactive_values_bg), set_gates_to_null, local_gate_reactive_values = gate_reactive_values_bg)
          }
          
          #update last-buttons tracker
          last_buttons_clicked_bg$second_to_last <- last_buttons_clicked_bg$last
          last_buttons_clicked_bg$last <- "clear_all_gates_button"
        })
      
      
      # ----- backgating data download handlers -----
      output$download_as_list_rds_bg <- downloadHandler(
        filename = "backgate_info_list.rds",
        content = function(file) {
          saveRDS(gate_list_bg(), file = file)
        }
      )
      output$download_as_df_rds_bg <- downloadHandler(
        filename = "backgate_info_df.rds",
        content = function(file) {
          saveRDS(reactive_gating_df_bg(), file = file)
        }
      )
      
    })  # belongs to OBSERVE WRAPPER for backgating tab
  })  # belongs to end of observe wrapper FOR ALL BACKEND 
}  # end of server/back end code 
