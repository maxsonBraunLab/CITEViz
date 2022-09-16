#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' 
#' @import shiny
#' @import ggplot2
#' @import magrittr
#'  
#' @importFrom arrow read_feather
#' @importFrom dplyr filter last left_join
#' @importFrom plotly add_histogram2dcontour add_markers config event_data event_register ggplotly layout plot_ly renderPlotly
#' @importFrom rlang is_empty
#' @importFrom SeuratObject Assays DefaultAssay Embeddings FetchData GetAssayData Reductions
#' @importFrom stats quantile
#' @importFrom tools file_ext
#'
#' @noRd
app_server <- function( input, output, session ) {

  #set max file upload size to 3gb (default is only 5mb) since rds files can be really big
  #takes in an integer argument for max filesize in megabytes
  #to improve: set max file upload size based on user's hardware limitations?
  options(shiny.maxRequestSize = 3000 * 1024^2)
  
  # Your application server logic
  
  # ---------- Filehandling ---------- 
  # Putting everything in an observe function will put everything in the server function into the same environment allowing for
  # a single read of the uploaded seurat object instead of a read every time myso is called in render* function.
  # This speeds up the code immensely
  # The main function of this initial observe is to allow for a single upload of a Seurat Object over all pages.
  # Filetype validation: 
  # Note that only RDS, Arrow, or Feather files can be inputted by the user due to the UI fileInput() argument `accept = c(".rds", ".arrow", ".feather")`. 
  # This works on an actual web browser but not in the RStudio viewer.
  
  observe({
    
    input_file_df <- input$file_input
    if (is.null(input_file_df)) { return(NULL) }
    
    # set boolean flag for whether user's RDS file was already converted to arrow files prior to app startup
    # this determines whether the server reads Seurat data from a Seurat object or from an arrow file for plots, tables, etc.
    arrow_flag <- reactiveVal() # can be TRUE or FALSE
    arrow_filename_list <- reactiveVal(NULL)
    myso <- reactiveVal(NULL) # initialize seurat object here so it is accessible to everything else in server side
    valid_file_input_flag <- reactiveVal() #set this flag so if invalid files are uploaded, the rest of the app doesn't render and throw errors due to invalid file input
    
    observeEvent(
      list(
        #list of input events that can trigger reactive flag
        input$file_input
        # include app startup? (for reading in seurat object from memory instead of file upload)
      ),
      {
        # code to execute when one of the above input events occurs

        # check if original names of uploaded input files have RDS extension or not (i.e., are feather or arrow files)
        file_extensions <- tolower(tools::file_ext(input_file_df$name))
        
        #if only 1 file is uploaded by user
        if (length(file_extensions) == 1) {
          # if only 1 rds file is uploaded
          if (file_extensions == "rds") {
            arrow_flag(FALSE)
            message("Arrow flag is FALSE")
            
            # reset arrow_filename_list to NULL so that there are no issues with getting dropdown menu choices fxn later
            arrow_filename_list(NULL)
            
            #read in RDS file
            rds_obj <- readRDS(input_file_df$datapath)
            
            if (typeof(rds_obj) == "list") {
              #check if integrated obj exists before retrieving it from RDS that was read in
              integrated_obj_index <- grep("integrated", names(rds_obj), ignore.case = TRUE)
              
              #if integrated obj is not in list of objs and each seurat obj in the sample list is also wrapped in a list
              if (length(integrated_obj_index) == 0) {
                # myso <- rds_obj[[1]]
                myso(rds_obj[[1]])
              }
              else {
                #read in integrated obj
                # myso <- rds_obj[[integrated_obj_index]]
                myso(rds_obj[[integrated_obj_index]])
              }
            }
            else {
              #if integrated obj is not in list of objs, and the obj from the RDS file is just a single Seurat obj and not a list
              # myso <- rds_obj
              myso(rds_obj)
            }
            # set valid_file_input_flag after reading in RDS so that a Seurat object is in memory before other parts of app that require a true valid_file_input_flag can run (that way a "true" flag doesn't prematurely trigger other events to happen before a valid seurat obj is read in)
            valid_file_input_flag(TRUE)
          }
          # else if the 1 file uploaded is not RDS 
          else {
            valid_file_input_flag(FALSE)
            output$file_validation_status <- renderText({ 
              "Please upload either a single RDS file or multiple Feather files."
            })
          }
        }
        # else if multiple files are uploaded by user
        else {
          if ("rds" %in% file_extensions) {
            # send error message to user saying they need to upload either 1 rds file or multiple arrow/feather files
            valid_file_input_flag(FALSE)
            output$file_validation_status <- renderText({ 
              "Please upload either a single RDS file or multiple Feather files."
            })
          }
          else if (("feather" %in% file_extensions) | ("arrow" %in% file_extensions)) {
            arrow_flag(TRUE)
            message("Arrow flag is TRUE")
            myso(NULL) # reset myso to NULL so that there are no issues with getting dropdown menu choices fxn later
          
            arrow_filename_list(input_file_df$name)
            
            # set valid_file_input_flag after reading in file list so that items are in memory before other parts of app that require a true valid_file_input_flag can run (that way a "true" flag doesn't prematurely trigger other events to happen before a valid file list is read in)
            valid_file_input_flag(TRUE)
          }
          else {
            # send error message to user saying they need to upload either 1 rds file or multiple arrow/feather files
            valid_file_input_flag(FALSE)
            output$file_validation_status <- renderText({ 
              "Please upload either a single RDS file or multiple Feather files."
            })
          }
        }
        
       
        
        
        # if (file_extensions == "rds") {
          # arrow_flag(FALSE)
          # #read in RDS file
          # rds_obj <- readRDS(input_file_df$datapath)
          # 
          # if (typeof(rds_obj) == "list") {
          #   #check if integrated obj exists before retrieving it from RDS that was read in
          #   integrated_obj_index <- grep("integrated", names(rds_obj), ignore.case = TRUE)
          #   
          #   #if integrated obj is not in list of objs and each seurat obj in the sample list is also wrapped in a list
          #   if (length(integrated_obj_index) == 0) {
          #     myso <- rds_obj[[1]]
          #   }
          #   else {
          #     #read in integrated obj
          #     myso <- rds_obj[[integrated_obj_index]]
          #   }
          # }
          # else {
          #   #if integrated obj is not in list of objs, and the obj from the RDS file is just a single Seurat obj and not a list
          #   myso <- rds_obj
          # }
          
          #check if Seurat object has "reductions" slot
          # file_validation <- reactive({ 
          #   validate(
          #     need(
          #       length(SeuratObject::Reductions(myso)) > 0,
          #       message = "Seurat object does not contain reductions. Please check input RDS file."
          #     )
          #   )
          # })

          # output$file_validation_status <- renderText({ 
          #   # file_validation() 
          #   })
          
          # message("Arrow flag is FALSE")        
          # }
        # else if input file extension is not RDS (i.e. is arrow or feather), set arrow_flag to true
        # else {
        #   arrow_flag(TRUE)
        #   message("Arrow flag is TRUE")
        #   
        #   #get list of input files to read and refer to downstream
        #   
        # }
        
        # if not RDS, set arrow_flag(TRUE)
        # else set arrow_flag(FALSE) and read RDS
        
        # triggered by file input button in UI
        #   ## give it option to accept multiple files (either one RDS or multiple feather files)
        #   ## if the user inputs multiple RDS files, throw an error to the user to tell them not to do that
        # what happens if the user inputs a mix of RDS and Feather files at the same time?
      }) 

    
    # ---------- ***** QA ***** ---------- 
    # QA plots generated from Maxson-Braun lab's CITE-seq data preprocessing pipeline
    
    observe({
      
      # require valid_file_input_flag to be TRUE in order to run rest of section in observe wrapper so that app doesn't crash if invalid file(s) are inputted
      req(valid_file_input_flag() == TRUE)
      
      updateSelectInput(
        session = session,
        inputId = "color_qa",
        choices = get_choices("metadata",
                              arrow_flag(),
                              myso(),
                              arrow_filename_list(),
                              input_file_df),
        selected = get_choices("metadata", 
                               arrow_flag(), 
                               myso(), 
                               arrow_filename_list(), 
                               input_file_df)[1]
        # choices = colnames(myso()[[]][lapply(myso()[[]], class) %in% c("factor", "character")]),
        # selected = colnames(myso()[[]][lapply(myso()[[]], class) %in% c("factor", "character")])[1]
      )
      
      # ----- QA distribution plot -----
      # reactive distribution plot
      #reactive function will rerun this expression every time distribution_plot is called, which should be only when QA or color_qa choice is changed
      distribution_plot <- reactive({
        
        req(input$file_input, input$QA, input$color_qa, valid_file_input_flag() == TRUE)
        color <- input$color_qa
        
        #This assigns the params variable a list of strings that act a varying parameters depending on input of QA
        params <- switch(input$QA,
                         "RNA Count Per Cell" = c("nCount_RNA", "Distribution of Counts per Cell", "Number of Counts"),
                         "Gene Count Per Cell" = c("nFeature_RNA", "Distribution of Genes Detected per Cell", "Number of Unique Genes"),
                         "Percent Mitochondria" = c("percentMito", "Distribution of Mito GE per Cell", "Mitochondrial Ratio"),
                         "ADT Count Per Cell" = c("nCount_ADT", "ADT Counts per Cell","Number of Counts"),
                         "Unique ADTs Per Cell" = c("nFeature_ADT", "Distribution of CITE-seq Antibodies per Cell", "Number of Unique Antibodies")
        )
        
        metadata_df <- get_data(category = "metadata",
                                arrow_flag = arrow_flag(), 
                                seurat_object = myso(), 
                                arrow_file_list = arrow_filename_list(), 
                                input_file_df = input_file_df, 
                                assay_name = NULL, 
                                reduction_name = NULL)
        
        #creates list of x-coordinates for quantiles of data. 
        quant <- stats::quantile(x = metadata_df[,params[1]],
                          probs = c(0.5,0.75,0.95),
                          na.rm = TRUE)
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(metadata_df[[color]])))
        
        #generate base plot template with features that all QA distribution plots will have
        base_distrib_plot <- metadata_df %>%
          ggplot2::ggplot(aes(x = eval(parse(text = params[1])), fill = eval(parse(text = color)), color = eval(parse(text = color)))) + #eval(parse(text = x)) necessary to turn string into variable name format
          ggplot2::labs(fill = color, color = color) + 
          ggplot2::theme(plot.title = element_text(hjust=0.5)) +
          ggplot2::ggtitle(params[2]) +
          ggplot2::xlab(params[3]) +
          ggplot2::geom_vline(xintercept = quant, size = 0.5, alpha = 0.5, linetype = "dashed", color = "grey30") +
          ggplot2::scale_color_manual(values = custom_palette) +
          ggplot2::scale_fill_manual(values = custom_palette)
        
        #initialize QA distribution plot before if/else statements below so that the plot object can be accessed outside of the if/else statements
        final_distrib_plot <- base_distrib_plot
        
        #create density/bar plot for selected input. If integrated object is uploaded, then the original identity of the cells will separate into graphs per sample
        if (input$QA %in% "ADT Count Per Cell") {
          final_distrib_plot <- base_distrib_plot + ggplot2::scale_x_log10() + ggplot2::geom_density(alpha = 0.25) + ggplot2::ylab("Density")
        } 
        else if (input$QA %in% "Unique ADTs Per Cell") {
          final_distrib_plot <- base_distrib_plot + ggplot2::geom_bar(alpha = 0.5, position = "dodge") + ggplot2::ylab("Frequency")
        }
        else {
          final_distrib_plot <- base_distrib_plot + ggplot2::geom_density(alpha = 0.25) + ggplot2::ylab("Density")
        }
        
        #show distribution plot
        final_distrib_plot %>% 
          plotly::ggplotly() %>% 
          plotly::config(toImageButtonOptions = list(format = "png", scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          plotly::layout(title = list(font = list(size = 14)), hovermode = FALSE) 
      })
      
      
      # ----- QA box plot -----
      # reactive box plot 
      box_plot <- reactive({
        
        req(input$file_input, input$QA, input$color_qa, valid_file_input_flag() == TRUE)
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
        
        metadata_df <- get_data(category = "metadata",
                                arrow_flag = arrow_flag(), 
                                seurat_object = myso(), 
                                arrow_file_list = arrow_filename_list(), 
                                input_file_df = input_file_df, 
                                assay_name = NULL, 
                                reduction_name = NULL)
        
        #creates list of x-coordinates for quantiles of data. 
        quant <- stats::quantile(x = metadata_df[,params[1]],
                          probs = c(0.5,0.75,0.95),
                          na.rm = TRUE)
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(metadata_df[[color]])))
        
        #generate base plot template with features that all QA boxplots will have
        base_box_plot <- metadata_df %>%
          ggplot2::ggplot(aes(x = eval(parse(text = color)), y = eval(parse(text = params[1])), fill = eval(parse(text = color)), color = eval(parse(text = color)))) + #eval(parse(text = x)) necessary to turn string into variable name format
          ggplot2::labs(fill = color, color = color) + 
          ggplot2::geom_boxplot(alpha = 0.5, width=0.5) + 
          ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
          ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
          ggplot2::ggtitle(params[2]) +
          ggplot2::xlab("Sample") +
          ggplot2::ylab(params[3]) +
          ggplot2::geom_violin(alpha = 0.2) +
          ggplot2::geom_hline(yintercept = quant, size = 0.5, alpha = 0.5, linetype = "dashed", color = "grey30") +
          ggplot2::scale_color_manual(values = custom_palette) +
          ggplot2::scale_fill_manual(values = custom_palette)
        
        #initialize QA box plot before if/else statements below so that the plot object can be accessed outside of the if/else statements
        final_box_plot <- base_box_plot
        
        #create box plot for selected input. If integrated object is uploaded, then the original identity of the cells will separate into boxes per sample
        if (input$QA %in% "ADT Count Per Cell") {
          final_box_plot <- base_box_plot + ggplot2::scale_y_log10()
        } 
        
        #show box plot
        final_box_plot %>% 
          plotly::ggplotly() %>% 
          plotly::config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          plotly::layout(title = list(font = list(size = 14)), hovermode = FALSE)
      })
      
      # ----- render QA plots -----
      # render the reactive plotly QA plots 
      output$distrib_plot <- plotly::renderPlotly({ distribution_plot() })
      output$box_plot <- plotly::renderPlotly({ box_plot() })
      
    }) #end of QA observe() wrapper
    
    
    # ---------- ***** Clustering ***** ---------- 
    observe({
      
      # require valid_file_input_flag to be TRUE in order to run rest of section in observe wrapper so that app doesn't crash if invalid file(s) are inputted
      req(valid_file_input_flag() == TRUE)
      
      #changes the selectInput "reduction" dropdown contents to include all reductions in Seurat Object
      updateSelectInput(
        session = session,
        inputId = "reduction",
        choices = get_choices("reductions",
                              arrow_flag(),
                              myso(),
                              arrow_filename_list(),
                              input_file_df),
        selected = dplyr::last(get_choices("reductions", 
                                           arrow_flag(), 
                                           myso(), 
                                           arrow_filename_list(), 
                                           input_file_df))
        # choices = sort(SeuratObject::Reductions(myso())),
        # selected = dplyr::last(sort(SeuratObject::Reductions(myso())))
      )
      
      #changes selectInput "color" dropdown contents to include all metadata in Seurat Object
      updateSelectInput(
        session = session,
        inputId = "color1",
        choices = get_choices("metadata",
                              arrow_flag(),
                              myso(),
                              arrow_filename_list(),
                              input_file_df),
        selected = get_choices("metadata", 
                               arrow_flag(), 
                               myso(), 
                               arrow_filename_list(), 
                               input_file_df)[1]
        # choices = colnames(myso()[[]][lapply(myso()[[]], class) %in% c("factor", "character")]),
        # selected = colnames(myso()[[]][lapply(myso()[[]], class) %in% c("factor", "character")])[1]
      )
      
      # ----- Reactive 2D reduction graph -----
      reduc_plot <- reactive({
        req(input$file_input, input$reduction, input$color1, valid_file_input_flag() == TRUE)
        
        #create string for reduction to plot
        reduc <- input$reduction
        
        #selected metadata to color clusters by
        color <- input$color1
        
        metadata_df <- get_data(category = "metadata",
                                arrow_flag = arrow_flag(), 
                                seurat_object = myso(), 
                                arrow_file_list = arrow_filename_list(), 
                                input_file_df = input_file_df, 
                                assay_name = NULL, 
                                reduction_name = NULL)
       
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(metadata_df[[color]])))
        
        #create dataframe from reduction selected
        cell_data <- get_data(category = "reductions",
                              arrow_flag = arrow_flag(), 
                              seurat_object = myso(), 
                              arrow_file_list = arrow_filename_list(), 
                              input_file_df = input_file_df, 
                              assay_name = NULL, 
                              reduction_name = reduc)
        
        # cell_data <- data.frame(eval(parse(text = paste0("SeuratObject::Embeddings(object = myso(), reduction = '", reduc, "')"))))
        
        #create list containing all column names of cell_data
        cell_col <- colnames(cell_data)
        
        #show plot
        plotly::plot_ly(cell_data, 
                        x = ~cell_data[,1], y = ~cell_data[,2],
                        customdata = rownames(cell_data), #customdata is printed in cell selection and used to find metadata
                        color = stats::as.formula(paste0("~metadata_df$", color)),
                        # color = ~eval(parse(text = paste0("myso()[[]]$", color))), #color by selected metadata in object; need to incorporate 
                        colors = custom_palette,
                        type = "scatter", 
                        mode = "markers",
                        marker = list(size = 3, width = 2),
                        source = "A") %>%
          
          plotly::config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          #Layout changes the aesthetic of the plot
          plotly::layout(
            title = toupper(reduc),
            xaxis = list(title = cell_col[1]),
            yaxis = list(title = cell_col[2]),
            dragmode = "select") %>% 
          #Determines the mode of drag interactions. "select" and "lasso" apply only to scatter traces with markers or text. "orbit" and "turntable" apply only to 3D scenes.
          
          plotly::event_register("plotly_selected")
      })
      
      # ----- Reactive 3D reduction graph -----
      reduc_plot_3d <- reactive({
        req(input$file_input, input$reduction, input$color1, valid_file_input_flag() == TRUE)
        
        #create string for reduction to plot
        reduc <- input$reduction
        
        #selected metadata to color clusters by
        color <- input$color1
        
        metadata_df <- get_data(category = "metadata",
                                arrow_flag = arrow_flag(), 
                                seurat_object = myso(), 
                                arrow_file_list = arrow_filename_list(), 
                                input_file_df = input_file_df, 
                                assay_name = NULL, 
                                reduction_name = NULL)
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(metadata_df[[color]])))
        
        #create dataframe from reduction selected
        cell_data <- get_data(category = "reductions",
                              arrow_flag = arrow_flag(), 
                              seurat_object = myso(), 
                              arrow_file_list = arrow_filename_list(), 
                              input_file_df = input_file_df, 
                              assay_name = NULL, 
                              reduction_name = reduc)
        # cell_data <- data.frame(eval(parse(text = paste0("SeuratObject::Embeddings(object = myso(), reduction = '", reduc, "')"))))
        
        #create list containing all column names of cell_data
        cell_col <- colnames(cell_data)
        
        #show plot
        plotly::plot_ly(cell_data, 
                        x = ~cell_data[,1], y = ~cell_data[,2], z = ~cell_data[,3],
                        customdata = rownames(cell_data), #customdata is printed in cell selection and used to find metadata
                        color = stats::as.formula(paste0("~metadata_df$", color)),
                        # color = ~eval(parse(text = paste0("myso()[[]]$", color))), #color by selected metadata in object; need to incorporate 
                        colors =  custom_palette,
                        type = "scatter3d", 
                        mode = "markers",
                        marker = list(size = 2, width = 1)) %>%
          
          plotly::config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          #Layout changes the aesthetic of the plot
          plotly::layout(
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
        metadata_df <- get_data(category = "metadata",
                                arrow_flag = arrow_flag(), 
                                seurat_object = myso(), 
                                arrow_file_list = arrow_filename_list(), 
                                input_file_df = input_file_df, 
                                assay_name = NULL, 
                                reduction_name = NULL)
        
        selected_metadata_df <- metadata_df[plotly::event_data("plotly_selected", source = "A")$customdata, ]
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
    
    
    # ------- ***** Expression (1D) ***** ----------
    observe({
      
      # require valid_file_input_flag to be TRUE in order to run rest of section in observe wrapper so that app doesn't crash if invalid file(s) are inputted
      req(valid_file_input_flag() == TRUE)
      
      # ----- update/render UI elements -----
      updateSelectInput(
        session = session,
        inputId = "reduction_expr_1d",
        choices = get_choices("reductions",
                              arrow_flag(),
                              myso(),
                              arrow_filename_list(),
                              input_file_df),
        selected = dplyr::last(get_choices("reductions", 
                                           arrow_flag(), 
                                           myso(), 
                                           arrow_filename_list(), 
                                           input_file_df))
        # choices = sort(SeuratObject::Reductions(myso())),
        # selected = dplyr::last(sort(SeuratObject::Reductions(myso())))
      )
      
      output$Assay_1d <- renderUI({
        menu_choices <- get_choices("assays",
                                    arrow_flag(),
                                    myso(),
                                    arrow_filename_list(),
                                    input_file_df)
        selectInput(inputId = "Assay_1d",
                    label = "Choose assay to color reduction plot by:",
                    choices = menu_choices,
                    selected = menu_choices[1]
                    # choices = sort(SeuratObject::Assays(object = myso())),
                    # selected = sort(SeuratObject::Assays(object = myso()))[1]
        )
      })
      
      output$feature_1d <- renderUI({
        req(input$Assay_1d)
        # feature_path <- paste0('SeuratObject::GetAssayData(object = myso(), slot = "data", assay = "', input$Assay_1d, '")')
        assay_name <- input$Assay_1d
        menu_choices <- get_choices(category = NULL, 
                                    arrow_flag(), 
                                    myso(), 
                                    arrow_filename_list(), 
                                    input_file_df, 
                                    assay_name)
        selectInput(inputId = "feature_1d",
                    label = "Choose feature to view expression levels for:",
                    choices = menu_choices,
                    selected = menu_choices[1]
                    # choices = rownames(eval(parse(text=feature_path))),
                    # selected = rownames(eval(parse(text=feature_path)))[1]
        )
      })
      
      # ----- 1D gene/ADT expression reactive reduction graph -----
      expr_reduc_plot_1d <- eventReactive(
        list(
          #list of input events that can trigger reactive 
          input$file_input,
          input$reduction_expr_1d,
          input$feature_1d
        ), 
        {
          req(input$file_input, input$reduction_expr_1d, input$Assay_1d, input$feature_1d, valid_file_input_flag() == TRUE)
          
          #create string for reduction to plot
          reduc <- input$reduction_expr_1d
          
          #selected feature to color clusters by
          color_x <- input$feature_1d
          
          # temp_myso <- myso()
          # SeuratObject::DefaultAssay(temp_myso) <- input$Assay_1d
          # # SeuratObject::DefaultAssay(myso()) <- input$Assay_1d
          # myso(temp_myso)
          # 
          # count_data <- SeuratObject::FetchData(object = myso(), vars = color_x, slot = "data")
          
          count_data <- get_data(category = "assays",
                                 arrow_flag = arrow_flag(), 
                                 seurat_object = myso(), 
                                 arrow_file_list = arrow_filename_list(), 
                                 input_file_df = input_file_df, 
                                 assay_name = input$Assay_1d, 
                                 reduction_name = NULL,
                                 assay_data_to_get = input$color_x)
          
          #create dataframe from reduction selected
          cell_data <- get_data(category = "reductions",
                                arrow_flag = arrow_flag(), 
                                seurat_object = myso(), 
                                arrow_file_list = arrow_filename_list(), 
                                input_file_df = input_file_df, 
                                assay_name = NULL, 
                                reduction_name = reduc)
          # cell_data <- data.frame(eval(parse(text = paste0("SeuratObject::Embeddings(object = myso(), reduction = '", reduc, "')"))))
          
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
                                        reversescale = TRUE),
                          source = "expression_1d_plot") %>%
            
            plotly::config(toImageButtonOptions = list(format = "png",
                                                      scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
            ) %>%
            #Layout changes the aesthetic of the plot
            plotly::layout(
              title = toupper(reduc),
              xaxis = list(title = cell_col[1]),
              yaxis = list(title = cell_col[2]),
              dragmode = "select") %>% #Determines the mode of drag interactions. "select" and "lasso" apply only to scatter traces with markers or text. "orbit" and "turntable" apply only to 3D scenes.
        
            plotly::event_register("plotly_selected")
          })
      
      
      # ----- render reactive reduction plots -----
      output$exploration_reduct_1d <- renderPlotly({
        expr_reduc_plot_1d()
      })
      
      
      # ----- datatable of expression for cells selected in plotly -----
      # `server = FALSE` helps make it so that user can copy entire datatable to clipboard, not just the rows that are currently visible on screen
      output$expression_pg_selected <- DT::renderDT(server = FALSE, {
        req(input$file_input, input$Assay_1d, input$feature_1d, valid_file_input_flag() == TRUE)
        
        #selected feature to color clusters by
        color_x <- input$feature_1d
        
        # SeuratObject::DefaultAssay(myso()) <- input$Assay_1d
        # temp_myso <- myso()
        # SeuratObject::DefaultAssay(temp_myso) <- input$Assay_1d
        # myso(temp_myso)
        # count_data <- SeuratObject::FetchData(object = myso(), vars = color_x, slot = "data")
        
        count_data <- get_data(category = "assays",
                               arrow_flag = arrow_flag(), 
                               seurat_object = myso(), 
                               arrow_file_list = arrow_filename_list(), 
                               input_file_df = input_file_df, 
                               assay_name = input$Assay_1d, 
                               reduction_name = NULL,
                               assay_data_to_get = color_x)

        # num_cells_selected <- nrow(count_data)
        num_cells_expressing <- count_data %>%
          dplyr::filter(eval(parse(text = paste0("`", color_x, "`"))) > 0) %>%
          nrow()
        
        # get total num of cells in sample
        num_cells_total <- nrow(count_data)
        
        #convert count_data for selected cells into a dataframe
        selected_counts_df <- data.frame(Feature = color_x,
                                         Num_Cells_Expressing = num_cells_expressing,
                                         # Percent_of_Selected = 100 * (num_cells_expressing_subset / num_cells_selected),
                                         Percent_of_Total_Sample = 100 * num_cells_expressing / num_cells_total)

        selected_counts_dt <- DT::datatable(selected_counts_df, 
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
        if (is.null(selected_counts_dt)) "Brushed points appear here (double-click to clear)" else selected_counts_dt
      })
      
    }) # belongs to OBSERVE WRAPPER for expression tab
    
    
    # ------- ***** Co-Expression ***** ----------
    
    ### choose assay for x and y axes and then display dropdowns
    
    observe({
      
      # require valid_file_input_flag to be TRUE in order to run rest of section in observe wrapper so that app doesn't crash if invalid file(s) are inputted
      req(valid_file_input_flag() == TRUE)
      
      # ----- update/render UI elements -----
      
      updateSelectInput(
        session = session,
        inputId = "reduction_expr_2d",
        choices = get_choices("reductions",
                              arrow_flag(),
                              myso(),
                              arrow_filename_list(),
                              input_file_df),
        selected = dplyr::last(get_choices("reductions", 
                                           arrow_flag(), 
                                           myso(), 
                                           arrow_filename_list(), 
                                           input_file_df))
        # choices = sort(SeuratObject::Reductions(myso())),
        # selected = dplyr::last(sort(SeuratObject::Reductions(myso())))
      )
      
      output$Assay_x_axis <- renderUI({
        menu_choices <- get_choices("assays",
                                    arrow_flag(),
                                    myso(),
                                    arrow_filename_list(),
                                    input_file_df)
        selectInput(inputId = "Assay_x_axis",
                    label = "Choose assay for x-axis colorscale:",
                    choices = menu_choices,
                    selected = menu_choices[1]
                    # choices = sort(SeuratObject::Assays(object = myso())),
                    # selected = sort(SeuratObject::Assays(object = myso()))[1]
        )
      })
      
      output$Assay_y_axis <- renderUI({
        menu_choices <- get_choices("assays",
                                    arrow_flag(),
                                    myso(),
                                    arrow_filename_list(),
                                    input_file_df)
        selectInput(inputId = "Assay_y_axis",
                    label = "Choose assay for y-axis colorscale:",
                    choices = menu_choices,
                    selected = menu_choices[1]
                    # choices = sort(SeuratObject::Assays(object = myso())),
                    # selected = sort(SeuratObject::Assays(object = myso()))[1]
        )
      })

      output$x_axis_feature <- renderUI({
        req(input$Assay_x_axis)
        assay_name <- input$Assay_x_axis
        menu_choices <- get_choices(category = NULL, 
                                    arrow_flag(), 
                                    myso(), 
                                    arrow_filename_list(), 
                                    input_file_df, 
                                    assay_name)
        selectInput(inputId = "x_axis_feature",
                    label = "Choose feature for x-axis colorscale:",
                    choices = menu_choices,
                    selected = menu_choices[1]
        )
        # feature_path <- paste0('SeuratObject::GetAssayData(object = myso(), slot = "data", assay = "', input$Assay_x_axis, '")')
        # selectInput(inputId = "x_axis_feature",
        #             label = "Choose feature for x-axis colorscale:",
        #             choices = rownames(eval(parse(text=feature_path))),
        #             selected = rownames(eval(parse(text=feature_path)))[1]
        # )
      })
      
      output$y_axis_feature <- renderUI({
        req(input$Assay_y_axis)
        assay_name <- input$Assay_y_axis
        menu_choices <- get_choices(category = NULL, 
                                    arrow_flag(), 
                                    myso(), 
                                    arrow_filename_list(), 
                                    input_file_df, 
                                    assay_name)
        selectInput(inputId = "y_axis_feature",
                    label = "Choose feature for y-axis colorscale:",
                    choices = menu_choices,
                    selected = menu_choices[2]
        )
        # feature_path <- paste0('SeuratObject::GetAssayData(object = myso(), slot = "data", assay = "', input$Assay_y_axis, '")')
        # selectInput(inputId = "y_axis_feature",
        #             label = "Choose feature for y-axis colorscale:",
        #             choices = rownames(eval(parse(text=feature_path))),
        #             selected = rownames(eval(parse(text=feature_path)))[2]
        # )
      })
      
      
      # ----- 2D gene/ADT coexpression reactive reduction graph -----
      
      expr_reduc_plot_2d <- eventReactive(
        list(
          #list of input events that can trigger reactive plot
          input$file_input,
          input$reduction_expr_2d,
          input$x_axis_feature,
          input$y_axis_feature
        ), 
        {
          req(input$file_input, input$reduction_expr_2d, input$x_axis_feature, input$y_axis_feature, valid_file_input_flag() == TRUE)
          
          #create string for reduction to plot
          reduc <- input$reduction_expr_2d
          
          #selected metadata to color clusters by
          color_x <- input$x_axis_feature
          color_y <- input$y_axis_feature
          
          # SeuratObject::DefaultAssay(myso()) <- input$Assay_x_axis
          # temp_myso <- myso()
          # SeuratObject::DefaultAssay(temp_myso) <- input$Assay_x_axis
          # myso(temp_myso)
          # count_data_x <- SeuratObject::FetchData(object = myso(), vars = color_x, slot = "data")
          
          count_data_x <- get_data(category = "assays",
                                   arrow_flag = arrow_flag(), 
                                   seurat_object = myso(), 
                                   arrow_file_list = arrow_filename_list(), 
                                   input_file_df = input_file_df, 
                                   assay_name = input$Assay_x_axis, 
                                   reduction_name = NULL,
                                   assay_data_to_get = color_x)
          # extract only the count values as a vector from the original count data dataframe
          count_data_x <- count_data_x[[color_x]]
          
          # SeuratObject::DefaultAssay(myso()) <- input$Assay_y_axis
          # temp_myso <- myso()
          # SeuratObject::DefaultAssay(temp_myso) <- input$Assay_y_axis
          # myso(temp_myso)
          # count_data_y <- SeuratObject::FetchData(object = myso(), vars = color_y, slot = "data")
          
          count_data_y <- get_data(category = "assays",
                                   arrow_flag = arrow_flag(), 
                                   seurat_object = myso(), 
                                   arrow_file_list = arrow_filename_list(), 
                                   input_file_df = input_file_df, 
                                   assay_name = input$Assay_y_axis, 
                                   reduction_name = NULL,
                                   assay_data_to_get = color_y)
          
          # extract only the count values as a vector from the original count data dataframe
          count_data_y <- count_data_y[[color_y]]
          
          #create dataframe from reduction selected
          cell_data <- get_data(category = "reductions",
                                arrow_flag = arrow_flag(), 
                                seurat_object = myso(), 
                                arrow_file_list = arrow_filename_list(), 
                                input_file_df = input_file_df, 
                                assay_name = NULL, 
                                reduction_name = reduc)
          # cell_data <- data.frame(eval(parse(text = paste0("SeuratObject::Embeddings(object = myso(), reduction = '", reduc, "')"))))
          
          #create list containing all column names of cell_data
          cell_col <- colnames(cell_data)
          
          # map gene expression values to 2d color grid
          ngrid <- 16
          color_matrix_df <- get_color_matrix_df(ngrid)
          
          # use range() instead of max() to account for negative count data (flawed input data?)???
          # and use count_data + abs(min(count_data)) instead of just count_data in the numerator so account for color mapping of negative counts???
          coexpression_df <- data.frame(x_value = round(ngrid * count_data_x / max(count_data_x)),
                                        y_value = round(ngrid * count_data_y / max(count_data_y)))
          coexpression_umap_df <- cbind(coexpression_df, cell_data) #combine umap reduction data with expression data
          mapped_df <- dplyr::left_join(coexpression_umap_df, color_matrix_df) # map hex color codes to interpolated gene expression values in merged data and create a new data frame
          
          
          # create UMAP that colors by expression levels
          plotly::plot_ly(mapped_df,
                          source = "expression_2d_plot",
                          x = ~cell_data[,1], y = ~cell_data[,2],
                          customdata = rownames(cell_data), #customdata is printed in cell selection and used to find metadata
                          type = "scatter",
                          mode = "markers",
                          marker = list(size = 3,
                                        color = ~mapped_df$hex_color_mix
                                        )) %>%
           
            plotly::config(toImageButtonOptions = list(format = "png",
                                               scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
            ) %>%
            #Layout changes the aesthetic of the plot
            plotly::layout(
              showlegend = FALSE,
              title = toupper(reduc),
              xaxis = list(title = cell_col[1]),
              yaxis = list(title = cell_col[2]),
              dragmode = "select") %>% #Determines the mode of drag interactions. "select" and "lasso" apply only to scatter traces with markers or text. "orbit" and "turntable" apply only to 3D scenes.
            
            plotly::event_register("plotly_selected")
      })
      
      
      # ----- render reactive reduction plots -----
      output$color_legend_2d <- renderPlot({ 
        req(valid_file_input_flag() == TRUE)
        create_2d_color_legend(input = input, 
                               arrow_flag = arrow_flag(),
                               seurat_object = myso(), 
                               arrow_file_list = arrow_filename_list(), 
                               input_file_df = input_file_df) 
        })
      output$exploration_reduct_2d <- renderPlotly({ expr_reduc_plot_2d() })
      
      
      # ----- datatable of expression for cells selected in plotly -----
      # `server = FALSE` helps make it so that user can copy entire datatable to clipboard, not just the rows that are currently visible on screen
      output$coexpression_pg_selected <- DT::renderDT(server = FALSE, {
        req(input$file_input, input$Assay_x_axis, input$Assay_y_axis, input$x_axis_feature, input$y_axis_feature, valid_file_input_flag() == TRUE)
        
        #selected metadata to color clusters by
        color_x <- input$x_axis_feature
        color_y <- input$y_axis_feature
        
        # SeuratObject::DefaultAssay(myso()) <- input$Assay_x_axis
        # temp_myso <- myso()
        # SeuratObject::DefaultAssay(temp_myso) <- input$Assay_x_axis
        # myso(temp_myso)
        # 
        # count_data_x <- SeuratObject::FetchData(object = myso(), vars = color_x, slot = "data")

        # SeuratObject::DefaultAssay(myso()) <- input$Assay_y_axis
        # temp_myso <- myso()
        # SeuratObject::DefaultAssay(temp_myso) <- input$Assay_y_axis
        # myso(temp_myso)
        # count_data_y <- SeuratObject::FetchData(object = myso(), vars = color_y, slot = "data")
        
        count_data_x <- get_data(category = "assays",
                                 arrow_flag = arrow_flag(), 
                                 seurat_object = myso(), 
                                 arrow_file_list = arrow_filename_list(), 
                                 input_file_df = input_file_df, 
                                 assay_name = input$Assay_x_axis, 
                                 reduction_name = NULL,
                                 assay_data_to_get = color_x)
        
        count_data_y <- get_data(category = "assays",
                                 arrow_flag = arrow_flag(), 
                                 seurat_object = myso(), 
                                 arrow_file_list = arrow_filename_list(), 
                                 input_file_df = input_file_df, 
                                 assay_name = input$Assay_y_axis, 
                                 reduction_name = NULL,
                                 assay_data_to_get = color_y)
       
        num_cells_expressing_x <- count_data_x %>%
          # dplyr::filter(eval(parse(text = paste0("`", color_x, "`"))) > 0) %>%
          dplyr::filter((!!color_x) > 0) %>%
          nrow()
          
        num_cells_expressing_y <- count_data_y %>%
          # dplyr::filter(eval(parse(text = paste0("`", color_y, "`"))) > 0) %>%
          dplyr::filter((!!color_y) > 0) %>%
          nrow()
        
        # get total num of cells in sample
        num_cells_total_x <- nrow(count_data_x)
        num_cells_total_y <- nrow(count_data_y)
        
        # convert count_data for selected cells into a dataframe
        selected_counts_df <- data.frame(Feature1 = color_x,
                                         Feature2 = color_y,
                                         Num_Cells_Expressing_Feat1 = num_cells_expressing_x,
                                         Num_Cells_Expressing_Feat2 = num_cells_expressing_y,
                                         Percent_of_Total_Sample_Feat1 = 100 * num_cells_expressing_x / num_cells_total_x,
                                         Percent_of_Total_Sample_Feat2 = 100 * num_cells_expressing_y / num_cells_total_y)

        selected_counts_dt <- DT::datatable(selected_counts_df, 
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
        if (is.null(selected_counts_dt)) "Brushed points appear here (double-click to clear)" else selected_counts_dt
      })
      
    })  # belongs to OBSERVE WRAPPER for co-expression tab
    
    
    # ---------- ***** Gating ***** ----------
    observe({
      
      # require valid_file_input_flag to be TRUE in order to run rest of section in observe wrapper so that app doesn't crash if invalid file(s) are inputted
      req(valid_file_input_flag() == TRUE)
      
      # ----- update/render UI elements -----
      output$Assay <- renderUI({
        menu_choices <- get_choices("assays",
                                    arrow_flag(),
                                    myso(),
                                    arrow_filename_list(),
                                    input_file_df)
        selectInput(inputId = "Assay",
                    label = "Choose assay:",
                    choices = menu_choices,
                    selected = menu_choices[1]
                    # choices = sort(SeuratObject::Assays(object = myso())),
                    # selected = sort(SeuratObject::Assays(object = myso()))[1])
        )
      })
      
      output$x_feature <- renderUI({
        req(input$Assay)
        assay_name <- input$Assay
        menu_choices <- get_choices(category = NULL, 
                                    arrow_flag(), 
                                    myso(), 
                                    arrow_filename_list(), 
                                    input_file_df, 
                                    assay_name)
        selectInput(inputId = "x_feature",
                    label = "Choose x-axis feature:",
                    choices = menu_choices,
                    selected = menu_choices[1]
        )
        # feature_path <- paste0('SeuratObject::GetAssayData(object = myso(), slot = "data", assay = "', input$Assay, '")')
        # selectInput(
        #   inputId = "x_feature",
        #   label = "Choose x-axis feature:",
        #   choices = rownames(eval(parse(text=feature_path))),
        #   selected = rownames(eval(parse(text=feature_path)))[1])
      })
      
      output$y_feature <- renderUI({
        req(input$Assay)
        assay_name <- input$Assay
        menu_choices <- get_choices(category = NULL, 
                                    arrow_flag(), 
                                    myso(), 
                                    arrow_filename_list(), 
                                    input_file_df, 
                                    assay_name)
        selectInput(inputId = "y_feature",
                    label = "Choose y-axis feature:",
                    choices = menu_choices,
                    selected = menu_choices[2]
        )
        # feature_path <- paste0('SeuratObject::GetAssayData(object = myso(), slot = "data", assay = "', input$Assay, '")')
        # selectInput(
        #   inputId = "y_feature",
        #   label = "Choose y-axis feature:",
        #   choices = rownames(eval(parse(text=feature_path))),
        #   selected = rownames(eval(parse(text=feature_path)))[2])
      })
      
      #changes the selectInput "reduction" dropdown contents to include all reductions in Seurat Object
      updateSelectInput(
        session = session,
        inputId = "reduction_g",
        choices = get_choices("reductions",
                              arrow_flag(),
                              myso(),
                              arrow_filename_list(),
                              input_file_df),
        selected = dplyr::last(get_choices("reductions", 
                                           arrow_flag(), 
                                           myso(), 
                                           arrow_filename_list(), 
                                           input_file_df))
        # choices = sort(SeuratObject::Reductions(myso())),
        # selected = dplyr::last(sort(SeuratObject::Reductions(myso())))
      )
      updateSelectInput(
        session = session,
        inputId = "color2",
        choices = get_choices("metadata",
                              arrow_flag(),
                              myso(),
                              arrow_filename_list(),
                              input_file_df),
        selected = get_choices("metadata", 
                               arrow_flag(), 
                               myso(), 
                               arrow_filename_list(), 
                               input_file_df)[1]
        # choices = colnames(myso()[[]][lapply(myso()[[]], class) %in% c("factor", "character")]),
        # selected = colnames(myso()[[]][lapply(myso()[[]], class) %in% c("factor", "character")])[1]
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
          input$file_input,
          input$gate,
          input$reset_adt_scatter,
          input$clear_all_gates,
          input$x_feature,
          input$y_feature,
          input$gating_pg_table_rows_selected
        ), 
        { 
          #code to execute when one of the above input events occurs
          req(input$x_feature, input$y_feature, valid_file_input_flag() == TRUE)
          
          # SeuratObject::DefaultAssay(myso()) <- input$Assay
          # temp_myso <- myso()
          # SeuratObject::DefaultAssay(temp_myso) <- input$Assay
          # myso(temp_myso)
          
          # count_data <- SeuratObject::FetchData(object = myso(), vars = c(input$x_feature, input$y_feature), slot = "data")
          
          count_data <- get_data(category = "assays",
                                arrow_flag = arrow_flag(), 
                                seurat_object = myso(), 
                                arrow_file_list = arrow_filename_list(), 
                                input_file_df = input_file_df, 
                                assay_name = input$Assay, 
                                reduction_name = NULL,
                                assay_data_to_get = c(input$x_feature, input$y_feature))
          
          #generate dataframe for custom colorscale for contour plot, where each hex color code is mapped to a specific z-value between 0 and 1 (inclusive)
          #colorscale needs to be in this format for Plotly's add_histogram2dcontour(colorscale = ...) parameter
          gating_color_scale <- data.frame(
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
              plotly::add_histogram2dcontour(showscale = FALSE, ncontours = 10, colorscale = gating_color_scale, 
                                     contours = list(coloring='heatmap')) %>%
              plotly::add_markers(x = count_data[,input$x_feature],
                          y = count_data[,input$y_feature],
                          marker = list(size=2),
                          color = I("black"),
                          alpha = 0.6) 
          }
          
          
          else {
            selected_cell_barcodes <- NULL
            count_data_subset <- NULL
            
            if (last_buttons_clicked$last == "gate_button" & is.null(input$gating_pg_table_rows_selected)) {
              selected_cell_barcodes <- plotly::event_data("plotly_selected", source = "C")$customdata
              count_data_subset <- count_data[rownames(count_data) %in% selected_cell_barcodes, ]
            }
            else if (!is.null(input$gating_pg_table_rows_selected)) {
              selected_cell_barcodes <- GetData(gate_list()[[selected_gate()]], "subset_cells")[[1]]
              count_data_subset <- count_data[rownames(count_data) %in% selected_cell_barcodes, ]
            }
            base_scatterplot <- plotly::plot_ly(count_data_subset,
                                                x = ~count_data_subset[,input$x_feature], y = ~count_data_subset[,input$y_feature],
                                                customdata = rownames(count_data_subset),
                                                mode = "markers",
                                                source = "C") %>% 
              plotly::add_histogram2dcontour(showscale=FALSE, ncontours=10, colorscale = gating_color_scale, 
                                     contours = list(coloring='heatmap')) %>%
              plotly::add_markers(x = count_data_subset[,input$x_feature], 
                          y = count_data_subset[,input$y_feature], 
                          marker = list(size=2.5), 
                          color = I("black"), 
                          alpha = 0.6)
          }
          
          # add configuration and layout options to base scatterplot, register selection events
          base_scatterplot %>%
            plotly::config(toImageButtonOptions = list(format = "png",
                                               scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
            ) %>%
            #Layout changes the aesthetic of the plot
            plotly::layout(
              title = "Normalized Feature Scatter Plot",
              xaxis = list(title = input$x_feature),
              yaxis = list(title = input$y_feature),
              showlegend = FALSE,
              dragmode = "select") %>% #Determines the mode of drag interactions. "select" and "lasso" apply only to scatter traces with markers or text. "orbit" and "turntable" apply only to 3D scenes.
            
            plotly::event_register("plotly_selected")
        })
      
      # ----- reactive gating 2D reduction graph -----
      gating_reduc_plot <- reactive({
        req(input$file_input, input$color2, input$reduction_g, valid_file_input_flag() == TRUE)
        
        #create string for reduction to plot
        reduc <- input$reduction_g
        
        #selected metadata to color clusters by
        color <- input$color2
        
        metadata_df <- get_data(category = "metadata",
                                arrow_flag = arrow_flag(), 
                                seurat_object = myso(), 
                                arrow_file_list = arrow_filename_list(), 
                                input_file_df = input_file_df, 
                                assay_name = NULL, 
                                reduction_name = NULL)
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(metadata_df[[color]])))
        
        #create dataframe from reduction selected
        cell_data <- get_data(category = "reductions",
                              arrow_flag = arrow_flag(), 
                              seurat_object = myso(), 
                              arrow_file_list = arrow_filename_list(), 
                              input_file_df = input_file_df, 
                              assay_name = NULL, 
                              reduction_name = reduc)
        # cell_data <- data.frame(eval(parse(text = paste0("SeuratObject::Embeddings(object = myso(), reduction = '", reduc, "')"))))
        
        #create list containing all column names of cell_data
        cell_col <- colnames(cell_data)
        
        # initialize selected_cells before if/else statements so it can be accessed outside of the statements
        selected_cells <- NULL
        
        if (!is.null(input$gating_pg_table_rows_selected)) {
          selected_cells <- GetData(gate_list()[[selected_gate()]], "subset_cells")[[1]]
        }
        else {
          selected_cells <- plotly::event_data("plotly_selected", source = "C")$customdata
        }
        
        if (is.null(selected_cells)) {
          plotly_color_list <- c(paste0("metadata_df$", color), custom_palette)
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
          
          plotly::config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          #Layout changes the aesthetic of the plot
          plotly::layout(
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
        
        # SeuratObject::DefaultAssay(myso()) <- input$Assay
        # temp_myso <- myso()
        # SeuratObject::DefaultAssay(temp_myso) <- input$Assay
        # myso(temp_myso)
        # 
        # count_data <- SeuratObject::FetchData(object = myso(), vars = c(input$x_feature, input$y_feature), slot = "data")
        
        count_data <- get_data(category = "assays",
                               arrow_flag = arrow_flag(), 
                               seurat_object = myso(), 
                               arrow_file_list = arrow_filename_list(), 
                               input_file_df = input_file_df, 
                               assay_name = input$Assay, 
                               reduction_name = NULL,
                               assay_data_to_get = c(input$x_feature, input$y_feature))
        
        # get plotly event data
        sel <- plotly::event_data("plotly_selected", source = "C")
        brushed_coords <- plotly::event_data("plotly_brushed", source = "C")

        # increment counter every time gate button is clicked
        counter <- as.integer(counter_reactive() + 1)
        counter_reactive(counter)

        # create gate object based on UI input
        gate_reactive_values[[paste0("gate_", counter)]] <- create_gate_from_input(input = input,
                                                                                   is_forward_gating = TRUE,
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
        gate_reactive_values[[gate_id]] <- SetSubsetName(gate_reactive_values[[gate_id]], cell_edit_data$value)
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
          selected = GetData(gate_list()[[selected_gate()]], "assay_name")
        )
        updateSelectInput(
          session = session,
          inputId = "x_feature",
          selected = GetData(gate_list()[[selected_gate()]], "x_axis")
        )
        updateSelectInput(
          session = session,
          inputId = "y_feature",
          selected = GetData(gate_list()[[selected_gate()]], "y_axis")
        )
      })
      
      # when clear all gates button is clicked or new rds file is uploaded, reset gating info
      observeEvent(
        list(
          #list of input events that can trigger resetting of gating info
          input$clear_all_gates,
          input$file_input
        ), 
        {
          # require valid_file_input_flag to be TRUE in order to run rest of section in observe wrapper so that app doesn't crash if invalid file(s) are inputted
          req(valid_file_input_flag() == TRUE)
          
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
      
      # require valid_file_input_flag to be TRUE in order to run rest of section in observe wrapper so that app doesn't crash if invalid file(s) are inputted
      req(valid_file_input_flag() == TRUE)
      
      # ----- update/render UI elements -----
      output$Assay_bg <- renderUI({
        menu_choices <- get_choices("assays",
                                    arrow_flag(),
                                    myso(),
                                    arrow_filename_list(),
                                    input_file_df)
        selectInput(inputId = "Assay_bg",
                    label = "Choose assay:",
                    choices = menu_choices,
                    selected = menu_choices[1]
                    # choices = sort(SeuratObject::Assays(object = myso())),
                    # selected = sort(SeuratObject::Assays(object = myso()))[1])
        )
      })
      
      output$x_feature_bg <- renderUI({
        req(input$Assay_bg)
        assay_name <- input$Assay_bg
        menu_choices <- get_choices(category = NULL, 
                                    arrow_flag(), 
                                    myso(), 
                                    arrow_filename_list(), 
                                    input_file_df, 
                                    assay_name)
        selectInput(inputId = "x_feature_bg",
                    label = "Choose x-axis feature:",
                    choices = menu_choices,
                    selected = menu_choices[1]
        )
        # feature_path <- paste0('SeuratObject::GetAssayData(object = myso(), slot = "data", assay = "', input$Assay_bg, '")')
        # selectInput(
        #   inputId = "x_feature_bg",
        #   label = "Choose x-axis feature:",
        #   choices = rownames(eval(parse(text=feature_path))),
        #   selected = rownames(eval(parse(text=feature_path)))[1])
      })
      
      output$y_feature_bg <- renderUI({
        req(input$Assay_bg)
        assay_name <- input$Assay_bg
        menu_choices <- get_choices(category = NULL, 
                                    arrow_flag(), 
                                    myso(), 
                                    arrow_filename_list(), 
                                    input_file_df, 
                                    assay_name)
        selectInput(inputId = "y_feature_bg",
                    label = "Choose y-axis feature:",
                    choices = menu_choices,
                    selected = menu_choices[2]
        )
        # feature_path <- paste0('SeuratObject::GetAssayData(object = myso(), slot = "data", assay = "', input$Assay_bg, '")')
        # selectInput(
        #   inputId = "y_feature_bg",
        #   label = "Choose y-axis feature:",
        #   choices = rownames(eval(parse(text=feature_path))),
        #   selected = rownames(eval(parse(text=feature_path)))[2])
      })
      
      #changes the selectInput "reduction" dropdown contents to include all reductions in Seurat Object
      updateSelectInput(
        session = session,
        inputId = "reduction_bg",
        choices = get_choices("reductions",
                              arrow_flag(),
                              myso(),
                              arrow_filename_list(),
                              input_file_df),
        selected = dplyr::last(get_choices("reductions", 
                                           arrow_flag(), 
                                           myso(),
                                           arrow_filename_list(), 
                                           input_file_df))
        # choices = sort(SeuratObject::Reductions(myso())),
        # selected = dplyr::last(sort(SeuratObject::Reductions(myso())))
      )
      updateSelectInput(
        session = session,
        inputId = "color2_bg",
        choices = get_choices("metadata",
                              arrow_flag(),
                              myso(),
                              arrow_filename_list(),
                              input_file_df),
        selected = get_choices("metadata", 
                               arrow_flag(), 
                               myso(), 
                               arrow_filename_list(), 
                               input_file_df)[1]
        # choices = colnames(myso()[[]][lapply(myso()[[]], class) %in% c("factor", "character")]),
        # selected = colnames(myso()[[]][lapply(myso()[[]], class) %in% c("factor", "character")])[1]
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
          input$file_input,
          input$gate_bg,
          input$reset_adt_scatter_bg,
          input$clear_all_gates_bg,
          input$x_feature_bg,
          input$y_feature_bg,
          input$gating_pg_table_bg_rows_selected,
          plotly::event_data("plotly_selected", source = "D")
        ), 
        {
          #code to execute when one of the above input events occurs
          req(input$x_feature_bg, input$y_feature_bg, valid_file_input_flag)
          

          gating_color_scale <- data.frame(
            z = c(0.0, 0.20, 0.40, 0.60, 0.80, 1.0),
            col = c("#FFFFFF", "#4564FE", "#76EFFF", "#FFF900", "#FFA300", "#FF1818")
          )
          
          # SeuratObject::DefaultAssay(myso()) <- input$Assay_bg
          # temp_myso <- myso()
          # SeuratObject::DefaultAssay(temp_myso) <- input$Assay_bg
          # myso(temp_myso)
          # 
          # count_data <- SeuratObject::FetchData(object = myso(), vars = c(input$x_feature_bg, input$y_feature_bg), slot = "data")
          
          count_data <- get_data(category = "assays",
                                 arrow_flag = arrow_flag(), 
                                 seurat_object = myso(), 
                                 arrow_file_list = arrow_filename_list(), 
                                 input_file_df = input_file_df, 
                                 assay_name = input$Assay_bg, 
                                 reduction_name = NULL,
                                 assay_data_to_get = c(input$x_feature_bg, input$y_feature_bg))
          
          
          selected_cell_barcodes <- NULL
          
          if (is.null(input$gating_pg_table_bg_rows_selected)) {
            selected_cell_barcodes <- plotly::event_data("plotly_selected", source = "D")$customdata
          }
          else {
            selected_cell_barcodes <- GetData(gate_list_bg()[[selected_gate_bg()]], "subset_cells")[[1]]
          }
          
          plotly::plot_ly(count_data,
                          x = ~count_data[,input$x_feature_bg], y = ~count_data[,input$y_feature_bg],
                          customdata = rownames(count_data),
                          mode = "markers",
                          color = rownames(count_data) %in% selected_cell_barcodes, #color cells by whether they're in the selection or not
                          colors = c("grey", "black")) %>% 
            plotly::add_histogram2dcontour(showscale = FALSE,
                                   ncontours = 10,
                                   colorscale = gating_color_scale,
                                   contours = list(coloring='heatmap'),
                                   color = I("black"),
                                   size = I(1.5)) %>%
            plotly::add_markers(x = count_data[,input$x_feature_bg],
                        y = count_data[,input$y_feature_bg],
                        marker = list(size=2),
                        alpha = 1) %>%
            plotly::config(toImageButtonOptions = list(format = "png", scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
            ) %>%
            #Layout changes the aesthetic of the plot
            plotly::layout(
              title = "Normalized Feature Scatter Plot",
              xaxis = list(title = input$x_feature_bg),
              yaxis = list(title = input$y_feature_bg),
              showlegend = FALSE)
        })
      
      # ----- backgate reactive reduction plot of selected cell features -----
      gating_reduc_plot_bg <- reactive({
        req(input$file_input, input$color2_bg, input$reduction_bg, valid_file_input_flag)
        
        #creates string for reduction to plot
        reduc <- input$reduction_bg
        
        #selected metadata to color clusters by
        color <- input$color2_bg
        
        metadata_df <- get_data(category = "metadata",
                                arrow_flag = arrow_flag(), 
                                seurat_object = myso(), 
                                arrow_file_list = arrow_filename_list(), 
                                input_file_df = input_file_df, 
                                assay_name = NULL, 
                                reduction_name = NULL)
        
        #interpolate the base color palette so that exact number of colors in custom palette is same as number of unique values for selected metadata category
        custom_palette <- get_palette(length(unique(metadata_df[[color]])))
        plotly_color_list <- c(paste0("metadata_df$", color), custom_palette)
        
        #creates dataframe from reduction selected
        cell_data <- get_data(category = "reductions",
                              arrow_flag = arrow_flag(), 
                              seurat_object = myso(), 
                              arrow_file_list = arrow_filename_list(), 
                              input_file_df = input_file_df, 
                              assay_name = NULL, 
                              reduction_name = reduc)
        # cell_data <- data.frame(eval(parse(text = paste0("SeuratObject::Embeddings(object = myso(), reduction = '", reduc, "')"))))
        
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
          plotly::config(toImageButtonOptions = list(format = "png",
                                             scale = 10) #scale title/legend/axis labels by this factor so that they are high-resolution when downloaded
          ) %>%
          #Layout changes the aesthetic of the plot
          plotly::layout(title = toupper(reduc),
                 xaxis = list(title = cell_col[1]),
                 yaxis = list(title = cell_col[2]),
                 dragmode = "select") %>%
          plotly::event_register("plotly_selected")
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
        
        # SeuratObject::DefaultAssay(myso()) <- input$Assay_bg
        # temp_myso <- myso()
        # SeuratObject::DefaultAssay(temp_myso) <- input$Assay_bg
        # myso(temp_myso)
        # 
        # count_data <- SeuratObject::FetchData(object = myso(), vars = c(input$x_feature_bg, input$y_feature_bg), slot = "data")
        
        count_data <- get_data(category = "assays",
                               arrow_flag = arrow_flag(), 
                               seurat_object = myso(), 
                               arrow_file_list = arrow_filename_list(), 
                               input_file_df = input_file_df, 
                               assay_name = input$Assay_bg, 
                               reduction_name = NULL,
                               assay_data_to_get = c(input$x_feature_bg, input$y_feature_bg))
        
        # get plotly event data
        sel <- plotly::event_data("plotly_selected", source = "D")
        brushed_coords <- plotly::event_data("plotly_brushed", source = "D")

        # increment counter every time gate button is clicked
        counter <- as.integer(counter_reactive_bg() + 1)
        counter_reactive_bg(counter)
        
        # create gate object based on UI input
        gate_reactive_values_bg[[paste0("gate_", counter)]] <- create_gate_from_input(input = input,
                                                                                      is_forward_gating = FALSE,
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
        gate_reactive_values_bg[[gate_id]] <- SetSubsetName(gate_reactive_values_bg[[gate_id]], cell_edit_data$value)
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
          selected = GetData(gate_list_bg()[[selected_gate_bg()]], "assay_name")
        )
        updateSelectInput(
          session = session,
          inputId = "x_feature_bg",
          selected = GetData(gate_list_bg()[[selected_gate_bg()]], "x_axis")
        )
        updateSelectInput(
          session = session,
          inputId = "y_feature_bg",
          selected = GetData(gate_list_bg()[[selected_gate_bg()]], "y_axis")
        )
      })
      
      # when clear all gates button is clicked or new rds file is uploaded, reset gating info
      observeEvent(
        list(
          #list of input events that can trigger resetting of gating info
          input$clear_all_gates_bg,
          input$file_input
        ), 
        {
          # require valid_file_input_flag to be TRUE in order to run rest of section in observe wrapper so that app doesn't crash if invalid file(s) are inputted
          req(valid_file_input_flag() == TRUE)
          
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
