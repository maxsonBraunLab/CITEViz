#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' 
#' @importFrom bslib bs_theme
#' @importFrom DT DTOutput
#' @importFrom plotly plotlyOutput
#' @importFrom vembedr embed_url
#' 
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    
    # Your application UI logic 
    navbarPage(
      id = "navbar_pg",
      title = div(icon("earlybirds", lib="font-awesome"), "CITE-Viz"),
      collapsible = TRUE,
      fluid = TRUE,
      theme = bslib::bs_theme(bootswatch = "flatly"),
      
      header = tags$header(
        div(
          # to remove
          #style = "margin-left: 0; padding-left: 1rem; padding-right: 1rem;",
          fileInput(
            inputId = "file_input", 
            placeholder = "Upload an RDS (.rds) file",
            accept = c(".rds"),
            label = NULL,
            multiple = FALSE)
        ),
        div(
          # to remove
          # style = "padding-left: 1rem; padding-right: 1rem;",
          textOutput(outputId = "file_validation_status")
        )
      ),
      
      footer = tags$footer(a(icon("github"), href="https://github.com/maxsonBraunLab/CITE-Viz"), "| Made with", icon("heart"), "at UO + OHSU", 
                           # to remove
                           # style = "
                           # background-image: linear-gradient(to right, #004F6E, #8D1D58);
                           # color: white;
                           # opacity: 0.97;
                           # font-size: 85%;
                           # width: 100%;
                           # text-align: center;
                           # padding: 0.4rem 0 0.4rem 0;
                           # position: fixed;
                           # bottom: 0;
                           # z-index: 1000;"
      ),
      
      
      # ---------- UI Landing/welcome tab ---------- 
      tabPanel("Getting Started", 
               # css style below can't easily be put into separate css file due to id of tab container being generated uniquely each time app is run
               style = "padding-left: 1.5rem; padding-right: 1.5rem; margin-left: auto; margin-right: auto; max-width: 1040px",
               
               #if user doesn't have javascript on their browser, display a message box near top of page
               tags$noscript(
                 div(
                   #to remove
                   # style = "
                   #     border: 2px solid gray; 
                   #     border-radius: 5px;
                   #     background-color: whitesmoke; 
                   #     margin-top: 0.1rem;
                   #     margin-bottom: 1.2rem;
                   #     padding: 1rem 1rem 0 1rem;",
                   
                   p("This app uses JavaScript. For the full user experience, please enable JavaScript in your browser.")
                 )
               ),
               
               h2(id = "welcome_h2",
                  strong("Welcome!"),
                  style = "margin-top: 0;"
               ),
               
               p("This app was developed through the", a("Knight Campus Graduate Internship Program", href= "https://internship.uoregon.edu/bioinformatics"), "at the University of Oregon, in collaboration with both the", a("Maxson", href="https://www.maxsonlab.org/"), "and", a("Braun", href="https://www.braunlab.org/team.html"), "Labs at Oregon Health and Science University (OHSU)  in Portland, Oregon."),
               p("The purpose of this app is to allow users to interactively engage with their CITE-seq data though interactive plotting, subsetting and flow cytometry-like cell gating features. Although this app was originally created to investigate acute myeloid leukemia (AML), one of the fastest-growing blood cancers, it can also be used to explore many other CITE-seq datasets."), 
               
               hr(style = "border-top: 3px solid whitesmoke;"),
               
               h3(strong("How to Use this App")), 
               p("To begin, upload an .rds file of a single Seurat object or a list of Seurat objects (with one of the objects being an 'integrated' object) in the file upload box at the top of any page."),
               
               h5(strong("Quality Assurance (QA)")),
               p("The Quality Assurance (QA) page allows the user to visualize data by several QA parameters including, but not limited to, mitochondrial ratio, RNA count per cell, and ADT count per cell, along different metadata separations. Dotted lines in each QA plot represent the values at which 50%, 75%, and 95% of the data falls at or below that value."),
               img(src = "www/QA_page.png", width = "100%"),
               
               br(),
               h5(strong("Clustering")),
               p("The Clustering page allows the user to select from a dropdown of dimensionality reduction plots and color cells by different metadata while viewing plots and cluster relationships in two- and three-dimensional space."),
               p("When the cursor hovers in the 2D reduction plot, the selection tools (box or lasso selection) appear for the user to choose from. The metadata for selected cells appears in the data table below the plots. The user also has the option to print or copy this data to clipboard."),
               img(src = "www/clustering_page.png", width = "100%"),
              
               br(),
               h5(strong("Gating")),
               
               p("The Gating page allows the user to select, gate, and subset cells by a variety of metadata parameters in the feature scatterplot. The user can input a name for a subset of cells to be gated by providing one in the text-input field of the sidebar panel prior to clicking the gate button. The user can also rename the selected cell subset by double-clicking in the Subset_Name cell of the data table, typing a new name, then pressing the Tab key or clicking outside of the table."),
               p("The user may also click on previous gates from the data table to view the feature scatter plot of that cell population. Using this feature, the user may also revert back to a previously created gate for further subsetting and gating."),
               p("To permanently delete all gating information in the data table, the user can click the 'Clear All Data' button. To view detailed gating data, the user can download this data in a variety of formats listed on the gating page."),
               p("To reset the colors of the cells in the reduction plot, double-click anywhere in the feature scatter plot."),
               
               vembedr::embed_url("https://youtu.be/JgaHlpWAKs4"),
               br(),
               
               h5(strong("Backgating")),
               p("The backgating feature works similarly to the Gating page, but allows the user to select a range of cells from the reduction plot (e.g. UMAP, PCA etc.) and locate them in the feature scatterplot for more extensive data exploration. Cell population selection and naming methods are identical to the gating page. To reset the colors of the cells in the feature scatter plot, double-click anywhere in the reduction plot."),
               vembedr::embed_url("https://youtu.be/KD8vH1l_Ju4"),
               
               hr(style = "border-top: 3px solid whitesmoke;"),
               
               h3(strong("Resources")),
               p("Consider using the", a("Maxson-Braun Lab's pipeline", href="https://github.com/maxsonBraunLab/cite_seq"), "to generate a Seurat object (saved as an .rds file) that can be uploaded into this app for further data analysis."), 
               p("If you are researching white blood cells, please see the following", a("Human Immune Cell Marker Guide", href="https://media.cellsignal.com/www/pdfs/science/pathways/Immune-Cell-Markers-Human.pdf"), "reproduced courtesy of", a("Cell Signaling Technology, Inc.", href="https://www.cellsignal.com")), 
               
               hr(style = "border-top: 3px solid whitesmoke;"),
               
               # uncomment this citation section and update the actual citation below once CITE-Viz has been published
               # h5(strong("How to Cite CITE-Viz")),
               # p("For citing this app in a publication, use:"),
               # tags$blockquote("Garth L. Kong, Thai T. Nguyen, Wesley K. Rosales, Anjali D. Panikar, John H. W. Cheney, Brittany M. Curtiss, Sarah A. Carratt, Theodore P. Braun, Julia E. Maxson. CITE-Viz: Replicating the Interactive Flow Cytometry Workflow in CITE-Seq.", style = "font-family: monospace;"),
               
               h5(strong("Acknowledgements")), 
               p("The development team would like to express sincere gratitude to this project's principal investigators, Dr. Julia Maxson and Dr. Ted Braun, and our project mentor, Garth Kong, for their guidance, encouragement and experience in the realm of blood cancers and CITE-seq data processing, not to mention their dedication to", a('"end cancer as we know it."', href ="https://ohsufoundation.org/stories/teaming-up-against-cancer-a-message-from-brian-druker-md/"), 'Their innumerable technical contributions helped bring this project to fruition.'),
               p("Furthermore, the development team would like to thank faculty and advisors of the Bioinformatics and Genomics Masters Program at the University of Oregon, expressly, Dr. Leslie Coonrod, Dr. Stacey Wagner, Pete Batzel and Jason Sydes for their extensive time spent providing countless resources and much-appreciated wisdom in all things bioinformatics."),
               br(),
               br(), 
               br()
      ), #end of tabPanel
      
      
      # ---------- UI QA tab ---------- 
      tabPanel("Quality Assurance",
               
               h3(strong("Quality Assurance"), style = "margin-top: 0;"),
               p("This page contains plots for quality assurance (QA) of Seurat objects that have been outputted by a CITE-seq data analysis pipeline. Dotted lines in each QA plot represent the values at which 50%, 75%, and 95% of the data falls at or below that value. Please see the \"How to Use\" guide on the Getting Started page for additional help on how to use this page."),
               
               #container for sidebar panel and main panel
               sidebarLayout(
                 
                 sidebarPanel(
                   width = 3,
                   selectInput("QA", 
                               label = "Select QA data to show",
                               choices = c("RNA Count Per Cell",
                                           "Gene Count Per Cell", 
                                           "Percent Mitochondria",
                                           "ADT Count Per Cell",
                                           "Unique ADTs Per Cell")
                   ),
                   selectInput("color_qa",
                               label = "Color QA data by:",
                               choices = c("")
                   )
                 ), #end of sidebarPanel
                 
                 mainPanel(
                   style = "padding-top: 0.6rem;",
                   width = 9,
                   fluidRow(
                     column(6, plotly::plotlyOutput("distrib_plot")),
                     column(6, plotly::plotlyOutput("box_plot"))
                   ),
                   br(), #html line break
                   br(),
                   br()
                 ) #end of mainPanel
               ) #end of sidebarLayout
      ), #end of tabPanel
      
      

      # ---------- UI Clustering tab ---------- 
      tabPanel("Clustering", 
               
               h3(strong("Clustering"), style = "margin-top: 0;"),
               p("This page contains cluster plots generated from single-cell data in the input Seurat object. Please see the \"How to Use\" guide on the Getting Started page for additional help on how to use this page."),
               
               #container for sidebar panel and main panel
               sidebarLayout(
                 
                 sidebarPanel(
                   width = 3,
                   
                   #dont hardcode options, but read in reductions slot from seurat object to detrmine which reductions are present 
                   #if no reductions present, tell user to check their file
                   #if reductions present, populate select options with what reductions the user wants to select for viewing
                   
                   selectInput("reduction", 
                               label = "Select reduction data to plot",
                               choices = c(""),
                               selected = ""
                   ),
                   selectInput("color1", 
                               label = "Color cells by:",
                               choices = c(""),
                               selected = ""
                   )
                   
                   #find a way to grab label of a cluster and 
                   #relate the x,y coords so that you can plot the label for that cluster in the center of the cluster?
                 ), #end of sidebarPanel
                 
                 mainPanel(
                   width = 9,
                   style = "padding-top: 0.6rem;",
                   
                   # Show plots here
                   fluidRow(
                     column(6,
                            plotly::plotlyOutput(outputId = "output_2dplot_1")),
                     column(6,
                            plotly::plotlyOutput(outputId = "output_3dplot_1"))
                   ),
                   
                   br(), #html line break
                   
                   #show plot selection coordinates here
                   fluidRow(
                     column(12,
                            DT::DTOutput("cluster_pg_selected")) #print metadata of selected cells 
                   ),
                   
                   br(), #html line break
                   br(),
                   br()
                 ) #end of mainPanel
               ) #end of sidebarLayout
      ), #end of tabPanel
      
      
      
      # ---------- UI Data exploration navbar menu ---------- 
      navbarMenu("Feature Expression",
                 
                 # ---------- UI Single expression tab ----------
                 tabPanel("Expression",
                          
                          h3(strong("Expression"), style = "margin-top: 0;"),
                          p("This page contains cluster plots colored by expression values generated from single-cell data in the input Seurat object. Please see the \"How to Use\" guide on the Getting Started page for additional help on how to use this page."),
                          
                          #container for sidebar panel and main panel
                          sidebarLayout(
                            sidebarPanel(
                              width = 3,
                              selectInput("reduction_expr_1d", 
                                          label = "Select reduction data to plot",
                                          choices = c(""),
                                          selected = ""),
                              uiOutput("Assay_1d"),
                              uiOutput("feature_1d")
                            ), #end of sidebarPanel
                            
                            mainPanel(
                              width = 9,
                              style = "padding-top: 0.6rem; padding-bottom: 12rem;",
                              fluidRow(
                                column(8,
                                       plotly::plotlyOutput(outputId = "exploration_reduct_1d"))
                              ),
                              
                              br(), #html line break
                              
                              #show plot selection coordinates here
                              fluidRow(
                                column(12,
                                       DT::DTOutput("expression_pg_selected")) #print metadata of selected cells 
                              ),
                              br(), #html line break
                              br(),
                              br()
                              
                            ) #end of mainPanel
                          ) #end of sidebarLayout
                          
                          ), #end of tabPanel
      
      
                # ---------- UI Co-expression tab ----------
                tabPanel("Co-Expression",
                         
                         h3(strong("Co-Expression"), style = "margin-top: 0;"),
                         p("This page contains cluster plots colored by co-expression values generated from single-cell data in the input Seurat object. Please see the \"How to Use\" guide on the Getting Started page for additional help on how to use this page."),
                         
                         #container for sidebar panel and main panel
                         sidebarLayout(
                           sidebarPanel(
                             width = 3,
                             selectInput("reduction_expr_2d", 
                                         label = "Select reduction data to plot",
                                         choices = c(""),
                                         selected = ""),
                             uiOutput("Assay_x_axis"),
                             uiOutput("x_axis_feature"),
                             uiOutput("Assay_y_axis"),
                             uiOutput("y_axis_feature")
                           ), #end of sidebarPanel
                           
                           mainPanel(
                             width = 9,
                             style = "padding-top: 0.6rem; padding-bottom: 12rem;",
                             fluidRow(
                               column(8,
                                      plotly::plotlyOutput(outputId = "exploration_reduct_2d")),
                               column(4,
                                      plotOutput(outputId = "color_legend_2d",
                                                 width = "225px", height = "225px"))
                             ),
                             
                             br(), #html line break
                             
                             #show plot selection coordinates here
                             fluidRow(
                               column(12,
                                      DT::DTOutput("coexpression_pg_selected")) #print metadata of selected cells 
                             ),
          
                             br(), #html line break
                             br(),
                             br()
                             
                           ) #end of mainPanel
                         ) #end of sidebarLayout
                ) #end of tabPanel
      ), # end of navbarMenu
      
      
      
      # ---------- UI Gating navbar menu ---------- 
      navbarMenu("Gating",
                 
                 # ---------- UI Forward Gating tab ---------- 
                 tabPanel("Forward Gating by Assay Features",
                          
                          h3(strong("Forward Gating by Assay Features"), style = "margin-top: 0;"),
                          p("This page is for gating by selected assay features such as gene transcripts (RNA) and antibody-derived tags (ADTs).
         Gating on this page is done via the lasso or selection tool in the feature scatter plot. To reset the colors of the cells in the reduction plot, double-click anywhere in the feature scatter plot.
         Please see the \"How to Use\" guide on the Getting Started page for additional help on how to use this page."),
                          #container for sidebar panel and main panel
                          sidebarLayout(
                            
                            sidebarPanel(
                              width = 3,
                              
                              uiOutput("Assay"),
                              uiOutput("x_feature"),
                              uiOutput("y_feature"),
                              selectInput("reduction_g", 
                                          label = "Select reduction data to plot",
                                          choices = c(""),
                                          selected = ""),
                              selectInput("color2", 
                                          label = "Color cells in reduction plot by:",
                                          choices = c(""),
                                          selected = ""),
                              textInput(inputId = "cell_subset_name",
                                        label = "Enter name for cell population to be gated:",
                                        value = "",
                                        placeholder = "e.g. HSCs, Blasts, CD34+CD38-"),
                              div(
                                actionButton(inputId = "reset_adt_scatter", label = "Reset Scatter"),
                                actionButton(inputId = "gate", label = "Gate")
                              )
                            ), #end of sidebarPanel
                            
                            mainPanel(
                              width = 9,
                              style = "padding-top: 0.6rem;",
                              fluidRow(
                                column(7,
                                       plotly::plotlyOutput(outputId = "featurescatter_2d")),
                                column(5,
                                       plotly::plotlyOutput(outputId = "gating_reduc_2d"))
                              ),
                              
                              br(), #html line break
                              
                              #show plot selection data here
                              fluidRow(
                                column(12,
                                       actionButton(inputId = "clear_all_gates", label = "Clear All Data"),
                                       DT::DTOutput("gating_pg_table")) #print gating data of selected cells 
                              ),
                              br(), #html line break
                              fluidRow(
                                column(12,
                                       p(strong("To view detailed gating data, download as:")),
                                       div(
                                         downloadButton("download_as_list_rds", "List (.rds)"),
                                         downloadButton("download_as_df_rds", "Dataframe (.rds)")
                                       )
                                )
                              ),
                              br(),
                              br(),
                              br()
                              
                            ) #end of main panel
                          ) #end of sidebarLayout
                 ), #end of tabPanel
                 
                 
                 # ---------- UI BackGating tab ---------- 
                 tabPanel("Backgating by Reduction",
                          
                          h3(strong("Backgating by Reduction"), style = "margin-top: 0;"),
                          p("This page is for backgating by selected assay features such as gene transcripts (RNA) and antibody-derived tags (ADTs).
                 Backgating on this page is done via the lasso or selection tool in the reduction plot. To reset the colors of the cells in the feature scatter plot, double-click anywhere in the reduction  plot.
                 Please see the \"How to Use\" guide on the Getting Started page for additional help on how to use this page."),
                          
                          #container for sidebar panel and main panel
                          sidebarLayout(
                            
                            sidebarPanel(
                              width = 3,
                              
                              selectInput("reduction_bg", 
                                          label = "Select reduction data to plot",
                                          choices = c(""),
                                          selected = ""),
                              selectInput("color2_bg", 
                                          label = "Color cells in reduction plot by:",
                                          choices = c(""),
                                          selected = ""),
                              
                              uiOutput("Assay_bg"),
                              uiOutput("x_feature_bg"),
                              uiOutput("y_feature_bg"),
                              
                              textInput(inputId = "cell_subset_name_bg",
                                        label = "Enter name for cell population to be backgated:",
                                        value = "",
                                        placeholder = "e.g. HSCs, Blasts, CD34+CD38-"),
                              div(
                                actionButton(inputId = "reset_adt_scatter_bg", label = "Reset Scatter"),
                                actionButton(inputId = "gate_bg", label = "Backgate")
                              )
                            ), #end of sidebarPanel
                            
                            mainPanel(
                              width = 9,
                              style = "padding-top: 0.6rem;",
                              fluidRow(
                                column(6,
                                       plotly::plotlyOutput(outputId = "gating_reduc_2d_bg")),
                                column(6,
                                       plotly::plotlyOutput(outputId = "featurescatter_2d_bg"))
                              ),
                              
                              br(), #html line break
                              
                              #show plot selection data here
                              fluidRow(
                                column(12,
                                       actionButton(inputId = "clear_all_gates_bg", label = "Clear All Data"),
                                       DT::DTOutput("gating_pg_table_bg")) #print gating data of selected cells 
                              ),
                              br(), #html line break
                              fluidRow(
                                column(12,
                                       p(strong("To view detailed backgating data, download as:")),
                                       div(
                                         downloadButton("download_as_list_rds_bg", "List (.rds)"),
                                         downloadButton("download_as_df_rds_bg", "Dataframe (.rds)")
                                       )
                                )
                              ),
                              br(),
                              br(),
                              br()
                              
                            ) #end of main panel
                          ) #end of sidebarLayout
                 ) # end of tabPanel
      ) # end of navbarMenu
    ) # end of navbarPage
  )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' 
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'CITEViz'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}
