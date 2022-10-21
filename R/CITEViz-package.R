#' CITEViz: CITEViz: Classifying Cell Clusters in CITE-Seq Data using the Flow Cytometry Gating Workflow
#'
#' \code{CITEViz} CITEViz is an RShiny app that allows users to
#' classify cell types in CITE-Seq data via cellular surface protein levels.
#' Similar to flow cytometry, users can physically draw gates and filter on a
#' mixture of cells, and the selected cells are immediately highlighted in
#' dimensional reduction space (e.g. UMAP, tSNE, PCA). This workflow provides
#' rapid feedback to the users to help classify cell types. To read more about
#' CITEViz, click (here. \url{https://github.com/maxsonBraunLab/CITEViz})
#'
#' The input of CITEViz is a pre-processed Seurat/SingleCelExperiment object
#' with assays containing normalized RNA/SCT (single-cell transform)
#' and ADT (antibody-derived tag) count data, as well as data for at least
#' one dimensional reduction (e.g. UMAP, PCA, etc.). The CITEViz
#' package allows users to quickly download an example dataset for the CITEViz
#' tutorial.
#'
#' Example data can be found in the (CITEVizTestData ExperimentHub \url{https://github.com/maxsonBraunLab/CITEVizTestData})
#' link. The dataset is a CITE-Seq assay of human peripheral bone and
#' mononuclear cells with 2500 cells and 228 antibodies from
#' (Hao et. al. 2021, \url{https://linkinghub.elsevier.com/retrieve/pii/S0092-8674(21)00583-3})
#' made available through the (Creative Commons BY 4.0, \url{https://creativecommons.org/licenses/by/4.0/})
#' license.
#' @author Garth Kong \email{kongga2017@gmail.com}
#' @author Thai Nguyen \email{tnguye14@uoregon.edu}
#' @author Wesley Rosales \email{wesleykrosales@gmail.com}
#' @author Anjali Panikar \email{anjali.panikar@gmail.com }
#' @author John Cheney \email{John.h.cheney@gmail.com}
#' @name CITEViz
#' @docType package
NULL
