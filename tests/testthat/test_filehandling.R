#---------------testing get_choices_from_sce()-----------------#

# create a model SingleCellExperiment object
gene_list <- c(
    "gene1", "gene2", "gene3", "gene4", "gene5",
    "gene6", "gene7", "gene8", "gene9", "gene10"
)
counts <- matrix(rpois(100, lambda = 10), ncol = 10, nrow = 10, dimnames = list(gene_list, LETTERS[1:10]))
sce <- SingleCellExperiment::SingleCellExperiment(list(logcounts = counts), mainExpName = "TEST")

sum_exp1 <- SummarizedExperiment::SummarizedExperiment(matrix(rpois(100, lambda = 5), ncol = 10, nrow = 10, dimnames = list(gene_list, LETTERS[1:10])))
# SingleCellExperiment::rowData(sum_exp1)$foo <- sample(LETTERS, nrow(sum_exp1), replace = FALSE)

test1_reduc <- matrix(rpois(20, lambda = 2), ncol = 2, nrow = 10)
SingleCellExperiment::reducedDim(sce, "pca") <- test1_reduc

SingleCellExperiment::altExp(sce, "test1") <- sum_exp1
SingleCellExperiment::colLabels(sce) <- sample(letters, ncol(sce), replace = TRUE)

# tests
test_that(desc = "get_choices_from_sce retrieves the correct choices", {
    test_choices_assays_sce <- get_choices_from_sce(category = "assays", sce_object = sce)
    test_choices_metadata_sce <- get_choices_from_sce(category = "metadata", sce_object = sce)
    test_choices_reduction_sce <- get_choices_from_sce(category = "reductions", sce_object = sce)

    expect_equal(test_choices_metadata_sce, "label")
    expect_equal(test_choices_reduction_sce, "pca")
    expect_equal(test_choices_assays_sce, c("TEST", "test1"))
})


#---------------testing get_data_from_sce()-----------------#

# tests
test_that(desc = "get_data_from_sce retrieves the correct data", {
    test_data_assays_sce <- get_data_from_sce(
        category = "assays",
        sce_object = sce, assay_name = "TEST",
        assay_data_to_get = c("gene4", "gene6")
    )
    test_data_metadata_sce <- get_data_from_sce(category = "metadata", sce_object = sce)
    test_data_reduction_sce <- get_data_from_sce(
        category = "reductions",
        sce_object = sce, reduction_name = "pca"
    )


    expect_equal(names(test_data_metadata_sce), "label")
    expect_equal(dim(test_data_metadata_sce), c(10, 1))

    expect_equal(names(test_data_reduction_sce), c("V1", "V2"))
    expect_equal(dim(test_data_reduction_sce), c(10, 2))

    expect_equal(rownames(test_data_assays_sce), LETTERS[1:10])
    expect_equal(colnames(test_data_assays_sce), c("gene4", "gene6"))
})

#---------------testing get_choices_from_seurat()-----------------#
# create fake data
gene_list <- c("gene1", "gene2", "gene3")
adt_list <- c("adt1", "adt2", "adt3")
counts <- matrix(rpois(30, lambda = 20),
    nrow = 3, ncol = 10,
    dimnames = list(gene_list, LETTERS[1:10])
)
test_so <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA")
test_adt <- Seurat::CreateAssayObject(counts = matrix(rpois(30, lambda = 20),
    nrow = 3, ncol = 10,
    dimnames = list(adt_list, LETTERS[1:10])
))
test_so[["ADT"]] <- test_adt
test_so[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rpois(20, lambda = 5),
        ncol = 2, nrow = 10,
        dimnames = list(LETTERS[1:10], c("pca_1", "pca_2"))
    ),
    key = "pca_",
    assay = "RNA"
)
test_so[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(rpois(20, lambda = 11),
        ncol = 2, nrow = 10,
        dimnames = list(LETTERS[1:10], c("umap_1", "umap_2"))
    ),
    key = "umap_",
    assay = "RNA"
)


# tests
test_that(desc = "get_choices_from_seurat retrieves the correct choices", {
    test_choices_assays_so <- get_choices_from_seurat(category = "assays", seurat_object = test_so)
    test_choices_metadata_so <- get_choices_from_seurat(category = "metadata", seurat_object = test_so)
    test_choices_reduction_so <- get_choices_from_seurat(category = "reductions", seurat_object = test_so)

    expect_equal(test_choices_metadata_so, "orig.ident")
    expect_equal(test_choices_reduction_so, c("pca", "umap"))
    expect_equal(test_choices_assays_so, c("ADT", "RNA"))
})

#---------------testing get_data_from_seurat()-----------------#
# tests
test_that(desc = "get_data_from_seurat retrieves the correct data", {
    test_data_assays_so <- get_data_from_seurat(
        category = "assays",
        seurat_object = test_so, assay_name = "RNA",
        assay_data_to_get = c("gene1", "gene2")
    )
    test_data_metadata_so <- get_data_from_seurat(category = "metadata", seurat_object = test_so)
    test_data_reduction_so <- get_data_from_seurat(
        category = "reductions",
        seurat_object = test_so,
        reduction_name = "pca"
    )

    expect_equal(names(test_data_metadata_so), c(
        "orig.ident", "nCount_RNA",
        "nFeature_RNA", "nCount_ADT",
        "nFeature_ADT"
    ))
    expect_equal(dim(test_data_metadata_so), c(10, 5))
    expect_equal(names(test_data_reduction_so), c("pca_1", "pca_2"))
    expect_equal(dim(test_data_reduction_so), c(10, 2))
    expect_equal(rownames(test_data_assays_so), LETTERS[1:10])
    expect_equal(colnames(test_data_assays_so), c("gene1", "gene2"))
})
