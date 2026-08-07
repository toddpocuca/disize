contraster <- function(dds, # should contain colData and design
                       group1, # list of character vectors each with 2 or more items
                       group2, # list of character vectors each with 2 or more items
                       weighted = FALSE) {
    mod_mat <- model.matrix(DESeq2::design(dds), SummarizedExperiment::colData(dds))

    grp1_rows <- list()
    grp2_rows <- list()


    for (i in seq_along(group1)) {
        grp1_rows[[i]] <- SummarizedExperiment::colData(dds)[[group1[[i]][1]]] %in% group1[[i]][2:length(group1[[i]])]
    }


    for (i in seq_along(group2)) {
        grp2_rows[[i]] <- SummarizedExperiment::colData(dds)[[group2[[i]][1]]] %in% group2[[i]][2:length(group2[[i]])]
    }

    grp1_rows <- Reduce(function(x, y) x & y, grp1_rows)
    grp2_rows <- Reduce(function(x, y) x & y, grp2_rows)

    mod_mat1 <- mod_mat[grp1_rows, , drop = F]
    mod_mat2 <- mod_mat[grp2_rows, , drop = F]

    if (!weighted) {
        mod_mat1 <- mod_mat1[!duplicated(mod_mat1), , drop = F]
        mod_mat2 <- mod_mat2[!duplicated(mod_mat2), , drop = F]
    }

    return(colMeans(mod_mat1) - colMeans(mod_mat2))
}
