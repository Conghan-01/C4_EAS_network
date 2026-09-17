# ============================================================
# C4A co-expression network preservation analysis
# Prenatal -> Postnatal
#
# Input:
#   C4A_Negative_Network_pre03.csv
#   C4A_Negative_Network_post03.csv
#   prenatal_expression.residual.txt
#   postnatal_expression.residual.txt
#
# Output:
#   Network preservation statistics
#   Permutation-based Z-summary score
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

# ============================================================
# 1. Input files
# ============================================================

network_pre_file  <- "C4A_Negative_Network_pre03.csv"
network_post_file <- "C4A_Negative_Network_post03.csv"

expr_pre_file  <- "prenatal_expression.residual.txt"
expr_post_file <- "postnatal_expression.residual.txt"


# ============================================================
# 2. Parameters
# ============================================================

N_PERM <- 1000

# minimum absolute correlation considered an edge
EDGE_THRESHOLD <- 0.3

# random gene sets should have the same size as the reference network
SEED <- 12345

set.seed(SEED)


# ============================================================
# 3. Read network results
# ============================================================

net_pre <- fread(network_pre_file)
net_post <- fread(network_post_file)

cat("Prenatal network genes:", nrow(net_pre), "\n")
cat("Postnatal network genes:", nrow(net_post), "\n")


# ============================================================
# 4. Remove Ensembl version numbers
# ============================================================

clean_gene <- function(x) {
    sub("\\..*$", "", x)
}

net_pre$Gene_clean  <- clean_gene(net_pre$Gene)
net_post$Gene_clean <- clean_gene(net_post$Gene)


# ============================================================
# 5. Read residual expression matrices
# ============================================================

expr_pre <- fread(
    expr_pre_file,
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE
)

expr_post <- fread(
    expr_post_file,
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE
)

# First column = gene
gene_pre  <- clean_gene(expr_pre[[1]])
gene_post <- clean_gene(expr_post[[1]])

rownames(expr_pre)  <- gene_pre
rownames(expr_post) <- gene_post

expr_pre  <- expr_pre[, -1, drop = FALSE]
expr_post <- expr_post[, -1, drop = FALSE]

expr_pre  <- as.matrix(expr_pre)
expr_post <- as.matrix(expr_post)

storage.mode(expr_pre)  <- "numeric"
storage.mode(expr_post) <- "numeric"


cat("Prenatal expression genes:", nrow(expr_pre), "\n")
cat("Postnatal expression genes:", nrow(expr_post), "\n")


# ============================================================
# 6. Define the reference network
# ============================================================

# Here we use prenatal C4A-negative network as reference

reference_genes <- unique(net_pre$Gene_clean)

# genes available in BOTH expression datasets
common_reference_genes <- intersect(
    reference_genes,
    intersect(
        rownames(expr_pre),
        rownames(expr_post)
    )
)

cat(
    "Reference network genes available in both cohorts:",
    length(common_reference_genes),
    "\n"
)


# ============================================================
# 7. Calculate network statistics
# ============================================================

calculate_network_stats <- function(expr, genes,
                                    threshold = 0.3) {

    genes <- intersect(genes, rownames(expr))
    if (length(genes) < 5) {
        return(
            data.frame(
                density = NA,
                mean_abs_cor = NA,
                connectivity = NA,
                n_genes = length(genes)
            )
        )
    }
    x <- t(expr[genes, , drop = FALSE])
    cor_mat <- cor(
        x,
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    diag(cor_mat) <- NA

    # upper triangle
    upper <- upper.tri(cor_mat)

    cor_values <- cor_mat[upper]

    # remove NA
    cor_values <- cor_values[!is.na(cor_values)]

    # network density:
    # proportion of gene pairs with |r| >= threshold

    density <- mean(
        abs(cor_values) >= threshold
    )

    # mean absolute correlation

    mean_abs_cor <- mean(
        abs(cor_values)
    )

    # connectivity of each gene
    connectivity <- rowMeans(
        abs(cor_mat),
        na.rm = TRUE
    )

    mean_connectivity <- mean(
        connectivity,
        na.rm = TRUE
    )

    return(
        data.frame(
            density = density,
            mean_abs_cor = mean_abs_cor,
            connectivity = mean_connectivity,
            n_genes = length(genes)
        )
    )
}


# ============================================================
# 8. Calculate observed prenatal network statistics
# ============================================================

obs_pre <- calculate_network_stats(
    expr_pre,
    common_reference_genes,
    EDGE_THRESHOLD
)

obs_post <- calculate_network_stats(
    expr_post,
    common_reference_genes,
    EDGE_THRESHOLD
)

cat("\n========================================\n")
cat("Observed network statistics\n")
cat("========================================\n")

print(obs_pre)
print(obs_post)


# ============================================================
# 9. Calculate edge preservation
# ============================================================

calculate_edge_preservation <- function(expr1,
                                        expr2,
                                        genes) {

    genes <- intersect(
        genes,
        intersect(
            rownames(expr1),
            rownames(expr2)
        )
    )

    x1 <- t(expr1[genes, , drop = FALSE])
    x2 <- t(expr2[genes, , drop = FALSE])

    cor1 <- cor(
        x1,
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    cor2 <- cor(
        x2,
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    upper <- upper.tri(cor1)

    r1 <- cor1[upper]
    r2 <- cor2[upper]

    keep <- !is.na(r1) & !is.na(r2)

    r1 <- r1[keep]
    r2 <- r2[keep]

    # correlation of edge weights
    edge_preservation <- cor(
        r1,
        r2,
        method = "pearson"
    )

    # preservation of absolute edge strength
    abs_edge_preservation <- cor(
        abs(r1),
        abs(r2),
        method = "pearson"
    )

    return(
        c(
            edge_preservation = edge_preservation,
            abs_edge_preservation = abs_edge_preservation
        )
    )
}


observed_edge <- calculate_edge_preservation(
    expr_pre,
    expr_post,
    common_reference_genes
)

cat("\n========================================\n")
cat("Observed edge preservation\n")
cat("========================================\n")

print(observed_edge)


# ============================================================
# 10. Permutation analysis
# ============================================================

# Universe:
# genes that are expressed in BOTH cohorts

gene_universe <- intersect(
    rownames(expr_pre),
    rownames(expr_post)
)

gene_universe <- setdiff(
    gene_universe,
    "C4A"
)

network_size <- length(common_reference_genes)

cat("\nRandom gene-set size:", network_size, "\n")
cat("Permutation number:", N_PERM, "\n")


perm_density_pre <- numeric(N_PERM)
perm_density_post <- numeric(N_PERM)

perm_edge <- numeric(N_PERM)
perm_abs_edge <- numeric(N_PERM)


# ------------------------------------------------------------
# permutation loop
# ------------------------------------------------------------

for (i in seq_len(N_PERM)) {

    random_genes <- sample(
        gene_universe,
        size = network_size,
        replace = FALSE
    )

    stat_pre <- calculate_network_stats(
        expr_pre,
        random_genes,
        EDGE_THRESHOLD
    )

    stat_post <- calculate_network_stats(
        expr_post,
        random_genes,
        EDGE_THRESHOLD
    )

    perm_density_pre[i] <- stat_pre$density
    perm_density_post[i] <- stat_post$density

    edge_tmp <- calculate_edge_preservation(
        expr_pre,
        expr_post,
        random_genes
    )

    perm_edge[i] <- edge_tmp["edge_preservation"]
    perm_abs_edge[i] <- edge_tmp["abs_edge_preservation"]

    if (i %% 100 == 0) {
        cat("Permutation:", i, "/", N_PERM, "\n")
    }
}


# ============================================================
# 11. Calculate preservation statistics
# ============================================================

# ------------------------------------------------------------
# Density preservation
# ------------------------------------------------------------

observed_density_ratio <-
    obs_post$density / obs_pre$density

random_density_ratio <-
    perm_density_post / perm_density_pre

random_density_ratio <-
    random_density_ratio[
        is.finite(random_density_ratio)
    ]


# ------------------------------------------------------------
# Edge preservation
# ------------------------------------------------------------

perm_edge <- perm_edge[
    is.finite(perm_edge)
]

perm_abs_edge <- perm_abs_edge[
    is.finite(perm_abs_edge)
]


# ============================================================
# 12. Z-score transformation
# ============================================================

density_Z <- (
    observed_density_ratio -
    mean(random_density_ratio)
) / sd(random_density_ratio)


edge_Z <- (
    observed_edge["edge_preservation"] -
    mean(perm_edge)
) / sd(perm_edge)


abs_edge_Z <- (
    observed_edge["abs_edge_preservation"] -
    mean(perm_abs_edge)
) / sd(perm_abs_edge)


# ============================================================
# 13. Composite Z-summary-like preservation score
# ============================================================

Zsummary <- mean(
    c(
        density_Z,
        edge_Z,
        abs_edge_Z
    ),
    na.rm = TRUE
)


# ============================================================
# 14. Empirical P values
# ============================================================

density_p <- (
    1 +
    sum(
        random_density_ratio >= observed_density_ratio
    )
) / (N_PERM + 1)


edge_p <- (
    1 +
    sum(
        perm_edge >= observed_edge["edge_preservation"]
    )
) / (N_PERM + 1)


abs_edge_p <- (
    1 +
    sum(
        perm_abs_edge >= observed_edge["abs_edge_preservation"]
    )
) / (N_PERM + 1)


# ============================================================
# 15. Final results
# ============================================================

result <- data.frame(

    Reference = "Prenatal_C4A_positive_gene",

    N_genes = length(common_reference_genes),

    Prenatal_density =
        obs_pre$density,

    Postnatal_density =
        obs_post$density,

    Density_ratio =
        observed_density_ratio,

    Density_Z =
        density_Z,

    Edge_preservation =
        observed_edge["edge_preservation"],

    Edge_Z =
        edge_Z,

    Absolute_edge_preservation =
        observed_edge["abs_edge_preservation"],

    Absolute_edge_Z =
        abs_edge_Z,

    Zsummary =
        Zsummary,

    Density_P =
        density_p,

    Edge_P =
        edge_p,

    Absolute_edge_P =
        abs_edge_p
)


print(result)

# ============================================================
# 16. Save results
# ============================================================

write.csv(
    result,
    "C4A_Negative_pre_as-reference_Network_Preservation_Zsummary.csv",
    row.names = FALSE
)

# ============================================================
# 17. Save permutation distributions
# ============================================================

perm_result <- data.frame(
    Density_Ratio = random_density_ratio,
    Edge_Preservation = perm_edge[
        seq_len(min(length(perm_edge),
                    length(random_density_ratio)))
    ]
)

write.csv(
    perm_result,
    "C4A_Negative_pre_as-reference_Network_Preservation_Permutations.csv",
    row.names = FALSE
)


cat("\n========================================\n")
cat("Network preservation analysis completed\n")
cat("========================================\n")
cat(
    "Number of reference genes:",
    length(common_reference_genes),
    "\n"
)

cat(
    "Zsummary-like preservation score:",
    round(Zsummary, 3),
    "\n"
)

cat(
    "Output:",
    "C4A_Negative_Network_Preservation_Zsummary.csv",
    "\n"
)


######################### Reverse analysis: postnatal as reference ################
# ============================================================
# C4A co-expression network preservation analysis
# Postnatal -> Prenatal (Reverse Direction)
#
# Input:
#   C4A_Negative_Network_post03.csv (or Positive)
#   prenatal_expression.residual.txt
#   postnatal_expression.residual.txt
#
# Output:
#   Network preservation statistics (Postnatal -> Prenatal)
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})


# ============================================================
# 1. Input files (Set Postnatal as Reference)
# ============================================================

# 注意：此处将参考网络换成产后文件
network_ref_file  <- "C4A_Positive_Network_post03.csv"
expr_pre_file  <- "prenatal_expression.residual.txt"
expr_post_file <- "postnatal_expression.residual.txt"


# ============================================================
# 2. Parameters
# ============================================================

N_PERM <- 1000

# minimum absolute correlation considered an edge
EDGE_THRESHOLD <- 0.3

SEED <- 12345
set.seed(SEED)


# ============================================================
# 3. Read network results
# ============================================================

net_ref <- fread(network_ref_file)
cat("Reference network genes (Postnatal):", nrow(net_ref), "\n")


# ============================================================
# 4. Remove Ensembl version numbers
# ============================================================

clean_gene <- function(x) {
    sub("\\..*$", "", x)
}

net_ref$Gene_clean <- clean_gene(net_ref$Gene)


# ============================================================
# 5. Read residual expression matrices
# ============================================================

expr_pre <- fread(
    expr_pre_file,
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE
)

expr_post <- fread(
    expr_post_file,
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE
)

# First column = gene
gene_pre  <- clean_gene(expr_pre[[1]])
gene_post <- clean_gene(expr_post[[1]])

rownames(expr_pre)  <- gene_pre
rownames(expr_post) <- gene_post

expr_pre  <- expr_pre[, -1, drop = FALSE]
expr_post <- expr_post[, -1, drop = FALSE]

expr_pre  <- as.matrix(expr_pre)
expr_post <- as.matrix(expr_post)

storage.mode(expr_pre)  <- "numeric"
storage.mode(expr_post) <- "numeric"


# ============================================================
# 6. Define the reference network
# ============================================================

# 使用产后 C4A 负相关网络作为参考
reference_genes <- unique(net_ref$Gene_clean)

# genes available in BOTH expression datasets
common_reference_genes <- intersect(
    reference_genes,
    intersect(
        rownames(expr_pre),
        rownames(expr_post)
    )
)

cat(
    "Reference network genes available in both cohorts:",
    length(common_reference_genes),
    "\n"
)


# ============================================================
# 7. Calculate network statistics function
# ============================================================

calculate_network_stats <- function(expr, genes, threshold = 0.3) {
    genes <- intersect(genes, rownames(expr))

    if (length(genes) < 5) {
        return(
            data.frame(
                density = NA,
                mean_abs_cor = NA,
                connectivity = NA,
                n_genes = length(genes)
            )
        )
    }

    x <- t(expr[genes, , drop = FALSE])

    cor_mat <- cor(
        x,
        method = "pearson",
        use = "pairwise.complete.obs"
    )

    diag(cor_mat) <- NA

    upper <- upper.tri(cor_mat)
    cor_values <- cor_mat[upper]
    cor_values <- cor_values[!is.na(cor_values)]

    density <- mean(abs(cor_values) >= threshold)
    mean_abs_cor <- mean(abs(cor_values))
    connectivity <- rowMeans(abs(cor_mat), na.rm = TRUE)
    mean_connectivity <- mean(connectivity, na.rm = TRUE)

    return(
        data.frame(
            density = density,
            mean_abs_cor = mean_abs_cor,
            connectivity = mean_connectivity,
            n_genes = length(genes)
        )
    )
}


# ============================================================
# 8. Calculate observed network statistics
# 注意：此时 Postnatal 是参考（Baseline），Prenatal 是待测目标
# ============================================================

obs_post_ref <- calculate_network_stats(
    expr_post,
    common_reference_genes,
    EDGE_THRESHOLD
)

obs_pre_target <- calculate_network_stats(
    expr_pre,
    common_reference_genes,
    EDGE_THRESHOLD
)

cat("\n========================================\n")
cat("Observed network statistics\n")
cat("========================================\n")
cat("Postnatal (Reference):\n")
print(obs_post_ref)
cat("Prenatal (Target):\n")
print(obs_pre_target)


# ============================================================
# 9. Calculate edge preservation function
# ============================================================

calculate_edge_preservation <- function(expr1, expr2, genes) {
    genes <- intersect(
        genes,
        intersect(
            rownames(expr1),
            rownames(expr2)
        )
    )

    x1 <- t(expr1[genes, , drop = FALSE])
    x2 <- t(expr2[genes, , drop = FALSE])

    cor1 <- cor(x1, method = "pearson", use = "pairwise.complete.obs")
    cor2 <- cor(x2, method = "pearson", use = "pairwise.complete.obs")

    upper <- upper.tri(cor1)
    r1 <- cor1[upper]
    r2 <- cor2[upper]

    keep <- !is.na(r1) & !is.na(r2)
    r1 <- r1[keep]
    r2 <- r2[keep]

    edge_preservation <- cor(r1, r2, method = "pearson")
    abs_edge_preservation <- cor(abs(r1), abs(r2), method = "pearson")

    return(
        c(
            edge_preservation = edge_preservation,
            abs_edge_preservation = abs_edge_preservation
        )
    )
}

# 传入顺序：expr1 = Postnatal (参考), expr2 = Prenatal (目标)
observed_edge <- calculate_edge_preservation(
    expr_post,
    expr_pre,
    common_reference_genes
)

cat("\n========================================\n")
cat("Observed edge preservation (Post -> Pre)\n")
cat("========================================\n")
print(observed_edge)


# ============================================================
# 10. Permutation analysis
# ============================================================

gene_universe <- intersect(
    rownames(expr_pre),
    rownames(expr_post)
)
gene_universe <- setdiff(gene_universe, "C4A")

network_size <- length(common_reference_genes)

cat("\nRandom gene-set size:", network_size, "\n")
cat("Permutation number:", N_PERM, "\n")

perm_density_post <- numeric(N_PERM)
perm_density_pre  <- numeric(N_PERM)
perm_edge         <- numeric(N_PERM)
perm_abs_edge     <- numeric(N_PERM)

for (i in seq_len(N_PERM)) {
    random_genes <- sample(
        gene_universe,
        size = network_size,
        replace = FALSE
    )

    stat_post <- calculate_network_stats(expr_post, random_genes, EDGE_THRESHOLD)
    stat_pre  <- calculate_network_stats(expr_pre, random_genes, EDGE_THRESHOLD)

    perm_density_post[i] <- stat_post$density
    perm_density_pre[i]  <- stat_pre$density

    edge_tmp <- calculate_edge_preservation(expr_post, expr_pre, random_genes)

    perm_edge[i]     <- edge_tmp["edge_preservation"]
    perm_abs_edge[i] <- edge_tmp["abs_edge_preservation"]

    if (i %% 100 == 0) {
        cat("Permutation:", i, "/", N_PERM, "\n")
    }
}


# ============================================================
# 11. Calculate preservation statistics
# ============================================================

observed_density_ratio <- obs_pre_target$density / obs_post_ref$density
random_density_ratio   <- perm_density_pre / perm_density_post
random_density_ratio   <- random_density_ratio[is.finite(random_density_ratio)]

perm_edge     <- perm_edge[is.finite(perm_edge)]
perm_abs_edge <- perm_abs_edge[is.finite(perm_abs_edge)]


# ============================================================
# 12. Z-score transformation
# ============================================================

density_Z <- (observed_density_ratio - mean(random_density_ratio)) / sd(random_density_ratio)
edge_Z    <- (observed_edge["edge_preservation"] - mean(perm_edge)) / sd(perm_edge)
abs_edge_Z <- (observed_edge["abs_edge_preservation"] - mean(perm_abs_edge)) / sd(perm_abs_edge)


# ============================================================
# 13. Composite Z-summary-like preservation score
# ============================================================

Zsummary <- mean(c(density_Z, edge_Z, abs_edge_Z), na.rm = TRUE)


# ============================================================
# 14. Empirical P values
# ============================================================

density_p   <- (1 + sum(random_density_ratio >= observed_density_ratio)) / (N_PERM + 1)
edge_p      <- (1 + sum(perm_edge >= observed_edge["edge_preservation"])) / (N_PERM + 1)
abs_edge_p  <- (1 + sum(perm_abs_edge >= observed_edge["abs_edge_preservation"])) / (N_PERM + 1)


# ============================================================
# 15. Final results
# ============================================================

result <- data.frame(
    Reference = "Postnatal_C4A_positive",
    N_genes = length(common_reference_genes),
    Postnatal_reference_density = obs_post_ref$density,
    Prenatal_target_density     = obs_pre_target$density,
    Density_ratio               = observed_density_ratio,
    Density_Z                   = density_Z,
    Edge_preservation           = observed_edge["edge_preservation"],
    Edge_Z                      = edge_Z,
    Absolute_edge_preservation  = observed_edge["abs_edge_preservation"],
    Absolute_edge_Z             = abs_edge_Z,
    Zsummary                    = Zsummary,
    Density_P                   = density_p,
    Edge_P                      = edge_p,
    Absolute_edge_P             = abs_edge_p
)

print(result)


# ============================================================
# 16. Save results
# ============================================================

write.csv(
    result,
    "C4A_Positive_post_as-reference_Network_Preservation_Zsummary.csv",
    row.names = FALSE
)

cat("\n========================================\n")
cat("Reverse network preservation analysis completed\n")
cat("========================================\n")
cat("Zsummary-like preservation score (Post -> Pre):", round(Zsummary, 3), "\n")


perm_result <- data.frame(
    Density_Ratio = random_density_ratio,
    Edge_Preservation = perm_edge[
        seq_len(min(length(perm_edge),
                    length(random_density_ratio)))
    ]
)

write.csv(
    perm_result,
    "C4A_Positive_post_as-reference_Network_Preservation_Permutations.csv",
    row.names = FALSE
)