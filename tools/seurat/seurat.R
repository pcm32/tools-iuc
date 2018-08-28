options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(dplyr)
    library(optparse)
})

option_list <- list(
    make_option(c("-ctrl","--ctrl"), type="character", help="Counts file"),
    make_option(c("-rctrl","--rctrl"), type="character", help="Counts file as RDS object"),
    make_option(c("-stim","--stim"), type="character", help="A second counts file (optional)"),
    make_option(c("-rstim","--rstim"), type="character", help="A second counts file as RDS object (optional)"),
    make_option(c("-numPCs","--numPCs"), type="integer", help="Number of PCs to use in plots"),
    make_option(c("-min.cells","--min.cells"), type="integer", help="Minimum cells to include"),
    make_option(c("-min.genes","--min.genes"), type="integer", help="Minimum genes to include"),
    make_option(c("-low.thresholds","--low.thresholds"), type="double", help="Low threshold for filtering cells"),
    make_option(c("-high.thresholds","--high.thresholds"), type="double", help="High threshold for filtering cells"),
    make_option(c("-x.low.cutoff","--x.low.cutoff"), type="double", help="X-axis low cutoff for variable genes"),
    make_option(c("-x.high.cutoff","--x.high.cutoff"), type="double", help="X-axis high cutoff for variable genes"),
    make_option(c("-y.cutoff","--y.cutoff"), type="double", help="Y-axis cutoff for variable genes"),
    make_option(c("-cells.use","--cells.use"), type="integer", help="Cells to use for PCHeatmap"),
    make_option(c("-resolution","--resolution"), type="double", help="Resolution in FindClusters"),
    make_option(c("-min.pct","--min.pct"), type="double", help="Minimum percent cells in FindClusters"),
    make_option(c("-logfc.threshold","--logfc.threshold"), type="double", help="LogFC threshold in FindClusters"),
    make_option(c("-cca","--cca"), type="logical", help="If have 2 samples run CCA"),
    make_option(c("-rdata","--rdata"), type="logical", help="Output Seurat RData file")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

runSeurat <- function(counts) {
    seuset <- CreateSeuratObject(raw.data = counts, min.cells = args$min.cells, min.genes = args$min.cells)
    p1 <- VlnPlot(object = seuset, features.plot = c("nGene", "nUMI"), nCol = 2)
    p2 <- GenePlot(object = seuset, gene1 = "nUMI", gene2 = "nGene")
    print(p1)
    print(p2)

    print("Filtering cells")
    if (!is.null(args$low.thresholds)){
        lowthresh <- args$low.thresholds
    } else {
        lowthresh <- "-Inf"
    }
    if (!is.null(args$high.thresholds)){
        highthresh <- args$high.thresholds
    } else {
        highthresh <- "Inf"
    }
    seuset <- FilterCells(object = seuset, subset.names = c("nUMI"), low.thresholds=c(lowthresh), high.thresholds = c(highthresh))

    print("Normalizing the data")
    seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", scale.factor = 10000)

    print("Finding variable genes")
    seuset <- FindVariableGenes(object = seuset, 
        mean.function = ExpMean,
        dispersion.function = LogVMR, 
        x.low.cutoff = args$x.low.cutoff, 
        x.high.cutoff = args$x.high.cutoff,
        y.cutoff = args$y.cutoff
    )

    print("Scaling the data and removing unwanted sources of variation")
    #temporarily removed , vars.to.regress = c("nUMI")
    seuset <- ScaleData(object = seuset)

    print("Performing PCA analysis")
    seuset <- RunPCA(object = seuset, pc.genes = seuset@var.genes)
    VizPCA(object = seuset, pcs.use = 1:2)
    PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2)
    PCHeatmap(
        object = seuset, 
        pc.use = 1:args$numPCs, 
        cells.use = args$cell.use, 
        do.balanced = TRUE, 
        label.columns = FALSE,
        use.full = FALSE
    )

    print("Determining statistically significant principal components")
    seuset <- JackStraw(object = seuset, num.replicate = 100, display.progress= FALSE)
    JackStrawPlot(object = seuset, PCs = 1:args$numPCs)
    p <- PCElbowPlot(object = seuset)
    print(p)

    print("Clustering the cells")
    seuset <- FindClusters(
        object = seuset, 
        reduction.type = "pca", 
        dims.use = 1:args$numPCs, 
        resolution = args$resolution,
        print.output = 0, 
        save.SNN = TRUE
    )

    print("Running non-linear dimensional reduction (tSNE)")
    # added check_duplicates = FALSE to get test to pass https://github.com/satijalab/seurat/issues/167
    seuset <- RunTSNE(object = seuset, dims.use = 1:args$numPCs, do.fast = TRUE, check_duplicates = FALSE)
    TSNEPlot(object = seuset)

    print("Finding differentially expressed genes (cluster biomarkers)")
    markers <- FindAllMarkers(object = seuset, only.pos = TRUE, min.pct = args$min.pct,
        logfc.threshold = args$logfc.threshold)
    top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
    p <- DoHeatmap(object = seuset, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
    print(p)
    return(seuset)
}

# Open PDF for plots
pdf("out.pdf")

if (!is.null(args$rctrl)) {
    ctrl <- readRDS(args$rctrl)
} else {
    ctrl <- read.delim(args$ctrl, row.names=1)
    ctrl <- runSeurat(ctrl)
}

if (!is.null(args$rstim)) {
    stim <- readRDS(args$rstim)
} else if (!is.null(args$stim)) {
    stim <- read.delim(args$stim, row.names=1)
    stim <- runSeurat(stim)
}


if (!is.null(args$cca)) {

    ctrl@meta.data$stim <- "CTRL"
    stim@meta.data$stim <- "STIM"

    print("Selecting genes for input to CCA")
    g.1 <- head(rownames(ctrl@hvg.info), 1000)
    g.2 <- head(rownames(stim@hvg.info), 1000)
    genes.use <- unique(c(g.1, g.2))
    genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
    genes.use <- intersect(genes.use, rownames(stim@scale.data))

    # remove cell cycle by running a reression on cell cycle genes
    # note that G2M/S phase cell cycle appears in PC2 in both datasets
    #cc.genes <- readLines("~/Downloads/IntegratedAnalysis_Examples/data/regev_lab_cell_cycle_genes.txt", warn = F)
    #mars <- RunPCA(mars, pc.genes = cc.genes)
    #ss2 <- RunPCA(ss2, pc.genes = cc.genes)
    #now we rescale the data, adding in the CC score as a regression covariate
    #mars <- ScaleData(mars, vars.to.regress = "PC2", genes.use = hvg.union)
    #ss2 <- ScaleData(ss2, vars.to.regress = "PC2", genes.use = hvg.union)

    print("Performing canonical correlation analysis (CCA)")
    combined <- RunCCA(ctrl, stim, genes.use = genes.use)

    print("Visualising CCA results")
    p1 <- DimPlot(object = combined, reduction.use = "cca", group.by = "stim", pt.size = 0.5, do.return = TRUE)
    p2 <- VlnPlot(object = combined, features.plot = "CC1", group.by = "stim", do.return = TRUE)
    p <- plot_grid(p1, p2)
    print(p)
    MetageneBicorPlot(combined, grouping.var = "stim", dims.eval = 1:args$numPCs, display.progress = FALSE)
    DimHeatmap(object = combined, reduction.type = "cca", cells.use = args$cells.use, dim.use = 1:args$numPCs, do.balanced = TRUE)

    # Discard cells where the variance explained by CCA is <2-fold  (ratio < 0.5) compared to PCA
    print("Discarding uninformative cells")
    combined <- CalcVarExpRatio(combined, reduction.type = "pca", grouping.var = "stim", dims.use = 1:args$numPCs)
    combined.all <- combined
    combined <- SubsetData(combined, subset.name = "var.ratio.pca", accept.low = 0.5)

    print("Aligning subspaces")
    combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "stim", dims.align = 1:args$numPCs)
    p1 <- VlnPlot(object = combined, features.plot = "ACC1", group.by = "stim", do.return = TRUE)
    p2 <- VlnPlot(object = combined, features.plot = "ACC2", group.by = "stim", do.return = TRUE)
    p <- plot_grid(p1, p2)
    print(p)


    print("Performing integrated analysis")
    # added check_duplicates = FALSE to get test to pass https://github.com/satijalab/seurat/issues/167
    combined <- RunTSNE(combined, reduction.use = "cca.aligned", dims.use = 1:args$numPCs, do.fast = T, check_duplicates = FALSE)
    p1 <- TSNEPlot(combined, do.return = T, pt.size = 0.5, group.by = "stim")
    p2 <- TSNEPlot(combined, do.label = T, do.return = T, pt.size = 0.5)
    p <- plot_grid(p1, p2)
    print(p)
}

# Close PDF for plots
dev.off()

if (!is.null(args$rdata)) {
    if (!is.null(args$stim)) {
        if (!is.null(args$cca)) {
            save(ctrl, stim, combined, file="Seurat_analysis.RData")
        } else {
            save(ctrl, stim, file="Seurat_analysis.RData")
        }
    } else {
        save(ctrl, file="Seurat_analysis.RData")
    }
} 

sessionInfo()