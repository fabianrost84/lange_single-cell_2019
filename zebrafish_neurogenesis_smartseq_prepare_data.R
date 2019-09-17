library(AnnotationHub)
library(data.table)
library(magrittr)
library(Matrix)
library(scater)
library(reticulate)
ad <- import("anndata", convert = FALSE)

file.list <- list.files("/scratch/rulands/zebrafish_brain_christian_lange/180120_SmartSeq/",full.names = T,pattern="\\.txt$")
dl <- lapply(file.list, 
             function(x) data.table::fread(x, header=T,skip=1L) %>%
               .[,-c("Chr","Start","End","Strand", "Length")]
)

# check data consistency
if (!identical(dl[[1]][, 1], dl[[2]][, 1])) {
  warning("Genes not identical")
}
if (sum(colnames(dl[[1]][, -1]) %in% colnames(dl[[2]][, -1])) > 0) {
  warning("Duplicate cell names")
}

d <- Reduce(function(x,y) merge(x,y, by="Geneid"),dl)
data.table::setnames(d, "Geneid","ID")

cell.names <- colnames(d)[-1]
count.matrix <- as.matrix(d[ ,-"ID", with=FALSE]) %>% Matrix(., sparse = TRUE)
rownames(count.matrix) <- (d[, ID])

file.name <- mapply(function(d, file.name) replicate(length(d)-1, file.name), dl, file.list) %>% Reduce(append, .) %>% substr(., 67, 72)

fish.name <- regexpr("_[0-9]{2}_", cell.names) %>% 
  regmatches(cell.names, .) %>% 
  substr(., 2, 3) %>% 
  paste(file.name, ., sep="_")

is.bulk = grepl("bulk", cell.names) | (fish.name=="855_30")
is.bulk[[which(cell.names == "RG-NBN30_01_C08")]] = TRUE

cell.type = regexpr('^[A-Z]{2,3}', cell.names) %>% 
  regmatches(cell.names, .) %>% 
  paste0(., lapply(is.bulk, function (x) if (x) {"_bulk"} else {""}))
cell.type[[which(cell.names == "RG-NBN30_01_C08")]] = "RG-NBN"

is.nc = (cell.type == "NC")
is.rna = (cell.type == "RNA")

file.summary.list <- list.files("/scratch/rulands/zebrafish_brain_christian_lange/180120_SmartSeq/",full.names = T,pattern="\\.summary$")
dl.summary <- lapply(file.summary.list, 
                     function(x) data.table::fread(x, header=T, key = "Status") 
)
d.summary <- Reduce(function(x, y) merge(x, y), dl.summary)  %>% 
  melt(., id="Status", variable.name = "cell.names", value.name = "Counts")  %>%
  dcast(., cell.names ~ Status, value.var = "Counts")
d.summary[, Total_Counts := Reduce("+", .SD), .SDcol = !"cell.names"]
d.summary[, Fraction_Mapped := Assigned / Total_Counts]

if(!unique(d.summary[, cell.names] == cell.names)) {
  warning("Something is wrong with the cell names in the summary files")
}

coldata <- data.frame(sample_id = cell.names, 
                      batch = file.name, 
                      cell.type = cell.type,
                      is.bulk = is.bulk,
                      is.nc = is.nc,
                      is.rna = is.rna,
                      file.name = file.name,
                      fish.name = fish.name,
                      fraction.mapped = d.summary[, Fraction_Mapped])

# pheno.data <- new("AnnotatedDataFrame",anno)
rownames(coldata) <- cell.names

sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = count.matrix), colData = coldata)

ah <- AnnotationHub()
EnsDb.Drerio.v87 <- query(ah, c("EnsDb", "Danio", "87"))[[1]]

rowData(sce)$ID <- rownames(sce)
rowData(sce)$Symbol <- c(mapIds(EnsDb.Drerio.v87, rowData(sce)$ID, "SYMBOL", "GENEID"), use.names=FALSE)
rowData(sce)$CHR <- c(mapIds(EnsDb.Drerio.v87, keys=rownames(sce), 
                             column="SEQNAME", keytype="GENEID"), use.names=FALSE)

rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
head(rownames(sce))

sce <- sce[, !(sce$is.bulk | sce$is.nc | sce$is.rna)]
sce <- scater::calculateQCMetrics(sce,
                                  feature_controls = list(ERCC=grepl("ERCC", rownames(sce)),
                                                          Mito=which(rowData(sce)$CHR == "MT")))

anndata.object <- ad$AnnData(X = t(counts(sce)), 
                             obs = as.data.frame(colData(sce)), 
                             var = as.data.frame(rowData(sce)))
anndata.object$write("/home/fabrost/pksCloud/projects/zebrafish_neurogenesis_smartseq/results/zebrafish_neurogenesis_smartseq.h5ad")
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = count.matrix), colData = coldata)
print("Done")