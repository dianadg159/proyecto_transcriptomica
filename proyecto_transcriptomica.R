# Cargamos la librería recount
library("recount3")
# Obtenemos el link de la base de datos
getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)
# Cambiamos el url
options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")
# Confirmando que se cambió el URL
getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)
# Gradamos los proyectos disponibles de humano
human_projects <- available_projects(organism = "human")
# Encontramos proyecto SRP107565 perteneciente al artículo
# "Multiomics Profiling Establishes the Polypharmacology of FDA-Approved CDK4/6 Inhibitors and the Potential for Differential Clinical Activity"
project_info <- subset(
  human_projects,
  project == "SRP107565"
)

# Cremos un objeto de tipo RangedSummarizedExperiment (RSE)
# con la información a nivel de genes
rse_prueba <- create_rse(project_info)
dim(rse_prueba)
## vemos dimensiones de 63856 (genes)  x 216 (muestras)
# hacer el objeto
# Sacamos las counts
assay(rse_prueba, "counts") <- compute_read_counts(rse_prueba)
prueba_matrix <- assay(rse_prueba, "counts")
# Expandemos tabla de atributos
rse_prueba <- expand_sra_attributes(rse_prueba)
# Crear tabla de atributos de las muestras
rse_prueba_data <- colData(rse_prueba)[
  ,grepl("^sra_attribute", colnames(colData(rse_prueba)))
]
# Sacamos solamente la línea celular MCF7 y de 24 horas
good_samples2 <- rse_prueba_data[,2] == "MCF7" & rse_prueba_data[,5] == "24 hr"
rse_prueba_data_MCF7_24hr <- rse_prueba_data[good_samples2,]
# Sacamos nombre de las muestras
dose <- as.character(rse_prueba_data_MCF7_24hr[,3])
dose <- substr(dose, 1, nchar(dose)-3)
sample_names <- paste(rse_prueba_data_MCF7_24hr[,1], dose, sep = "_")

# Hacemos el objeto DEGList
library("edgeR")
dge_prueba <- DGEList(
  counts = prueba_matrix,
  genes = rowData(rse_prueba)
)
# Normalizamos datos
dge_prueba <- calcNormFactors(dge_prueba)
# Quitamos los genes que solo están muy poco expresados
cutoff <- 1
drop <- which(apply(cpm(dge_prueba), 1, max) < cutoff)
d_prueba <- dge_prueba[-drop,]
d_prueba <- d_prueba[,good_samples2]

# Visualizamos con un MDS plot
library("RColorBrewer")
# Convirtiendo los grupos de tratamiento a color para mejor visualización
agent2 <- rse_prueba_data_MCF7_24hr[,1]
col.group <- factor(agent2)
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
# Hacemos el plot
plotMDS(d_prueba, labels = sample_names, col = col.group)

# Visualizamos con un PCA
library(factoextra)
library(FactoMineR)

pca.raw.d <- log2(d_prueba$counts + 0.5)
colnames(pca.raw.d) <- sample_names
pca.d <- PCA(t(pca.raw.d), graph = F)
fviz_pca_ind(pca.d, col.ind = agent2)

# Transformar los datos
library("limma")
sample_names <- factor(sample_names)
#factor(agent2)
mm <- model.matrix(~ 0 + sample_names)
# voom: Mean-variance trend debe asumir un patrón conforme a la y
voom.y.d <- voom(d_prueba, mm, plot = TRUE)

# Modificar los datos a un modelo lineal
fit <- lmFit(voom.y.d, mm)
coef.fit <- fit$coefficients
head(coef(fit))

# Hacemos el análisis diferencial
# Hacemos la estadística de Bayes Empírica para Expresión Diferencial
eb_results <- eBayes(fit)
# sacar los genes del modelo lineal
de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_prueba),
  sort.by = "none"
)
head(de_results)

# Comparamos con coeficientes
coef_DEG <- coef.fit[rownames(coef.fit) %in% rownames(de_results)[1:5],]

# Visualizamos los genes diferencialmente expresados más significativos
# Agarramos los genes que tengan un p-valor menor a 0.05
deg_pval <- de_results$adj.P.Val < 0.05
deg_20 <- de_results[deg_pval,]
deg_20 <- deg_20[order(deg_20$adj.P.Val),]
deg_20 <- as.data.frame(deg_20)
deg_20 <- de_results[rank(de_results$adj.P.Val) <= 20,]
ggplot(deg_20 , aes(x = adj.P.Val, y = logFC)) +
  geom_text(label = deg_20$gene_name)

# Agarramos los primeros cincuenta genes de los 28740
# que tienen un p-valor menor a 0.05
heatmap <- voom.y.d$E[rank(de_results$adj.P.Val) <= 50, ]
de_results[rank(de_results$adj.P.Val) <= 50, ]
# Creemos una tabla con las muestras
df_heatmap <- as.data.frame(rse_prueba_data_MCF7_24hr[,-4])
colnames(df_heatmap) <- c("Agent", "Cell line", "Dose (uM)", "Time (hr)")
# Extraer los nombres de los genes
gnames <- rownames(heatmap)
rownames(heatmap) <- de_results$gene_name[
  match(rownames(heatmap), rownames(de_results))
]

## Heatmap
library("pheatmap")
pheatmap(
  heatmap2,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df_heatmap
)

#Sacamos los gene_ids de los genes para el background
library(strex)
file.create("gene_ids.txt")
gene_ids <- de_results$gene_id
gene_ids <- str_before_last_dot(gene_ids)
write.table(gene_ids, file = "gene_ids.txt", sep = "\n", row.names = FALSE, quote = FALSE)
# Sacamos la lista de los DEGsde mayor a menor expresados
gene_names <- noquote(deg_20[order(abs(deg_20$logFC), decreasing = TRUE),]$gene_name)
