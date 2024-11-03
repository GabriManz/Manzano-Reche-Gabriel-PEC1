if (!(require(limma))){
  source("http://bioconductor.org/biocLite.R")
  biocLite("limma")
}

library(SummarizedExperiment)
library(readxl)
library(tidyverse)

# Carga de los datos de la hoja 1
original_data <- read_excel("TIO2+PTYR-human-MSS+MSIvsPD.XLSX")
head(original_data)

# Carga de los datos de la hoja 1
targets <- read_excel("TIO2+PTYR-human-MSS+MSIvsPD.XLSX", sheet = 2)
head(targets)

# Nombres de las columnas
colnames(original_data)

# Seleccionar las columnas de abundancia y almacenarlas en una variable
abundancias <- original_data %>% select(5:16)
abundancias <- as.data.frame(abundancias)

# Asignarle nuevos nombres a las columnas de la matriz de abundancia
accesion_rownames <- make.names(original_data$Accession, unique = TRUE)
rownames(abundancias) <- accesion_rownames
head(abundancias)

# Crear el DataFrame de metadatos de las muestras
sample_info <- DataFrame(
  SampleID = colnames(abundancias),
  Group = c(rep("MSS", 6), rep("PD", 6)),
  Replicate = rep(c(1, 2), times = 6)
)

# Crear el DataFrame con metadatos de los fosfopéptidos
row_data <- DataFrame(
  SequenceModifications = original_data$SequenceModifications,
  Description = original_data$Description,
  row.names = rownames(abundancias)
)

# Crear el objeto SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(abundancias)),  # Matriz de abundancia
  colData = sample_info,  # Metadata de las muestras
  rowData = row_data  # Metadata de los fosfopéptidos
)

# Resumen estadístico de la matriz de abundancia
summary(assay(se, "counts"))

# Convertir los datos de abundancia a formato largo
abundancias_long <- abundancias %>% 
  pivot_longer(cols = everything(), names_to = "Muestra", values_to = "Abundancia")

# Crear el histograma para cada muestra
ggplot(abundancias_long, aes(x = Abundancia)) +
  geom_histogram(bins = 20, fill = "red", color = "black") +
  facet_wrap(~ Muestra, scales = "free_x") +  
  labs(title = "Distribución de abundancias en cada muestra",
       x = "Abundancia",
       y = "Frecuencia") +
  theme_minimal()

# Visualizar distribuciones de abundancia en escala logarítmica
boxplot(log10(assay(se, "counts") + 1), 
        main = "Distribución de abundancias en escala log10",
        xlab = "Muestras",
        ylab = "Log10(Abundancia + 1)",
        las = 2,  
        col = "lightblue",  
        border = "darkblue")  

# Extraer la matriz de abundancia del objeto SummarizedExperiment
abundancias <- assay(se, "counts")

# Convertir la matriz de abundancia a formato largo y aplicar transformación logarítmica
logDat <- as.data.frame(abundancias) %>%
  pivot_longer(cols = everything(), names_to = "Muestra", values_to = "Abundancia") %>%
  mutate(log_abundance = log10(Abundancia + 1))

# Extraer información de Grupo y Réplica desde el nombre de la muestra
covs <- str_split(logDat$Muestra, "_", simplify = TRUE)
colnames(covs) <- c("Sample", "Replicate", "Group")

# Añadir las nuevas columnas (Sample, Replicate, Group) al dataframe logDat
logDat2 <- cbind(logDat, covs)

# Crear el boxplot con ggplot2 usando el grupo y la réplica como variables estéticas
ggplot(logDat2, aes(x = Muestra, y = log_abundance, fill = Group, colour = Replicate)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Distribución de abundancias de fosfopéptidos (escala log10)",
       x = "Muestras",
       y = "Log10(Abundancia + 1)")

# Análisis de Componentes Principales
source("https://raw.githubusercontent.com/uebvhir/UEB_PCA/master/UEB_plotPCA3.R")
plotPCA3(datos=as.matrix(log10(abundancias+1)), labels=colnames(abundancias), 
         factor=targets$Phenotype, title ="Phosphoproteomic data",
         scale=FALSE, colores=1:2, size = 3.5, glineas = 2.5)

# Definir grupos como factor y crear matriz de diseño
targets <- as.data.frame(targets)  # Asegúrate de que 'targets' esté correctamente cargado
groups <- as.factor(targets$Phenotype)  # Ajusta 'Phenotype' según la columna que contiene los grupos MSS y PD

# Crear la matriz de diseño para Limma
designMat <- model.matrix(~ -1 + groups)
print(designMat)

# Calcular la correlación media entre réplicas técnicas
if (!require(statmod)) install.packages("statmod")
dupcor <- duplicateCorrelation(abundancias, designMat, block=targets$Individual)
print(dupcor$consensus.correlation)

# Crear la matriz de contraste para comparar PD y MSS
contMat <- makeContrasts(mainEff = groupsPD - groupsMSS, levels = designMat)
show(contMat)

fit <- lmFit(abundancias, designMat, block=targets$Individual, correlation=dupcor$consensus)
fit2 <- contrasts.fit(fit, contMat)
fit2 <- eBayes(fit2)
results <- topTableF(fit2, adjust="BH", number=nrow(abundancias))
head(results)

# Graficar el Volcano Plot
volcanoplot(fit2, highlight=10, names=rownames(abundancias), cex=0.75, xlim=c(-1e+06, 1e+06))

# Guardar el objeto SummarizedExperiment como .RDS
saveRDS(se, file = "summarize_data/summarized_experiment.rds")

# Guardar el objeto en formato .RData
save(se, file = "summarize_data/summarized_experiment.RData")

# Guardar los datos en forma de texto
# Exportar la matriz de datos de abundancia
write.csv(as.data.frame(assay(se, "counts")), "datos_texto/abundancias.csv", row.names = TRUE)

# Exportar los metadatos de las muestras
write.csv(as.data.frame(colData(se)), "datos_texto/sample_info.csv", row.names = TRUE)

# Exportar los metadatos de los fosfopéptidos
write.csv(as.data.frame(rowData(se)), "datos_texto/row_info.csv", row.names = TRUE)


