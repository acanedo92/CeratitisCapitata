#-----------------------------------------------------------------------------#
#-------------------------CERATITIS CAPITATA COLIMA CODIGO----------------------------#
#-----------------------------------------------------------------------------#
library("hierfstat")
library("adegenet")
library("corrplot")
library("poppr")
library("PopGenReport")


##Leer el archivo tipo genealex
  ceratitis.ind <- read.genalex(here::here("00.Archivos/3.Deschepper_RemocionAlelos.csv"), ploidy = 2, geo = FALSE,region = FALSE, sep = "\t", recode = FALSE)
ceratitis.genid <- genclone2genind(ceratitis.ind) 

##PopGeneReport

popgenreport(ceratitis.genid , mk.counts = TRUE, mk.locihz = TRUE, mk.hwe = TRUE,
             mk.fst = TRUE, mk.pcoa = TRUE,              
             mk.allele.dist = TRUE, mk.null.all = TRUE, mk.allel.rich = TRUE,
             mk.differ.stats = TRUE, mk.custom = FALSE, fname = "RemocionAlelos",
             foldername = "Reporte", 
              path.pgr = "/home/usuario/Escritorio/CeratitisDeschepper/03.Resultados_RemocionAlelos/", 
             mk.Rcode = FALSE,
             mk.complete = TRUE, mk.pdf = FALSE)


#Diversidad Genetica


Cols <- c("yellow", "chocolate4", "red4", "black", "blue")
color_pallete_function <- colorRampPalette(colors = Cols, space="Lab")
num_colors <- nlevels(pop(ceratitis.ind))
Colores <- color_pallete_function((num_colors))


##Plot
BasicPop <- basic.stats(ceratitis.ind)
Alelos.r <- allelic.richness(ceratitis.ind)

Sitios <- levels(pop(ceratitis.ind))

png("03.Resultados_RemocionAlelos/DiversidadGenetica.png")

par(mfrow=c(1,3))     
barplot(colMeans(BasicPop$Hs, na.rm = T), col = Colores, names.arg = Sitios, las=2, axis.lty = 6, axisnames = T, ylab = "H. Esperada")
barplot(colMeans(BasicPop$Ho, na.rm = T), col = Colores, names.arg = Sitios, las=2, axis.lty = 6, axisnames = T, ylab = "H. observada")
barplot(colSums(Alelos.r$Ar),  las=2, col = Colores, ylab = "N. de alelos")

dev.off()

## Estimar Fst Pareada
ceratitis.fstat <- genind2hierfstat(ceratitis.ind)
Fstat.metrics <- pairwise.neifst(ceratitis.fstat)


library("corrplot")

png("03.Resultados_RemocionAlelos/CorrelacionFst.png")
corrplot(as.matrix(Fstat.metrics), is.corr = F, type = "lower", diag=F)      
dev.off()


## Hardy-Weinberg
library("lattice")
ceratitis.hwe.pop <- hw.test(ceratitis.genid ,pop=NULL,permut=FALSE,nsim=1999,hide.NA=TRUE,res.type=c("full","matrix"))
ceratitis.hwe.full <- (ceratitis.hwe.pop <- seppop(ceratitis.genid) %>% lapply(hw.test, B = 0))

ceratitis.hwe.mat <- sapply(ceratitis.hwe.pop, "[", i = TRUE, j = 3) # Take the third column with all rows
alpha  <- 0.05
ceratitis.hwe.mat[ceratitis.hwe.mat > alpha] <- 1


png("03.Resultados_RemocionAlelos/Equilibrio.png")
levelplot(t(ceratitis.hwe.mat), main="Equilibrio HW")
dev.off()

## Generar una matriz de distancias entre las poblaciones
library(usedist)

png("03.Resultados_RemocionAlelos/Dendograma.png")
nei.1 <- read.csv("03.Resultados_RemocionAlelos/Reporte/RemocionAlelos-pairwise_Gst_Nei.csv")
nei.2 <- dist(nei.1)
nm <- levels(pop(ceratitis.ind))
nei.3 <- dist_setNames(nei.2, nm)
nei.4 <- hclust(nei.3, method="average")
plot(as.dendrogram(nei.4),  horiz = F, ylab="Distancia")
dev.off()




