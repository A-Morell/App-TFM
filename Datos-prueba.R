# Datos de prueba

# Para crear una lista de genes
data(geneList, package="DOSE")
genes <- names(geneList)[abs(geneList) > 2]

write.table(genes, file = "genes.txt", sep = "\t", col.names = FALSE, row.names = FALSE)

# Para crear el universo de genes
write.table(geneList, file = "universe.txt", sep = "\t", col.names = FALSE, row.names = FALSE)

# Para crear un lista de picos
files_NGS <- getSampleFiles()
peaks <- readPeakFile(files_NGS[[4]])

write.table(peaks, file = "peaks.txt", sep = "\t", dec = ".", col.names = FALSE, row.names = FALSE)
