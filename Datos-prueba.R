# Datos de prueba

# Para crear una lista de genes
data(geneList, package="DOSE")
genes <- names(geneList)[abs(geneList) > 2]

write.table(genes, file = "genes.txt", sep = "\t", dec = ".", col.names = FALSE, row.names = FALSE)


# Para crear un lista de picos
files_NGS <- getSampleFiles()
peak <- readPeakFile(files_NGS[[4]])
genes_NGS <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)

write.table(genes_NGS, file = "genes_NGS.txt", sep = "\t", dec = ".", col.names = FALSE, row.names = FALSE)
