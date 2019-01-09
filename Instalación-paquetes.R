#Instalación de paquetes

install.packages(c("shiny", "shinydashboard", "utils"), dep = TRUE)

source("http://bioconductor.org/biocLite.R")
biocLite(c("clusterProfiler", "DOSE", "GOSemSim", "ChIPseeker", "pathview", "org.Hs.eg.db", "org.Mm.eg.db",
           "org.Rn.eg.db", "org.Sc.sgd.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.At.eg.db", "org.Bt.eg.db",
           "org.Ce.eg.db", "org.Gg.eg.db", "org.Cf.eg.db", "org.Ss.eg.db", "org.Mmu.eg.db", "org.Eck12.eg.db",
           "org.Ag.eg.db", "org.Xl.eg.db", "org.Pt.eg.db", "org.Pf.plasmo.db", "org.EcSakai.eg.db",
           "TxDb.Hsapiens.UCSC.hg19.knownGene"))
