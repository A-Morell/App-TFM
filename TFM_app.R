# Herramienta en línea para la realización de análisis de significación biológica
#
# Máster en Bioinformática y Bioestadística - UOC
#
# Antonio Morell Bennasser








## Carga de los paquetes

require(shiny)
require(shinydashboard)
require(utils)

require(clusterProfiler)
require(DOSE)
require(GOSemSim)
require(ChIPseeker)
require(pathview)

require(org.Hs.eg.db)
require(org.Mm.eg.db)
require(org.Rn.eg.db)
require(org.Sc.sgd.db)
require(org.Dm.eg.db)
require(org.Dr.eg.db)
require(org.At.eg.db)
require(org.Bt.eg.db)
require(org.Ce.eg.db)
require(org.Gg.eg.db)
require(org.Cf.eg.db)
require(org.Ss.eg.db)
require(org.Mmu.eg.db)
require(org.Eck12.eg.db)
require(org.Ag.eg.db)
require(org.Xl.eg.db)
require(org.Pt.eg.db)
require(org.Pf.plasmo.db)
require(org.EcSakai.eg.db)

require(TxDb.Hsapiens.UCSC.hg19.knownGene)





## UI & Server

# Definimos la UI
ui <- dashboardPage(
  
  dashboardHeader(title = "Análisis de significación biológica", titleWidth = 350),
     
     # Generamos la barra lateral
     dashboardSidebar(width = 350,
       
       # Definimos la barra lateral y los menús              
       sidebarMenu(
           
           menuItem(text = "Información", tabName = "info", icon = icon("fas fa-info-circle")),
           menuItem(text = "Carga de datos", tabName = "upload", icon = icon("upload")),
           menuItem(text = "GO Analysis", tabName = "goanalysis", icon = icon("fas fa-globe"),
                    menuSubItem(text = "GO Gene Enrichment Analysis (con enrichGO)", tabName = "GOgea"),
                    menuSubItem(text = "GO Gene Set Enrichment Analysis (con gseGO)", tabName = "GOgsea"),
                    menuSubItem(text = "GO Semantic Similarity Analysis (con mgoSim)", tabName = "GOssa")
                    ),
           menuItem(text = "KEGG Analysis", tabName = "kegganalysis", icon = icon("fas fa-globe"),
                    menuSubItem(text = "KEGG Gene Enrichment Analysis (con enrichKEGG)", tabName = "KEGGgea"),
                    menuSubItem(text = "KEGG Gene Set Enrichment Analysis (con gseKEGG)", tabName = "KEGGgsea")
                    ),
           menuItem(text = "NGS Analysis", tabName = "ngsanalysis", icon = icon("fas fa-globe"),
                    menuSubItem(text = "Preparación de los datos", tabName = "NGSdata")
           )
         )
       ),
           
     # Generamos los paneles y definimos cada apartado
     dashboardBody(
       tabItems(
         tabItem(tabName = "info",
                 p("Los análisis de significación biológica o análisis funcionales proporcionan resultados
                   dentro de un contexto biológico permitiendo interpretar resultados de análisis de expresión
                   diferencial de una forma más clara y concisa, relacionando las listas de genes
                   diferencialmente expresados obtenidas en estos análisis con bases de datos de anotaciones
                   funcionales. Esta herramienta en línea permite realizar los principales análisis de
                   significación biológica así como otros análisis complementarios de interés. Seguidamente se
                   describen las distintas partes de esta herramienta en línea:"),
                 tags$ol(
                   tags$li(p(strong("Carga de los datos: "), " Cargar el conjunto de genes de interés contenidos
                             en un archivo en formato .txt con una sola columna que disponga de una cabecera y
                             en cuyas filas se encuentren los IDs de los genes. Cargar el universo de genes
                             (opcional). Cargar el conjunto de datos de recuentos. Cargar los conjuntos de
                             términos de GO para el análisis de similitudes entre términos.")
                           ),
                   tags$li(p(strong("GO analysis: "), " Análisis de enriquecimiento (Gen Enrichment Analysis -
                             GEA y Gen Set Enrichment Analysis - GSEA) a partir de los conjuntos de genes de
                             interés utilizando las anotaciones de Gene Ontology (GO) para la obtención de
                             términos sobre-representados o sub-representados. También es posible realizar
                             análisis de similitudes semánticas (SSA) para estudiar la presencia o ausencia de
                             similitudes entre los aspectos funcionales relacionados entre dos o más términos de
                             GO.")
                           ),
                   tags$li(p(strong("KEGG analysis: "), " Análisis de enriquecimiento (GEA y GSEA) a partir de
                             los conjuntos de genes de interés utilizando las anotaciones de la Kyoto
                             Encyclopedia of Genes and Genomes (KEGG) para la obtención de términos
                             sobre-representados o sub-representados así como obtener las rutas KEGG
                             (KEGG-pathways) sobre-representadas.")
                   ),
                   tags$li(p(strong("NGS analysis: "), " Análisis de datos obtenidos a partir de Next Generation
                             Sequencing (NGS). En primer lugar se debe realizar una preparación de los datos,
                             que, posteriormente, deberán descargarse para poder aplicar los análisis anteiores
                             deseados.")
                   )
                 )
                 ),
         
         # Panel de carga de datos
         tabItem(tabName = "upload",
                 h2("Cargar los archivos de IDs de los genes"),
                 br(),
                 fileInput(inputId = "geneFile",
                           label = "Cargar los ID de los genes",
                           accept = c("text", ".txt"),
                           buttonLabel = "Buscar...",
                           placeholder = "Ningún archivo"),
                 fileInput(inputId = "universeFile",
                           label = "Cargar el universo de genes (opcional)",
                           accept = c("text", ".txt"),
                           buttonLabel = "Buscar...",
                           placeholder = "Ningún archivo"),
                 fileInput(inputId = "GO1",
                           label = "Conjunto 1 de términos de GO (GO SSA)",
                           accept = c("text", ".txt"),
                           buttonLabel = "Buscar...",
                           placeholder = "Ningún archivo"),
                 fileInput(inputId = "GO2",
                           label = "Conjunto 2 de términos de GO (GO SSA)",
                           accept = c("text", ".txt"),
                           buttonLabel = "Buscar...",
                           placeholder = "Ningún archivo"),
                 fileInput(inputId = "peakFile",
                           label = "Cargar el archivo de NGS",
                           accept = c("text", ".txt"),
                           buttonLabel = "Buscar...",
                           placeholder = "Ningún archivo")
                 ),
         
         # GO Analysis
         ## Panel para el análisis de enriquecimiento GEA de GO (enrichGO)
         tabItem(tabName = "GOgea",
                 h2("GO Gene Enrichment Analysis"),
                 fluidRow(
                   box(title = "Configurar los parámetros del análisis",
                       status = "warning",
                       selectInput(inputId = "OrgDb",
                                   label = "OrgDb",
                                   choices = c("org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Sc.sgd.db",
                                               "org.Dm.eg.db", "org.Dr.eg.db", "org.At.eg.db", "org.Bt.eg.db",
                                               "org.Ce.eg.db", "org.Gg.eg.db", "org.Cf.eg.db", "org.Ss.eg.db",
                                               "org.Mmu.eg.db", "org.Eck12.eg.db", "org.Ag.eg.db",
                                               "org.Xl.eg.db", "org.Pt.eg.db", "org.Pf.plasmo.db",
                                               "org.EcSakai.eg.db"),
                                   selected = "org.Hs.eg.db"),
                       selectInput(inputId = "keyType",
                                   label = "Keytype of input gene",
                                   choices = c(keytypes(org.Hs.eg.db),
                                               "MGI", "COMMON", "DESCRIPTION", "INTERPRO", "ORF", "SGD", "SMART"),
                                   selected = "ENTREZID"),
                       selectInput(inputId = "ont",
                                   label = "Subontologies",
                                   choices = c("MF", "BP", "CC", "ALL"),
                                   selected = "MF"),
                       numericInput(inputId = "pvalueCutoff",
                                    label = "p-value cutoff",
                                    value = 0.05,
                                    step = 0.01),
                       selectInput(inputId = "pAdjMethod",
                                   label = "p-value adjust method",
                                   choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                   selected = "BH"),
                       numericInput(inputId = "qvalueCutoff",
                                    label = "q-value cutoff",
                                    value = 0.2,
                                    step = 0.1),
                       numericInput(inputId = "minGSSize",
                                    label = "Minimal size of genes annotated by Ontology term for testing",
                                    value = 10,
                                    step = 1),
                       numericInput(inputId = "maxGSSize",
                                    label = "Maximal size of genes annotated for testing",
                                    value = 500,
                                    step = 1),
                       selectInput(inputId = "readable",
                                   label = "Readable: whether mapping gene ID to gene Name",
                                   choices = c("TRUE", "FALSE"),
                                   selected = "TRUE"),
                       selectInput(inputId = "pool",
                                   label = "Pool: if ont='ALL', whether pool 3 GO sub-ontologies",
                                   choices = c("TRUE", "FALSE"),
                                   selected = "FALSE"),
                       actionButton(inputId = "submit1",
                                    label = "Submit")
                       )
                   ),
                 fluidRow(
                   box(title = "Tabla de resultados de GO GEA",
                       status = "primary",
                       tableOutput(outputId = "enrichGO_results"),
                       downloadButton(outputId = "download_GOGEA_res", label = "Download")
                       ),
                   box(title = "Diagrama de puntos de GO GEA",
                       status = "success",
                       plotOutput(outputId = "dotplot_GOgea"),
                       downloadButton(outputId = "download_GOGEA_dot", label = "Download")
                       ),
                   box(title = "Mapa de enriquecimiento de GO GEA",
                       status = "danger",
                       plotOutput(outputId = "emapplot_GOgea"),
                       downloadButton(outputId = "download_GOGEA_map", label = "Download")
                       ),
                   box(title = "Representación de asociaciones",
                       status = "info",
                       plotOutput(outputId = "cnetplot_GOgea"),
                       downloadButton(outputId = "download_GOGEA_cnet", label = "Download")
                       )
                 )
         ),
         
         ## Panel para el análisis de enriquecimiento GSEA de GO (gseGO)
         tabItem(tabName = "GOgsea",
                 h2("GO Gene Set Enrichment Analysis"),
                 fluidRow(
                   box(title = "Configurar los parámetros del análisis",
                       status = "warning",
                       selectInput(inputId = "ont_gse",
                                   label = "Subontologies",
                                   choices = c("BP", "MF", "CC", "GO"),
                                   selected = "MF"),
                       selectInput(inputId = "OrgDb_gse",
                                   label = "OrgDb",
                                   choices = c("org.Hs.eg.db", "org.Mm.eg.db", "org.Sc.sgd.db"),
                                   selected = "org.Hs.eg.db"),
                       textInput(inputId = "keyType_gse",
                                 label = "Keytype of input gene",
                                 value = "ENTREZID"),
                       numericInput(inputId = "exponent",
                                    label = "Exponent: weight of each step",
                                    value = 1,
                                    step = 0.5),
                       numericInput(inputId = "nPerm",
                                    label = "Permutation numbers",
                                    value = 1000,
                                    step = 100),
                       numericInput(inputId = "minGSSize_gse",
                                    label = "Minimal size of each geneSet for analyzing",
                                    value = 10,
                                    step = 1),
                       numericInput(inputId = "maxGSSize_gse",
                                    label = "Maximal size of genes annotated for testing",
                                    value = 500,
                                    step = 1),
                       numericInput(inputId = "pvalueCutoff_gse",
                                    label = "p-value cutoff",
                                    value = 0.05,
                                    step = 0.01),
                       selectInput(inputId = "pAdjMethod_gse",
                                   label = "p-value adjust method",
                                   choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                   selected = "BH"),
                       selectInput(inputId = "verbose",
                                   label = "Print message or not",
                                   choices = c("TRUE", "FALSE"),
                                   selected = "TRUE"),
                       selectInput(inputId = "seed",
                                   label = "Seed",
                                   choices = c("TRUE", "FALSE"),
                                   selected = "FALSE"),
                       selectInput(inputId = "by",
                                   label = "By one of 'fgsea' or 'DOSE'",
                                   choices = c("fgsea", "DOSE"),
                                   selected = "fgsea"),
                       actionButton(inputId = "submit2",
                                    label = "Submit")
                   )
                 ),
                 fluidRow(
                   box(title = "Tabla de resultados de GO GSEA",
                       status = "primary",
                       tableOutput(outputId = "gseGO_results"),
                       downloadButton(outputId = "download_GOGSEA_res", label = "Download")
                   ),
                   box(title = "Diagrama de puntos de GO GSEA",
                       status = "success",
                       plotOutput(outputId = "dotplot_GOgsea"),
                       downloadButton(outputId = "download_GOGSEA_dot", label = "Download")
                       ),
                   box(title = "Mapa de enriquecimiento de GO GSEA",
                       status = "danger",
                       plotOutput(outputId = "emapplot_GOgsea"),
                       downloadButton(outputId = "download_GOGSEA_map", label = "Download")
                   ),
                   box(title = "Representación de asociaciones",
                       status = "info",
                       plotOutput(outputId = "cnetplot_GOgsea"),
                       downloadButton(outputId = "download_GOGSEA_cnet", label = "Download")
                   ),
                   box(title = "Puntuación de GSEA y asociación de fenotipo",
                       status = "warning",
                       plotOutput(outputId = "gseaplot_GOgsea"),
                       downloadButton(outputId = "download_GOGSEA_fen", label = "Download")
                   )
                 )
         ),
         
         ## Panel para el test de similitudes semánticas de GO (mgoSim)
         tabItem(tabName = "GOssa",
                 h2("GO Semantic Similarity Analysis"),
                 fluidRow(
                   box(title = "Configurar los parámetros del análisis",
                       status = "warning",
                       selectInput(inputId = "measure",
                                   label = "Choose one measure method",
                                   choices = c("Resnik", "Lin", "Rel", "Jiang", "Wang"),
                                   selected = "Wang"),
                       selectInput(inputId = "combine",
                                   label = "Choose one combine method",
                                   choices = c("max", "avg", "rcmax", "BMA"),
                                   selected = "BMA"),
                       actionButton(inputId = "submit3",
                                    label = "Submit")
                   )
                 ),
                 fluidRow(
                   box(title = "Resultados de GO SSA",
                       status = "primary",
                       plotOutput(outputId = "mgoSim_results"),
                       downloadButton(outputId = "download_GOSSA_res", label = "Download")
                   )
                 )
         ),
         
         # KEEG Analysis
         ## Panel para el test de sobre-representación de KEGG (enrichKEGG)
         tabItem(tabName = "KEGGgea",
                 h2("KEGG Gene Enrichment Analysis"),
                 fluidRow(
                   box(title = "Configurar los parámetros del análisis",
                       status = "warning",
                       textInput(inputId = "organism",
                                 label = "Organism from 'https://www.genome.jp/kegg/catalog/org_list.html'",
                                 value = "hsa"),
                       textInput(inputId = "keyType",
                                 label = "Keytype of input gene",
                                 value = "kegg"),
                       numericInput(inputId = "pvalueCutoff",
                                    label = "p-value cutoff",
                                    value = 0.05,
                                    step = 0.01),
                       selectInput(inputId = "pAdjMethod",
                                   label = "p-value adjust method",
                                   choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                   selected = "BH"),
                       numericInput(inputId = "minGSSize",
                                    label = "minimal size of genes annotated by Ontology term for testing",
                                    value = 10,
                                    step = 1),
                       numericInput(inputId = "maxGSSize",
                                    label = "maximal size of genes annotated for testing",
                                    value = 500,
                                    step = 1),
                       numericInput(inputId = "qvalueCutoff",
                                    label = "q-value cutoff",
                                    value = 0.2,
                                    step = 0.1),
                       actionButton(inputId = "submit4",
                                    label = "Submit")
                   )
                 ),
                 fluidRow(
                   box(title = "Tabla de resultados de KEGG GEA",
                       status = "primary",
                       verbatimTextOutput(outputId = "enrichKEGG_results"),
                       downloadButton(outputId = "download_KEGGGEA_res", label = "Download")
                   ),
                   box(title = "Diagrama de puntos de KEGG GEA",
                       status = "success",
                       plotOutput(outputId = "dotplot_KEGGgea"),
                       downloadButton(outputId = "download_KEGGGEA_dot", label = "Download")
                   ),
                   box(title = "Mapa de enriquecimiento de KEGG GEA",
                       status = "danger",
                       plotOutput(outputId = "emapplot_KEGGgea"),
                       downloadButton(outputId = "download_KEGGGEA_map", label = "Download")
                   ),
                   box(title = "Representación de asociaciones",
                       status = "info",
                       plotOutput(outputId = "cnetplot_KEGGgea"),
                       downloadButton(outputId = "download_KEGGGEA_cnet", label = "Download")
                   ),
                   box(title = "Representación de rutas KEGG",
                       status = "warning",
                       textInput(inputId = "pathwayId",
                                 label = "The KEGG pathway ID"),
                       textInput(inputId = "species",
                                 label = "Target species",
                                 value = "hsa"),
                       actionButton(inputId = "submit5",
                                    label = "Submit")
                       )
                 )
         ),
         
         ## Panel para el análisis de enriquecimiento GSEA de KEGG (gseKEGG)
         tabItem(tabName = "KEGGgsea",
                 h2("KEGG Gene Set Enrichment Analysis"),
                 fluidRow(
                   box(title = "Configurar los parámetros del análisis",
                       status = "warning",
                       textInput(inputId = "organism_gse",
                                 label = "Organism from 'https://www.genome.jp/kegg/catalog/org_list.html'",
                                 value = "hsa"),
                       selectInput(inputId = "keyType_KEGGgse",
                                 label = "Keytype of input gene",
                                 choices = c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot"),
                                 selected = "kegg"),
                       numericInput(inputId = "exponent_gse",
                                    label = "Exponent: weight of each step",
                                    value = 1,
                                    step = 0.5),
                       numericInput(inputId = "nPerm_gse",
                                    label = "Permutation numbers",
                                    value = 1000,
                                    step = 100),
                       numericInput(inputId = "minGSSize_KEGGgse",
                                    label = "Minimal size of each geneSet for analyzing",
                                    value = 10,
                                    step = 1),
                       numericInput(inputId = "maxGSSize_KEGGgse",
                                    label = "Maximal size of genes annotated for testing",
                                    value = 500,
                                    step = 1),
                       numericInput(inputId = "pvalueCutoff_KEGGgse",
                                    label = "p-value cutoff",
                                    value = 0.05,
                                    step = 0.01),
                       selectInput(inputId = "pAdjMethod_KEGGgse",
                                   label = "p-value adjust method",
                                   choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                   selected = "BH"),
                       selectInput(inputId = "verbose_gse",
                                   label = "Print message or not",
                                   choices = c("TRUE", "FALSE"),
                                   selected = "TRUE"),
                       selectInput(inputId = "use_internal_data",
                                   label = "Use KEGG.db or latest online KEGG data",
                                   choices = c("TRUE", "FALSE"),
                                   selected = "FALSE"),
                       selectInput(inputId = "seed_gse",
                                   label = "Seed",
                                   choices = c("TRUE", "FALSE"),
                                   selected = "FALSE"),
                       selectInput(inputId = "by_gsea",
                                   label = "By one of 'fgsea' or 'DOSE'",
                                   choices = c("fgsea", "DOSE"),
                                   selected = "fgsea"),
                       actionButton(inputId = "submit6",
                                    label = "Submit")
                   )
                 ),
                 fluidRow(
                   box(title = "Tabla de resultados de KEGG GSEA",
                       status = "primary",
                       tableOutput(outputId = "gseKEGG_results"),
                       downloadButton(outputId = "download_KEGGGSEA_res", label = "Download")
                   ),
                   box(title = "Diagrama de puntos de KEGG GSEA",
                       status = "success",
                       plotOutput(outputId = "dotplot_KEGGgsea"),
                       downloadButton(outputId = "download_KEGGGSEA_dot", label = "Download")
                   ),
                   box(title = "Mapa de enriquecimiento de KEGG GSEA",
                       status = "danger",
                       plotOutput(outputId = "emapplot_KEGGgsea"),
                       downloadButton(outputId = "download_KEGGGSEA_map", label = "Download")
                   ),
                   box(title = "Representación de asociaciones",
                       status = "info",
                       plotOutput(outputId = "cnetplot_KEGGgsea"),
                       downloadButton(outputId = "download_KEGGGSEA_cnet", label = "Download")
                   ),
                   box(title = "Puntuación de GSEA y asociación de fenotipo",
                       status = "warning",
                       plotOutput(outputId = "gseaplot_KEGGgsea"),
                       downloadButton(outputId = "download_KEGGGSEA_fen", label = "Download")
                   )
                 )
         ),
         
         ## Panel para la preparación de los datos de NGS
         tabItem(tabName = "NGSdata",
                 h2("Análisis de datos de NGS. Preparación de los datos"),
                 fluidRow(
                   box(title = "Preparación de los datos de NGS",
                       status = "warning",
                       numericInput(inputId = "tssRegion_down",
                                    label = "TSS region down",
                                    value = -1000,
                                    step = 500),
                       numericInput(inputId = "tssRegion_up",
                                    label = "TSS region up",
                                    value = 1000,
                                    step = 500),
                       numericInput(inputId = "flankDistance",
                                    label = "Flanking search radius",
                                    value = 3000,
                                    step = 500),
                       textInput(inputId = "TxDb",
                                    label = "TranscriptDb object",
                                    value = "TxDb.Hsapiens.UCSC.hg19.knownGene"),
                       selectInput(inputId = "sameStrand",
                                   label = "Whether find nearest/overlap gene in the same strand",
                                   choices = c("TRUE", "FALSE"),
                                   selected = "FALSE"),
                       actionButton(inputId = "submit7",
                                    label = "Submit"),
                       downloadButton(outputId = "downloadGeneFile", label = "Download")
                   )
                 )
         )
       )
     )
  )
         


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Cargamos los datos del usuario
  ## Carga del archivo de datos de los ID de los genes
  genesID <- reactive({
    genesFile <- input$genesFile
    if (!is.null(genesFile)) {
      genesID <- read.table(genesFile$datapath)
      genesID <- as.character(genesID[, 1])
    }
    return(genesID)
  })
  
  ## Carga del archivo de datos del universo de genes
  universe <- reactive({
    universeFile <- input$universeFile
    if (!is.null(universeFile)) {
      universe <- read.table(universeFile$datapath)
      universe <- as.character(universe[, 1])
    }
    return(universe)
  })
  
  ## Carga de los archivos de datos de los términos GO para la función mgoSim
  GOTerms1 <- reactive({
    GOFile1 <- input$GO1
    if (!is.null(GOFile1)) {
      GOTerms1 <- read.table(GOFile1$datapath)
      GOTerms1 <- as.character(GOTerms1[, 1])
    }
    return(GOTerms1)
  })
  
  GOTerms2 <- reactive({
    GOFile2 <- input$GO2
    if (!is.null(GOFile2)) {
      GOTerms2 <- read.table(GOFile2$datapath)
      GOTerms2 <- as.character(GOTerms2[, 1])
    }
    return(GOTerms2)
  })
  
  ## Carga de los archivos de datos de NGS
  NGSpeaks <- reactive({
    NGSpeakFile <- input$peakFile
    if (!is.null(NGSpeakFile)) {
      NGSpeaks <- read.table(NGSpeakFile$datapath)
      NGSpeaks <- as.character(NGSpeaks[, 1])
    }
    return(NGSpeaks)
  })
  
  # GO Analysis
  ## Output de la función enrichGO
  e_GO <- eventReactive(input$submit1, {
    e_GO <- enrichGO(
      gene = genesID(), OrgDb = input$OrgDb, keyType = input$keyType, ont = input$ont,
      pvalueCutoff = input$pvalueCutoff, pAdjustMethod = input$pAdjMethod,
      universe = universe(), qvalueCutoff = input$qvalCutoff,
      minGSSize = input$minGSSize, maxGSSize = input$maxGSSize,
      readable = input$readable, pool = input$pool
      )
    })
  
  ### Resultados de la función enrichGO
  results_enrichGO <- reactive({
    if (!is.null(e_GO())) {
      results_enrichGO <- e_GO@result
    } else {NULL}
  })
  
  ### Output de la tabla con los resultados de la función enrichGO
  output$enrichGO_results <- renderTable(
    results_enrichGO()
    )
  
  #### Descarga de la tabla de resultados de enrichGO
  output$download_GOGEA_res <- downloadHandler(
    filename = "GO_GEA_results.txt",
    content = function(file) {
      write.csv(results_enrichGO(), file, sep = ",", row.names = FALSE, col.names = TRUE)
    }
  )
  
  ### Output con el diagrama de puntos
  output$dotplot_GOgea <- renderPlot(
    if (!is.null(e_GO())) {
      DOSE::dotplot(e_GO())
    }
  )
  
  #### Descarga del diagrama de puntos de GO GEA
  output$download_GOGEA_dot <- downloadHandler(
    filename = "GO_GEA_dotplot.txt",
    content = function(file) {
      ggsave(file, DOSE::dotplot(e_GO()))
    }
  )
  
  ### Output con el mapa de enriquecimiento
  output$emapplot_GOgea <- renderPlot(
    if (!is.null(e_GO())) {
      DOSE::emapplot(e_GO())
    }
  )
  
  #### Descarga del mapa de enriquecimiento
  output$download_GOGEA_map <- downloadHandler(
    filename = "GO_GEA_emapplot.txt",
    content = function(file) {
      ggsave(file, DOSE::emapplot(e_GO()))
    }
  )
  
  ### Output con la representación de las asociaciones
  output$cnetplot_GOgea <- renderPlot(
    if (!is.null(e_GO())) {
      DOSE::cnetplot(e_GO())
    }
  )
  
  #### Descarga de la representación de asociaciones
  output$download_GOGEA_cnet <- downloadHandler(
    filename = "GO_GEA_cnetplot.txt",
    content = function(file) {
      ggsave(file, DOSE::cnetplot(e_GO()))
    }
  )
  
  ## Output de la función gseGO
  gse_GO <- eventReactive(input$submit2, {
    gse_GO <- gseGO(
      genesID(), ont = input$ont_gse, OrgDb = input$OrgDb_gse, keyType = input$keyType_gse,
      exponent = input$exponent, nPerm = input$nPerm, minGSSize = input$minGSSize_gse,
      maxGSSize = input$maxGSSize_gse, pvalueCutoff = input$pvalueCutoff_gse,
      pAdjustMethod = input$pAdjMethod_gse, verbose = input$verbose, seed = input$seed, by = input$by
    )
  })
  
  ### Resultados de la función gseGO
  results_gseGO <- reactive({
    if (!is.null(gse_GO())) {
      results_gseGO <- gse_GO@result
    } else {NULL}
  })
  
  ### Output de la tabla con los resultados de la función gseGO
  output$gseGO_results <- renderTable(
    results_gseGO()
  )
  
  #### Descarga de la tabla de resultados de gseGO
  output$download_GOGSEA_res <- downloadHandler(
    filename = "GO_GSEA_results.txt",
    content = function(file) {
      write.csv(results_gseGO(), file, sep = ",", row.names = FALSE, col.names = TRUE)
    }
  )
  
  ### Output con el diagrama de puntos
  output$dotplot_GOgsea <- renderPlot(
    if (!is.null(gse_GO())) {
      DOSE::dotplot(gse_GO())
    }
  )
  
  #### Descarga del diagrama de puntos de GO GSEA
  output$download_GOGSEA_dot <- downloadHandler(
    filename = "GO_GSEA_dotplot.txt",
    content = function(file) {
      ggsave(file, DOSE::dotplot(gse_GO()))
    }
  )
  
  ### Output con el mapa de enriquecimiento
  output$emapplot_GOgsea <- renderPlot(
    if (!is.null(gse_GO())) {
      DOSE::emapplot(gse_GO())
    }
  )
  
  #### Descarga del mapa de enriquecimiento
  output$download_GOGSEA_map <- downloadHandler(
    filename = "GO_GSEA_emapplot.txt",
    content = function(file) {
      ggsave(file, DOSE::emapplot(gse_GO()))
    }
  )
  
  ### Output con la representación de las asociaciones
  output$cnetplot_GOgsea <- renderPlot(
    if (!is.null(gse_GO())) {
      DOSE::cnetplot(gse_GO())
    }
  )
  
  #### Descarga de la representación de las asociaciones
  output$download_GOGSEA_cnet <- downloadHandler(
    filename = "GO_GSEA_cnetplot.txt",
    content = function(file) {
      ggsave(file, DOSE::cnetplot(gse_GO()))
    }
  )
  
  ### Output con la puntuación de GSEA y asociación de fenotipo
  output$gseaplot_GOgsea <- renderPlot(
    if (!is.null(gse_GO())) {
      DOSE::gseaplot(gse_GO())
    }
  )
  
  #### Descarga de la puntuación de GSEA y asociación de fenotipo
  output$download_GOGSEA_fen <- downloadHandler(
    filename = "GO_GSEA_gseaplot.txt",
    content = function(file) {
      ggsave(file, DOSE::gseaplot(gse_GO()))
    }
  )
  
  ## Output de la función mgoSim
  ssa_GO <- eventReactive(input$submit3, {
    ssa_GO <- mgoSim(
      GO1 = GOTerms1(), GO2 = GOTerms2(), semData = (hsGO <- godata('org.Hs.eg.db', ont = "MF")),
      measure = input$measure, combine = input$combine
    )
  })
  
  ### Resultados de la función mgoSim
  results_mgoSimGO <- reactive({
    if (!is.null(ssa_GO())) {
      results_mgoSimGO <- ssa_GO@result
    } else {NULL}
  })
  
  ### Output de la tabla con los resultados de la función mgoSim
  output$mgoSim_results <- renderTable(
    results_mgoSimGO()
  )
  
  #### Descarga de la tabla de resultados de gseGO
  output$download_GOSSA_res <- downloadHandler(
    filename = "GO_SSA_results.txt",
    content = function(file) {
      write.csv(results_mgoSimGO(), file, sep = ",", row.names = FALSE, col.names = TRUE)
    }
  )
  
  # KEGG Analysis
  ## Output de la función enrichKEGG
  e_KEGG <- eventReactive(input$submit4, {
    e_KEGG <- enrichKEGG(
      gene = genesID(), organism = input$organism, keyType = input$keyType,
      pvalueCutoff = input$pvalueCutoff, pAdjustMethod = input$pAdjMethod,
      universe = universe(), minGSSize = input$minGSSize,
      maxGSSize = input$maxGSSize, qvalueCutoff = input$qvalCutoff
    )
  })
  
  ### Resultados de la función enrichKEGG
  results_enrichKEGG <- reactive({
    if (!is.null(e_KEGG())) {
      results_enrichKEGG <- e_KEGG@result
    } else {NULL}
  })
  
  ### Output de la tabla con los resultados de la función enrichKEGG
  output$enrichKEGG_results <- renderPrint(
    results_enrichKEGG()
  )
  
  #### Descarga de la tabla de resultados de enrichKEGG
  output$download_KEGGGEA_res <- downloadHandler(
    filename = "KEGG_GEA_results.txt",
    content = function(file) {
      write.csv(results_enrichKEGG(), file, sep = ",", row.names = FALSE, col.names = TRUE)
    }
  )
  
  ### Output con el diagrama de puntos
  output$dotplot_KEGGgea <- renderPlot(
    if (!is.null(e_KEGG())) {
      DOSE::dotplot(e_KEGG())
    }
  )
  
  #### Descarga del diagrama de puntos
  output$download_KEGGGEA_dot <- downloadHandler(
    filename = "KEGG_GEA_dotplot.txt",
    content = function(file) {
      ggsave(file, DOSE::dotplot(e_KEGG()))
    }
  )
  
  ### Output con el mapa de enriquecimiento
  output$emapplot_KEGGgea <- renderPlot(
    if (!is.null(e_KEGG())) {
      DOSE::emapplot(e_KEGG())
    }
  )
  
  #### Descarga del mapa de enriquecimiento
  output$download_KEGGGEA_map <- downloadHandler(
    filename = "KEGG_GEA_emapplot.txt",
    content = function(file) {
      ggsave(file, DOSE::emapplot(e_KEGG()))
    }
  )
  
  ### Output con la representación de las asociaciones
  output$cnetplot_KEGGgea <- renderPlot(
    if (!is.null(e_KEGG())) {
      DOSE::cnetplot(e_KEGG())
    }
  )
  
  #### Descarga de la representación de las asociaciones
  output$download_KEGGGEA_cnet <- downloadHandler(
    filename = "KEGG_GEA_cnetplot.txt",
    content = function(file) {
      ggsave(file, DOSE::cnetplot(e_KEGG()))
    }
  )

  ## Output de la función gseKEGG
  gse_KEGG <- eventReactive(input$submit6, {
    gse_KEGG <- gseKEGG(
      genesID(), organism = input$organism_gse, keyType = input$keyType_KEGGgse,
      exponent = input$exponent_gse, nPerm = input$nPerm_gse, minGSSize = input$minGSSize_KEGGgse,
      maxGSSize = input$maxGSSize_KEGGgse, pvalueCutoff = input$pvalueCutoff_KEGGgse,
      pAdjustMethod = input$pAdjMethod_KEGGgse, verbose = input$verbose_gse,
      use_internal_data = input$use_internal_data, seed = input$seed_gse, by = input$by_gse
    )
  })
  
  ### Resultados de la función gseKEGG
  results_gseKEGG <- reactive({
    if (!is.null(gse_KEGG())) {
      results_gseKEGG <- gse_KEGG@result
    } else {NULL}
  })
  
  ### Output de la tabla con los resultados de la función gseKEGG
  output$gseKEGG_results <- renderTable(
    results_gseKEGG()
  )
  
  #### Descarga de la tabla de resultados de gseKEGG
  output$download_KEGGGSEA_res <- downloadHandler(
    filename = "KEGG_GSEA_results.txt",
    content = function(file) {
      write.csv(results_gseKEGG(), file, sep = ",", row.names = FALSE, col.names = TRUE)
    }
  )
  
  ### Output con el diagrama de puntos
  output$dotplot_KEGGgsea <- renderPlot(
    if (!is.null(gse_KEGG())) {
      DOSE::dotplot(gse_KEGG())
    }
  )
  
  #### Descarga del diagrama de puntos
  output$download_KEGGGSEA_dot <- downloadHandler(
    filename = "KEGG_GSEA_dotplot.txt",
    content = function(file) {
      ggsave(file, DOSE::dotplot(gse_KEGG()))
    }
  )
  
  ### Output con el mapa de enriquecimiento
  output$emapplot_KEGGgsea <- renderPlot(
    if (!is.null(gse_KEGG())) {
      DOSE::emapplot(gse_KEGG())
    }
  )
  
  #### Descarga del mapa de enriquecimiento
  output$download_KEGGGSEA_map <- downloadHandler(
    filename = "KEGG_GSEA_emapplot.txt",
    content = function(file) {
      ggsave(file, DOSE::emapplot(gse_KEGG()))
    }
  )
  
  ### Output con la representación de las asociaciones
  output$cnetplot_KEGGgsea <- renderPlot(
    if (!is.null(gse_KEGG())) {
      DOSE::cnetplot(gse_KEGG())
    }
  )
  
  #### Descarga de la representación de las asociaciones
  output$download_KEGGGSEA_cnet <- downloadHandler(
    filename = "KEGG_GSEA_cnetplot.txt",
    content = function(file) {
      ggsave(file, DOSE::cnetplot(gse_KEGG()))
    }
  )
  
  ### Output con la puntuación de GSEA y asociación de fenotipo
  output$gseaplot_KEGGgsea <- renderPlot(
    if (!is.null(gse_KEGG())) {
      DOSE::gseaplot(gse_KEGG())
    }
  )
  
  #### Descarga de la puntuación de GSEA y asociación de fenotipo
  output$download_KEGGGSEA_fen <- downloadHandler(
    filename = "KEGG_GSEA_gseaplot.txt",
    content = function(file) {
      ggsave(file, DOSE::gseaplot(gse_KEGG()))
    }
  )
  
  # NGS Analysis
  ## Output de la función seq2gene
  s2g <- eventReactive(input$submit7, {
    s2g <- seq2gene(
      seq = NGSpeaks(), tssRegion = c(input$tssRegion_down, input$tssRegion_up),
      flankDistance = input$flankDistance, TxDb = input$TxDb, sameStrand = input$sameStrand
    )
  })
  
  output$downloadGeneFile <- downloadHandler(
    filename = "geneIDs.txt",
    content = function(file) {
      write.table(s2g(), file, sep = "", row.names = FALSE, col.names = FALSE)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

