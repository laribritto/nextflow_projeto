params.results = 'resultados'
params.cogdata_file = '$projectDir/cogs (1).rda'
params.cogs_of_interest_file = '$projectDir/cogs_of_interest (1).csv'

process METRICAS {
    publishDir "$params.results/resultados", mode: 'copy'
    container 'laribritto/geneplast:v1.1'

    input: path cogs_of_interest
           patg cogs

    output: path 'grafico.pdf'
            path 'cogs.rda'
            path 'phyloTree.rda'
            path 'sspids.rda'
    script:
    """
    #!/usr/bin/Rscript

    # Importando dados
    library(geneplast)
    library(tidyr)
    library(ggplot2)

    ogp <- gplast.preprocess(cogdata=cogs, sspids=sspids, cogids=cogs_of_interest, verbose=TRUE)
    ogp <- gplast(ogp, verbose=FALSE)
    res <- gplast.get(ogp, what="results")
    head(res)

    subset <- head(res)
    subset\$cog_id <- row.names(subset)

    subset <- pivot_longer(subset,-cog_id)

   pdf(file="grafico.pdf")
   ggplot(subset) + 
    geom_bar(aes(x = cog_id, y = value, fill = name), stat = "identity", position = 'dodge') + 
    scale_fill_grey() +
    theme_bw()
   dev.off()
   save(cogs, file = "cogs.rda")
   save(phyloTree, file = "phyloTree.rda")
   save(sspids, file = "sspids.rda")
    """
}

process RAIZ {
    publishDir "$params.results/resultados", mode: 'copy'
    container 'laribritto/geneplast:v1.1'

    input: path cogs
           path phyloTree
           path sspids
        
    output: path 'gproot_COG0085_9606LCAs.pdf'
            path 'gproot_9606LCAs.pdf'

    script:
    """
    #!/usr/bin/Rscript

    # Análise evolutiva
    library(geneplast)
    library(tidyr)
    library(ggplot2)
    load('cogs.rda')
    load('phyloTree.rda')
    load('sspids.rda')

    ogr <- groot.preprocess(cogdata=cogs, phyloTree=phyloTree, spid="9606", verbose=FALSE)

    set.seed(1)
    ogr <- groot(ogr, nPermutations=100, verbose=FALSE)

    res <- groot.get(ogr, what="results")
    head(res)

    groot.plot(ogr, whichOG="COG0085")

    groot.plot(ogr, plot.lcas = TRUE)
    """
}

workflow {
    METRICAS(params.cogdata_file,params.cogs_of_interest_file) 
    RAIZ(METRICAS.out[1],METRICAS.out[2],METRICAS.out[3]) 
}

workflow.onComplete {
    log.info ( workflow.success ? "\nVocê é uma máquina de vencer! :)" : "\nO fracasso é inevitável! :(" )
}
