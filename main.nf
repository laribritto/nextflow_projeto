params.results = "resultados"
params.cogdata = "$projectDir/cogdata.rda"
params.sspids = "$projectDir/sspids.rda"
params.cogids = "$projectDir/cogids.rda"
params.phyloTree = "$projectDir/phyloTree.rda"

process METRICAS {
    publishDir "$params.results/resultados", mode: 'copy'
    container 'laribritto/geneplast:v1.1'

    input:
    path cogdata
    path sspids
    path cogids


    output:
    path 'grafico.pdf'

    script:
    """
    #!/usr/bin/Rscript

    # Importando dados
    library(geneplast)
    library(tidyr)
    library(ggplot2)
    load("$cogdata")
    load("$sspids")
    load("$cogids")
   

    ogp <- gplast.preprocess(cogdata = cogdata, sspids = sspids, cogids = cogids, verbose = TRUE)
    ogp <- gplast(ogp, verbose = FALSE)
    res <- gplast.get(ogp, what = "results")
    head(res)

    subset <- head(res)
    subset\$cog_id <- row.names(subset)

    subset <- pivot_longer(subset, -cog_id)

    pdf(file = "grafico.pdf")
    ggplot(subset) +
        geom_bar(aes(x = cog_id, y = value, fill = name), stat = "identity", position = 'dodge') +
        scale_fill_grey() +
        theme_bw()
    dev.off()
    """
}

process RAIZ {
    publishDir "$params.results/resultados", mode: 'copy'
    container 'laribritto/geneplast:v1.1'

    input:
    path cogdata
    path phyloTree
    path sspids

    output:
    path 'gproot_NOG40170_9606LCAs.pdf'
    path 'gproot_9606LCAs.pdf'

    script:
    """
    #!/usr/bin/Rscript

    # Análise evolutiva
    library(geneplast)
    library(tidyr)
    library(ggplot2)
    load("$cogdata")
    load("$phyloTree")
    load("$sspids")

    ogr <- groot.preprocess(cogdata = cogdata, phyloTree = phyloTree, spid = "9606", verbose = FALSE)

    set.seed(1)
    ogr <- groot(ogr, nPermutations = 100, verbose = FALSE)

    res <- groot.get(ogr, what = "results")
    head(res)

    groot.plot(ogr, whichOG = "NOG40170")

    groot.plot(ogr, plot.lcas = TRUE)
    """
}

workflow {
    METRICAS(params.cogdata,params.sspids,params.cogids)
    RAIZ(params.cogdata,params.sspids,params.phyloTree)
}

workflow.onComplete {
    log.info(workflow.success ? "\nVocê é uma máquina de vencer! :)" : "\nO fracasso é inevitável! :(")
}
