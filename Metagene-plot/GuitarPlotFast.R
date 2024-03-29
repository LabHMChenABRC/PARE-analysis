#!/usr/bin/env Rscript
# This is a command-line program modified from "https://github.com/LabHMChenABRC/XRNs-DegradomeAnalysis/blob/master/GuitarPlotFast.R"
# Input files (Bed format) should be Single-stranded reads(strand-specific reads coming from the forward strand of RNAs)
# The distribution is calculated and sizes of 5'UTR, CDS, 3'UTR, upstream, and downstream are scaled into 1:2:4:2:1 ratio by Guitar package.
# The program takes 5'Ends of forward alignments on transcripts whose 5'UTR, CDS, 3'UTR, upstream, and downstream sequences all are >= 100 bp
# To reduce RAM requirements and improve performance, some functions of Guitar are modified with vector, matrix, and data.table objects for vectorized computing.
# Those functions are modified based on the source of Guitar v2.10.0.
# written by Bo-Han Hou 2024/03/29

suppressMessages({
  library(argparser, quietly = TRUE)
})
required.packages=c("Guitar","data.table","R.utils","ggplot2","cowplot")

process.reprot <- function(text){
  message(format(Sys.time(), "%Y/%m/%d %H:%M:%S"), " ... ", text)
}
Stop.execution <- function(argv, text = NA) {
  if (!is.na(text) & nchar(text)) {
    warning(text)
  }
  print(argv)
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
check.packages.intalled<-function(required.packages){
  check.res <- lapply( required.packages,function(x) try(find.package(x), silent = TRUE) )
  Not.installed <-unlist(check.res[sapply(check.res, function(x) inherits(x, "try-error"))])
  if( ! is.null(Not.installed)){
    message("The Rscript requires: ", paste(c("Guitar","data.table","R.utils","ggplot2","cowplot"),collapse=", "))
    message(Not.installed)
    stop()
  }
}
check.packages.intalled(required.packages)

###### Main ########
###### Args setting #####
arg.obj   <- arg_parser(
description="
Useage:

Create density file:

  GuitarPlotFast.R -b <bed_folder> -g <gtf_files> -o <output folder>
  
Create metagene plot:

  GuitarPlotFast.R -m -d <metagene plot info> -p <prefix> -o <output folder>
  
  <metagene plot info> is a tabular file providing this info with header
  
  Below example will make two groups of  WT/MutA, WT/MutB in a metagene plot
  
  file polygroup
  
  <PathofDensity.WT> <1,2> 
  
  <PathofDensity.MutA> <1> 
  
  <PathofDensity.MutB> <2> ", hide.opts = TRUE)
arg.obj   <- add_argument(arg.obj, arg = "--bed_folder"    ,type = "character", help = "folder contains bed.gz files")
arg.obj   <- add_argument(arg.obj, arg = "--gtf"           ,type = "character", help = "gene model annotation file in gtf format")
arg.obj   <- add_argument(arg.obj, arg = "--metaplot"      ,flag = TRUE, help = "switch to create metagene plot")
arg.obj   <- add_argument(arg.obj, arg = "--density_info"  ,type = "character", help = "tablur file contains <File: Denisty File Location > and <polygroup: Group the density> per row")
arg.obj   <- add_argument(arg.obj, arg = "--prefix"        ,type = "character", help = "a prefix of metagene plot file")
arg.obj   <- add_argument(arg.obj, arg = "--output"        ,type = "character", help = "ouput folder")
argv      <- parse_args(arg.obj)
if ( argv$metaplot == TRUE){
  #### make metagene Plot ####
  requiredArgs=c("desity_folder","density_info","prefix", "output")
  #### Check args ####
  argv.na.check <- unlist(sapply(requiredArgs, function(x) is.na(argv[[x]]) | argv[[x]] == FALSE))
  if ( any(argv.na.check) ) {
    Stop.execution(arg.obj)
  }
  #### make metagene Plot functions####
  suppressWarnings({
    suppressMessages({
      library(data.table)
      library(ggplot2)
      library(cowplot)
    })
  })
  load_mRNA_density<-function(density_info_file){
    density_info<-fread(file = density_info_file,colClasses = c("character","character"))
    if ( ! all(file.exists(density_info$file)) ){
      stop("Some density files are not exists")
    }
    Grouplevels <-sub(".mrna.density","",basename(unique(density_info$file)))
    mRNA.density.files<-unique(density_info$file)
    names(mRNA.density.files) <- basename(mRNA.density.files)
    if( any(grep(",",density_info$plotgroup)) ){
      density_info<-rbind(
        density_info[grep(',',plotgroup), .(file,plotgroup= strsplit(plotgroup,split = ',')), by = .I][,.(file,plotgroup)][,lapply(.SD, unlist),file],
        density_info[!grep(',',plotgroup)]
      )
    }
    mRNA.density <-rbindlist(lapply(mRNA.density.files,fread),idcol = "ID"
    )[,group:=factor(group,levels = Grouplevels)
    ][density_info[,file:=basename(file)],on="ID==file",allow.cartesian=TRUE][]
    return(mRNA.density)
  }
  
  get.Model <- function(Feature.Width,height,label=c("1kb","5'UTR","CDS","3'UTR","1kb")){
    if(length(Feature.Width)!=5){
      stop("Feature.Width should have 5 values!")
    }
    if( !is.null(label) && length(label)!=5){
      stop("label should have 5 values!")
    }
    Feature.Width     <- Feature.Width/sum(Feature.Width)
    model.height      <- height*1.06
    model.height.diff <- model.height-height
    Feature.model <- data.table(
      comp  = c("promoter","5'UTR","CDS","3'UTR","tail"),
      label = label,
      width = Feature.Width,
      xmin  = c(0,head(cumsum(Feature.Width),-1)) + 0.001,
      xmax  = cumsum(Feature.Width),
      ymin  = model.height,
      ymax  = model.height,
      label_pos = model.height
    )[,`:=`(mid = (xmax+xmin)/2,
            label_pos = label_pos + 12/10 * model.height.diff)
    ][grep("UTR",comp),`:=`(ymin=ymin-3/10*model.height.diff,ymax=ymax+3/10*model.height.diff)
    ][grep("CDS",comp),`:=`(ymin=ymin-4/10*model.height.diff,ymax=ymax+4/10*model.height.diff)
    ][width!=0
    ][,bottom:=min(ymin)*0.975][]
    return(Feature.model)
  }
  
  Density.plot <- function(Density.dt,ncol=NULL,colorset=NULL){
    if(is.null(colorset)){
      if(is.factor(Density.dt$group)){
        colorset <- setNames(getColors(length(levels(Density.dt$group))),levels(Density.dt$group))
      }else{
        colorset <- setNames(getColors(length(unique(Density.dt$group))),unique(Density.dt$group))
      }
    }
    ggplot()+
      geom_area(data = Density.dt,alpha=0.5,position = "identity",aes(x=x,y=density,fill=group,color=group)) +
      scale_fill_manual(name = "group",values = colorset) +
      scale_colour_manual(name = "group",values = colorset) +
      scale_y_continuous(expand = expansion(mult=c(0,0.05)))+
      labs(x=NULL)+
      facet_wrap(~ID,ncol = ncol)+
      theme_half_open()+
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            legend.title = element_blank(),
            legend.position = 'top')
  }
  add.model<-function(componentStructure,fontsize=10){
    list(geom_segment(data = componentStructure[grep("promoter|tail",comp),],aes(x=xmin,xend=xmax,y=ymin,yend=ymax),color='black'),
         geom_rect(data = componentStructure[grep("UTR",comp),],aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='grey70'),
         geom_rect(data = componentStructure[grep("CDS",comp),],aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='grey50'),
         geom_text(data = componentStructure,aes(x=mid,y=label_pos,label=label),size=fontsize/ggplot2::.pt),
         geom_segment(data = head(componentStructure,-1),aes(x=xmax,xend=xmax,y=0,yend=bottom),linetype='dashed',color="grey")
    )
  }
  #### Main: make metagene Plot ####
  # Fetch and merge *.density files to make plot.
  mRNA.density <- load_mRNA_density(density_info_file = argv$density_info)
  model.position=max(mRNA.density$density) * 1.05
  if ( !dir.exists(argv$output) ) {
    dir.create(argv$output,recursive = TRUE,showWarnings = FALSE)
  }
  setwd(dir = argv$output)
  mRNA.p <- ggplot()+
    geom_area(data = mRNA.density,alpha=0.3,position = "identity",aes(x=x,y=density,fill=group,color=group)) +
    add.model(componentStructure = get.Model(c(1,2,4,2,1),
                                             height = model.position,
                                             label = c("","5'UTR","CDS","3'UTR","")),fontsize = 8)+
    scale_y_continuous(expand = expansion(mult=c(0,0.05)))+
    guides(color = guide_legend(nrow = 1),fill = guide_legend(nrow = 1))+
    labs(x=NULL,y="density")+
    facet_grid(~plotgroup,scales = "free_y")+
    cowplot::theme_half_open(font_size = 8)+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          legend.position = 'top')
  figure_file=sprintf("%s.mRNA-metaplot.pdf",argv$prefix)
  ggsave(filename = figure_file,plot = mRNA.p,width = 0.2+7* length(unique(mRNA.density$plotgroup)), height = 6,units = "cm")
  message("Output figure: ",normalizePath(figure_file))
}else{
  #### make density ####
  #### Check args ####
  requiredArgs=c("bed_folder", "gtf", "output")
  argv.na.check <- unlist(sapply(requiredArgs, function(x) is.na(argv[[x]]) | argv[[x]] == FALSE))
  if ( any(argv.na.check) ) {
    Stop.execution(arg.obj)
  }
  if ( ! dir.exists(argv$bed_folder) ) {
    Stop.execution(arg.obj)
  }
  if ( ! file.exists(argv$gtf) ) {
    Stop.execution(arg.obj)
  }
  if ( !dir.exists(argv$output) ) {
    dir.create(argv$output,recursive = TRUE,showWarnings = FALSE)
  }
  if ( !dir.exists(argv$output) ) {
    message("output folder can't be created.")
    Stop.execution(arg.obj)
  }
  #### make density functions ####
  message("Load packages")
  suppressWarnings({
    suppressMessages({
      library(data.table)
      library(Guitar)
    })
  })
  
  GuitarPlotFast <- function(txGTF = NULL, txGFF = NULL, txGenomeVer = NULL, txTxdb = NULL,
                             alignment_direction=c("forward","reverse","ignore"),
                             txGuitarTxdb = NULL, txGuitarTxdbSaveFile = NA, stBedFiles = NULL, 
                             stGRangeLists = NULL, stGroupName = NULL, stAmblguity = 5, 
                             stSampleNum = 3, stSampleModle = c("Equidistance","random"), txfiveutrMinLength = 100, 
                             txcdsMinLength = 100, txthreeutrMinLength = 100, txlongNcrnaMinLength = 100, 
                             txlncrnaOverlapmrna = FALSE, txpromoterLength = 1000, txtailLength = 1000, 
                             txAmblguity = 5, txPrimaryOnly = FALSE, txTxComponentProp = NULL, 
                             txMrnaComponentProp = NULL, txLncrnaComponentProp = NULL, 
                             mapFilterTranscript = TRUE, headOrtail = TRUE, enableCI = TRUE, 
                             pltTxType = c("tx", "mrna", "ncrna"), overlapIndex = 1, 
                             siteLengthIndex = 1, adjust = 1, CI_ResamplingTime = 1000, 
                             CI_interval = c(0.025, 0.975), miscOutFilePrefix = NA) 
  {
    stSampleModle<-match.arg(stSampleModle)
    # stSampleModle<-"Equidistance"
    pltTxType    <-match.arg(pltTxType,several.ok=TRUE)
    
    genomeVersion2Txdb <- list(hg18 = "TxDb.Hsapiens.UCSC.hg18.knownGene", 
                               hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene", hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene", 
                               mm9 = "TxDb.Mmusculus.UCSC.mm9.knownGene", mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene")
    if (headOrtail) {
      txpromoterLength <- txpromoterLength
      txtailLength <- txtailLength
    }  else {
      txpromoterLength <- 0
      txtailLength <- 0
    }
    
    process.reprot("Start to make density file")
    
    guitarTxdb <- getGuitarTxdb.fast(txGTF = txGTF, txGFF = txGFF, 
                                     txGenomeVer = txGenomeVer, txTxdb = txTxdb, txGuitarTxdb = txGuitarTxdb, 
                                     txfiveutrMinLength = txfiveutrMinLength, txcdsMinLength = txcdsMinLength, 
                                     txthreeutrMinLength = txthreeutrMinLength, txlongNcrnaMinLength = txlongNcrnaMinLength, 
                                     txlncrnaOverlapmrna = txlncrnaOverlapmrna, txpromoterLength = txpromoterLength, 
                                     txtailLength = txtailLength, txAmblguity = txAmblguity, 
                                     txTxComponentProp = txTxComponentProp, txMrnaComponentProp = txMrnaComponentProp, 
                                     txLncrnaComponentProp = txLncrnaComponentProp, txPrimaryOnly = txPrimaryOnly, 
                                     pltTxType = pltTxType, genomeVersion2Txdb)
    
    
    if (!(is.na(txGuitarTxdbSaveFile))) {
      txGuitarTxdbSaveFile <- paste(txGuitarTxdbSaveFile,"GuitarTxdb", sep = "-")
      saveRDS(guitarTxdb, file = txGuitarTxdbSaveFile)
    }
    
    if (any(!pltTxType %in% guitarTxdb$txTypes)) {
      pltTxType.notfound <- pltTxType[!pltTxType %in% guitarTxdb$txTypes]
      process.reprot(paste(pltTxType.notfound,"is not found in tx annotation"))
      stop(paste("Cannot plot distribution for", pltTxType.notfound))
    }
    
    process.reprot("import BED file")
    sitesGroup <- Guitar:::.getStGroup(stBedFiles = stBedFiles, 
                                       stGRangeLists = stGRangeLists, 
                                       stGroupName = stGroupName)
    GroupNames <- names(sitesGroup)
    sitesGroupNum <- length(sitesGroup)
    sitesPointsNormlize <- list()
    sitesPointsRelative <- list()
    pointWeight         <- list()
    # if is 
    if ( alignment_direction=="reverse" ){
      for (i in seq_len(sitesGroupNum)) {
        sitesGroup[[i]]<-invertStrand(sitesGroup[[i]])
      }
    }
    ignore.strand<-ifelse(alignment_direction=="none",TRUE,FALSE)
    
    for (i in seq_len(sitesGroupNum)) {
      GroupName = GroupNames[[i]]
      process.reprot(sprintf("Sample points (%s) for %s", stSampleNum, GroupName))
      
      sitesPoints <- samplePoints.fast(sitesGrangelists = sitesGroup[i], stSampleNum = stSampleNum, 
                                       stAmblguity = stAmblguity, pltTxType = pltTxType, 
                                       stSampleModle = stSampleModle, mapFilterTranscript = mapFilterTranscript,
                                       ignore.strand=ignore.strand,
                                       guitarTxdb = guitarTxdb)
      
      for (txType in pltTxType) {
        process.reprot(paste("normalize for",txType,sep = " "))
        sitesPointsNormlize[[txType]][[GroupName]] <- normalize.fast(sitesPoints, guitarTxdb, txType, overlapIndex, siteLengthIndex)
        sitesPointsRelative[[txType]][[GroupName]] <- sitesPointsNormlize[[txType]][[GroupName]][[1]]
        pointWeight[[txType]][[GroupName]] <- sitesPointsNormlize[[txType]][[GroupName]][[2]]
      }
    }
    p.list=list()
    
    for (txType in pltTxType) {
      process.reprot(sprintf("Calucate %s density", txType) )
      title <- switch(txType,
                      tx="Distribution on Transcript",
                      mrna="Distribution on mRNA", 
                      ncrna="Distribution on ncRNA")
      densityDataframe_CI <- Guitar:::.generateDensity_CI(sitesPointsRelative[[txType]],pointWeight[[txType]], CI_ResamplingTime, adjust = adjust, enableCI = enableCI)
      p.list[[txType]] <- Guitar:::.plotDensity_CI(densityDataframe_CI, componentWidth = guitarTxdb[[txType]]$componentWidthAverage_pct, headOrtail, title, enableCI = enableCI)
    }
    process.reprot("End")
    return(p.list)
  } 
  getGuitarTxdb.fast<-function (txGTF = NULL, txGFF = NULL, 
                                txGenomeVer = NULL, txTxdb = NULL, 
                                txGuitarTxdb = NULL, txfiveutrMinLength = 100, 
                                txcdsMinLength = 100, txthreeutrMinLength = 100, 
                                txlongNcrnaMinLength = 100, txlncrnaOverlapmrna = FALSE, 
                                txpromoterLength = 1000, txtailLength = 1000, 
                                txAmblguity = txAmblguity, txTxComponentProp = NULL, 
                                txMrnaComponentProp = NULL, txLncrnaComponentProp = NULL, 
                                txPrimaryOnly = FALSE, pltTxType = c("tx", "mrna", "ncrna"), genomeVersion2Txdb) 
  {
    if (!(is.null(txGuitarTxdb))) {
      if (!file.exists(txGuitarTxdb)) {
        stop(paste0("\"txGuitarTxdb\" is not exist, please double check!"))
      }
      process.reprot(sprintf("Using existed guitarTxdb file: %s", basename(txGuitarTxdb)))
      guitarTxdb <- readRDS(txGuitarTxdb)
      if(!is.null(txTxComponentProp)){
        guitarTxdb$tx$componentWidthAverage_pct      <- setNames(txTxComponentProp/sum(txTxComponentProp),c("promoter","rna","tail"))
        guitarTxdb$tx$componentStartAverage_pct      <- setNames(c(0,head(cumsum(guitarTxdb$tx$componentWidthAverage_pct),-1)),c("promoter","rna","tail"))
      }
      if(!is.null(txLncrnaComponentProp)){
        guitarTxdb$ncrna$componentWidthAverage_pct   <- setNames(txLncrnaComponentProp/sum(txLncrnaComponentProp),c("promoter","ncrna","tail"))
        guitarTxdb$ncrna$componentStartAverage_pct   <- setNames(c(0,head(cumsum(guitarTxdb$ncrna$componentWidthAverage_pct),-1)),c("promoter","ncrna","tail"))
      }
      if(!is.null(txMrnaComponentProp)){
        guitarTxdb$mrna$componentWidthAverage_pct    <- setNames(txMrnaComponentProp/sum(txMrnaComponentProp),c("promoter","utr5","cds","utr3","tail"))
        guitarTxdb$mrna$componentStartAverage_pct    <- setNames(c(0,head(cumsum(guitarTxdb$mrna$componentWidthAverage_pct),-1)),c("promoter","utr5","cds","utr3","tail"))
      }
    }else{
      process.reprot("make GuitarTxdb")
      txdb <- Guitar:::.getTxdb(txGTF = txGTF, txGFF = txGFF, txGenomeVer = txGenomeVer, 
                                txTxdb = txTxdb, genomeVersion2Txdb)
      guitarTxdb <- makeGuitarTxdb(txdb, txfiveutrMinLength = txfiveutrMinLength, 
                                   txcdsMinLength = txcdsMinLength, txthreeutrMinLength = txthreeutrMinLength, 
                                   txlongNcrnaMinLength = txlongNcrnaMinLength, txlncrnaOverlapmrna = txlncrnaOverlapmrna, 
                                   txpromoterLength = txpromoterLength, txtailLength = txtailLength, 
                                   txAmblguity = txAmblguity, txTxComponentProp = txTxComponentProp, 
                                   txMrnaComponentProp = txMrnaComponentProp, txLncrnaComponentProp = txLncrnaComponentProp, 
                                   txPrimaryOnly = txPrimaryOnly, pltTxType = pltTxType)
    }
    return(guitarTxdb)
  }
  
  samplePoints.fast <- function(sitesGrangelists, stSampleNum = 5, stAmblguity = 5, 
                                pltTxType = c("tx", "mrna", "ncrna"), 
                                stSampleModle = "Equidistance",
                                ignore.strand = FALSE,
                                mapFilterTranscript = FALSE, guitarTxdb) 
  {
    mkSitesPoints <- function(sitesWidth,stSampleModle,stSampleNum){
      Sample.fun <- switch(stSampleModle,
                           Equidistance=function(x, i) {
                             if (i == 1) {
                               round(x/2)
                             }    else {
                               round(seq(1, x - 1, length.out = i))
                             }
                           },
                           random=function(x, i) {
                             if (i == 1) {
                               round(x/2)
                             }
                             else {
                               sort(sample(x, i, replace = FALSE))
                             }
                           })
      # build sitesPoints matrix based on unique sitesWidth
      sitesPointsVector <- vector(mode = "numeric",length(sitesWidth))
      sitesWidth.uni    <- sort(unique(sitesWidth))
      sitesWidth.idx    <- setNames(1:length(sitesWidth.uni), sitesWidth.uni)
      sitesPoints.uni   <- t(vapply(sitesWidth.uni,Sample.fun,i=stSampleNum,FUN.VALUE = numeric(stSampleNum)))
      if( all(dim(sitesPoints.uni) == rev(c(length(sitesWidth.idx),stSampleNum))) ){
        sitesPoints.uni<-t(sitesPoints.uni)
      }
      
      # Covert sitesWidth to sitesPointsVector
      sitesPointsVector <- as.numeric(t(sitesPoints.uni[sitesWidth.idx[as.character(sitesWidth)],]))
      
      return(sitesPointsVector)
    }
    
    sitesGRangeDataframe <- list()
    stSampleNum <- 2 * stSampleNum - 1
    
    for (txType in pltTxType) {
      
      process.reprot(sprintf("Step1: Map to Transcripts (%s)",txType))
      mapsiteGRanges <- GRangesListmapToTranscripts.fast(site = sitesGrangelists[[1]], 
                                                         mapFilterTranscript = mapFilterTranscript,   
                                                         transcripts = guitarTxdb[[txType]]$tx,
                                                         ignore.strand=ignore.strand)
      sitesWidth <- width(mapsiteGRanges)
      process.reprot(paste("Step2: Build sitesPoints by",stSampleModle,sep=" "))
      sitesPointsVector <- mkSitesPoints(sitesWidth = sitesWidth,stSampleModle = stSampleModle,stSampleNum = stSampleNum)
      
      process.reprot("Step3: Make sitesPoints in GRanges")
      sitesGRangeDataframe[[txType]]         <- mapsiteGRanges[rep(seq_len(length(mapsiteGRanges)), each=stSampleNum),]
      start(sitesGRangeDataframe[[txType]])  <- sitesPointsVector+start(sitesGRangeDataframe[[txType]]) - 1L
      end(sitesGRangeDataframe[[txType]])    <- start(sitesGRangeDataframe[[txType]]) + 1L
      strand(sitesGRangeDataframe[[txType]]) <- "*"
      
      # use table() to count and sign them
      xHits.times     <- table(mapsiteGRanges$xHits)
      df <- data.frame(sitesLength     = sitesWidth, 
                       xHits           = mapsiteGRanges$xHits, 
                       pointsOverlapTx = as.vector(xHits.times[as.character(mapsiteGRanges$xHits)]))
      mcols(sitesGRangeDataframe[[txType]]) <- df[rep(seq_len(nrow(df)), each=stSampleNum), ]
      
    }
    return(sitesGRangeDataframe)
  }
  

  
  collapseGRanges <- function(x){
    return(
      as.data.table(x)[,.(score=as.integer(.N),qnames=list(.I)),by=.(seqnames,start,end,strand)
      ][,GRanges(.SD)]
    )
  }
  
  uncollapseGRanges<-function(gr.un){
    gr<-gr.un[rep(seq_along(gr.un$score),gr.un$score),]
    mcols(gr)$score <-NULL
    mcols(gr)$qnames <-NULL
    gr$xHits<-unlist(gr.un$qnames)
    return(gr)
  }
  GRangesListmapToTranscripts.fast <- function(site, mapFilterTranscript = FALSE,transcripts,ignore.strand=FALSE)
  { 
    
    names(site) <- 1:length(site)
    if(mapFilterTranscript) {
      xWidthes <- sum(width(site))
      names(xWidthes) <- names(site)
    }
    if(is(site,"CompressedGRangesList")){  
      site <- unlist(site)
    }
    site <- collapseGRanges(site)
    # Give new levels of seqnames which are exist in transcripts, because some ref sequence did not annotated any gene features
    site <- site[seqnames(site) %in% levels(seqnames(unlist(transcripts)))]
    seqlevels(site) <- as.character(levels(seqnames(unlist(transcripts))))
    
    tx_coord <- mapToTranscripts(site, transcripts,ignore.strand=ignore.strand)
    hit_txHit.DF <-as.data.table(mcols(tx_coord))[,.(.GRP),by=c("xHits","transcriptsHits")
    ][,.(xHits,txHits=transcriptsHits,score=site[xHits]$score,qnames=site[xHits]$qnames,GRP)
    ][,DataFrame(.SD[,],row.names=.SD$GRP)][,c("xHits","txHits","score","qnames")]
    
    tx_coord_grouped <- split(tx_coord, rownames(hit_txHit.DF))
    mapping_reduced  <- reduce(tx_coord_grouped)
    
    # some sites map to tx has multiple regions after reduce because of isoform of tx
    tx.region.num            <- table(names(unlist(width(mapping_reduced))))
    tx.region.selected       <- names(tx.region.num[tx.region.num == 1])
    tx_coord_filtered        <- unlist(mapping_reduced[tx.region.selected])
    mcols(tx_coord_filtered) <- hit_txHit.DF[names(tx_coord_filtered),]
    
    
    # remove hits whose length smaller than sites because of isoform of tx
    if(mapFilterTranscript) {
      tx_coord_filtered_width <- width(tx_coord_filtered)
      tx_coord_filtered <- tx_coord_filtered[tx_coord_filtered_width == xWidthes[tx_coord_filtered$xHits]]
    }
    idx <- GenomicRanges:::get_out_of_bound_index(tx_coord_filtered)
    if (length(idx) != 0) {
      tx_coord_filtered <- tx_coord_filtered[-idx]
    }
    
    return(sort(uncollapseGRanges(tx_coord_filtered)))
  }
  
  normalize.fast <- function(sitesGRanges, guitarTxdb, txType, overlapIndex, siteLengthIndex) 
  {
    sitesPointsPositionNormalize <- list()
    sitesPointsPosTx        <- list()
    sitesPointsPosTx        <- end(sitesGRanges[[txType]])
    names(sitesPointsPosTx) <- seqnames(sitesGRanges[[txType]])
    startPointMat           <- guitarTxdb[[txType]]$startPoint[names(sitesPointsPosTx), ]
    startPointDiffer        <- startPointMat - sitesPointsPosTx
    
    max_which.fast<-function(x){
      x[x>=0] <- -Inf # assign -Inf to prevent from max.col() catch them
      setNames(max.col(x),rownames(x))
    }
    sitesPointsComponet <- max_which.fast(startPointDiffer)
    sitesPointsPositionComponet <- startPointDiffer[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)] * -1
    
    # normalize
    sitesPointsComponetWidthAvg  <- guitarTxdb[[txType]]$componentWidthAverage_pct[sitesPointsComponet]
    sitesPointsComponetStart_pct <- guitarTxdb[[txType]]$componentStartAverage_pct[sitesPointsComponet]
    componentWidthMat            <- guitarTxdb[[txType]]$componentWidth[names(sitesPointsPosTx),]
    sitesPointsComponetWidth     <- componentWidthMat[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)]
    sitesPointsPositionNormalize <- sitesPointsPositionComponet/sitesPointsComponetWidth * sitesPointsComponetWidthAvg + sitesPointsComponetStart_pct
    names(sitesPointsPositionNormalize) <- sitesGRanges[[txType]]$xHits
    
    
    # weight
    sitesComponet_pct       <- guitarTxdb[[txType]]$componentWidthPtc[names(sitesPointsComponet), ]
    col.levels              <- sort(unique(sitesPointsComponet))
    sitesPointsComponet_pct <- setNames(numeric(length(sitesPointsComponet)),colnames(sitesComponet_pct)[sitesPointsComponet])
    for (i in col.levels){
      sitesPointsComponet_pct[sitesPointsComponet == i] <- sitesComponet_pct[sitesPointsComponet == i,i]
    }
    
    sitesPointsWeight        <- sitesPointsComponetWidthAvg/(sitesGRanges[[txType]]$pointsOverlapTx^overlapIndex)/sitesPointsComponet_pct * (sitesGRanges[[txType]]$sitesLength^siteLengthIndex)
    names(sitesPointsWeight) <- sitesGRanges[[txType]]$xHits
    return(list(sitesPointsPositionNormalize, sitesPointsWeight))
  }
  
  bedAlignemntTo5Pend <- function(file,outdir,MAPQ=NULL){
    # input file is produced via bedtools bam2bed
    # Unique hit with 255 MAPQ if STAR aligner is used to map.
    process.reprot("start to covert alignments to 5'Ends")
    if (requireNamespace('data.table')){
      for (file.i in unlist(file)){
        message("input:",file.i) 
        NewFile <- sprintf("%s/%s",outdir,basename(sub(".bed.gz$",".5P.bed.gz",file.i)))
        if (!file.exists(NewFile) ){
          Ailgnments <- data.table::fread(file.i,sep="\t",header = FALSE)
          if( !is.null(MAPQ)) {
            Ailgnments<-Ailgnments[V5>=MAPQ] # filter by MAPQ
          }
          Ailgnments[V6=="+",V3:=V2+1]
          Ailgnments[V6=="-",V2:=V3-1]
          data.table::fwrite(Ailgnments,file=NewFile,sep="\t",col.names=F)
          message("output:",NewFile) 
        }else{
          message(paste("Skip!","output file",NewFile, "is already existed.",sep=" ")) 
        }
      }
    }
  }
  txNum.guitarTxdb<-function(guitarTxdb,pltTxTypes){
    if ( 'tx' %in%  pltTxTypes){
      message( sprintf("\t%s selected transcripts(tx) have minimum size(bp) of features:",
                       length(D[['tx']]$tx)))
      message(paste0(paste("\t",
                           capture.output(setNames(c(txpromoterLength,txtailLength),
                                                   c("Upstream","Downstream"))
                           )
      ), collapse = "\n"))
    }
    if ( 'mrna' %in%  pltTxTypes){
      
      message( sprintf("\t%s selected transcripts(mrna) have minimum size(bp) of features:",
                       length(guitarTxdb[['mrna']]$tx)))
      message(paste0(paste("\t",
                           capture.output(setNames(c(txpromoterLength,txfiveutrMinLength,txcdsMinLength,txthreeutrMinLength,txtailLength),
                                                   c("Upstream","5'UTR","CDS","3'UTR","Downstream"))
                           )
      ), collapse = "\n"))
    }
    if ( 'ncrna' %in%  pltTxTypes){
      
      message( sprintf("\t%s selected transcripts(ncrna) have minimum size(bp) of features:",
                       length(guitarTxdb[['ncrna']]$tx)))
      message(paste0(paste("\t",
                           capture.output(setNames(c(txpromoterLength,txlongNcrnaMinLength,txtailLength),
                                                   c("Upstream","non-coding RNA","Downstream"))
                           )
      ), collapse = "\n"))
    }
  }
  
  
  #### Main: make density file####
  bed_folder <-normalizePath(argv$bed_folder)
  gtf_file   <-normalizePath(argv$gtf)
  output     <-normalizePath(argv$output)
  setwd(dir = output)
  
  # Build a GuitarTxdb, and txMrnaComponentProp/txTxComponentProp/txLncrnaComponentProp should be given at this time. 
  # GuitarTxdb object can be saved by saveRDS()
  # txdb is created by makeTxDbFromXXXX() function of GenomicFeatures package.
  txGuitarTxdb.file <-sub("gtf$","guitarTxdb",gtf_file)
  pltTxType = c("mrna") #"tx","mrna","ncrna"
  txMrnaComponentProp <- c(1,2,4,2,1)
  if ( file.exists(txGuitarTxdb.file) ){
    message("use exists GuitarTxdb file:",txGuitarTxdb.file)
    guitarTxdb <- readRDS(txGuitarTxdb.file)
  }else{
    message("make GuitarTxdb file:")
    txdb <- makeTxDbFromGFF(file = gtf_file,format = "gtf")
    guitarTxdb <- makeGuitarTxdb(txdb,
                                 txpromoterLength = 100,
                                 txcdsMinLength = 100,
                                 txtailLength = 100,
                                 txTxComponentProp =NULL,
                                 txLncrnaComponentProp = NULL,
                                 txMrnaComponentProp = txMrnaComponentProp)
    saveRDS(guitarTxdb,txGuitarTxdb.file)
    rm(txdb)
    
    
    for (pltTxType.i in pltTxType){
      txNum.guitarTxdb(guitarTxdb,pltTxType.i)
    }
    rm(guitarTxdb)
  }
  
  # Make 5'Ends of primary alignments file (5P.bed.gz)
  Bed5PEnd_dir=sprintf("%s/%s",getwd(),"5PEnds")
  dir.create(path = Bed5PEnd_dir,recursive = TRUE,showWarnings = FALSE)
  AlignemntFiles  <- list.files(bed_folder,pattern = "*.bed.gz$",full.names = TRUE)
  bedAlignemntTo5Pend(file=AlignemntFiles,outdir = Bed5PEnd_dir,MAPQ = 255) # Primary alignment has 255 MAPQ in STAR aligner
  bed.Files        <- list.files(Bed5PEnd_dir,pattern = "*.5P.bed.gz$",full.names = TRUE)
  names(bed.Files) <- basename(sub(".5P.bed.gz","",bed.Files))
  
  # Start to calculate density:
  # Less than 15 min for mRNA density calculation for stranded single ends 20M library.
  density.dir=sprintf("%s/%s",getwd(),"density")
  dir.create(density.dir,showWarnings = FALSE,recursive = TRUE)
  run.time=list()
  for ( ID in names(bed.Files) ){
    run.time[[ID]]<-system.time({
      suppressWarnings({
      count <- GuitarPlotFast(txGuitarTxdb = txGuitarTxdb.file,
                              stBedFiles = as.list(bed.Files[ID]),
                              headOrtail = TRUE,
                              enableCI = FALSE,
                              txpromoterLength = 100,
                              txcdsMinLength = 100,
                              txtailLength = 100,
                              stSampleNum = 1,
                              mapFilterTranscript = TRUE,
                              pltTxType = pltTxType,
                              txMrnaComponentProp = txMrnaComponentProp,
                              alignment_direction = "forward",
                              stGroupName = ID)
      })
    })
    message(sprintf("elapsed time: %.0f sec",run.time[[ID]][[3]]))
    for (pltTxType.i in pltTxType){
      density_file=sprintf("%s/%s.%s.density",density.dir,ID,pltTxType.i)
      fwrite(count[[pltTxType.i]]$data,file=density_file,sep="\t")
      message( sprintf("Output density file: %s",density_file) )
    }
  }
  
  time.log <- t(data.table::as.data.table(run.time))
  colnames(time.log)<-names(run.time[[1]])
  message(sprintf("mean of elapsed time: %.2f min",mean(time.log[,"elapsed"])/60))
}
