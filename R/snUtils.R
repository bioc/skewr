#################################################
##
## Functions for getting MethyLumiSet objects
## and for preprocessing of them
##
#################################################

getBarcodes <- function(path = getwd(), recurse = FALSE){
  unique <- vector(mode = 'character')
  for(i in 1:length(path)){
    unique <- c(unique, grep('_Red.idat$', list.files(path[i], recursive = recurse), value = TRUE))
    unique <- gsub('_Red.idat$', '', unique)
  }

  unique
}


getMethyLumiSet <- function(path = getwd(), barcodes = NULL,
                            norm = c("none", "illumina", "SWAN", "dasen"), bg.corr = TRUE){
  norm <- match.arg(norm)
  if(is.null(barcodes)) barcodes <- getBarcodes(path)
  object <- methylumIDAT(barcodes = barcodes, idatPath = path)
  if(norm != "none") object <- preprocess(object, norm = norm, bg.corr = bg.corr)
  object
}

preprocess <- function(MethyLumiSet,
                       norm = c("none", "illumina", "SWAN", "dasen"), bg.corr = TRUE){

  norm <- match.arg(norm)
  if(norm == "none"){
    if(bg.corr == TRUE){
      object <- methylumi.bgcorr(MethyLumiSet, method = 'illumina', controls = negctls(controlData(MethyLumiSet)))
    }
  }
  if(norm == "illumina"){
    object <- normalizeMethyLumiSet(MethyLumiSet)
    if(bg.corr == TRUE){
      object <- methylumi.bgcorr(object, method = 'illumina', controls = negctls(controlData(object)))
    }
  }
  else if(norm == "SWAN"){
    object <- swan(MethyLumiSet, return.MethylSet = TRUE)
  }
  else if(norm == "dasen"){
    object <- dasen(MethyLumiSet, ret2 = TRUE)
  }
  object
}

subsetProbes <- function(object, allele = c('M', 'U'), type = c('I-red', 'I-green', 'II'),
                         cg = TRUE, snps = TRUE, idmr = TRUE, ch = FALSE){
  data(IlluminaHumanMethylation450kmanifest, envir = environment())
  allele <- match.arg(allele)
  type <- match.arg(type)

  if(cg){
    I.red.names <- minfi::getProbeInfo(IlluminaHumanMethylation450kmanifest, "I-Red")$Name
    I.grn.names <- minfi::getProbeInfo(IlluminaHumanMethylation450kmanifest, "I-Green")$Name
    II.names <- minfi::getProbeInfo(IlluminaHumanMethylation450kmanifest, "II")$Name
  }
  else I.red.names <- I.grn.names <- II.names <- NULL

  if(snps){
    snpI.red.names <- S4Vectors::subset(getProbeInfo(IlluminaHumanMethylation450kmanifest, "SnpI"), Color == "Red")$Name
    snpI.grn.names <- S4Vectors::subset(getProbeInfo(IlluminaHumanMethylation450kmanifest, "SnpI"), Color == "Grn")$Name
    snpII.names <- minfi::getProbeInfo(IlluminaHumanMethylation450kmanifest, "SnpII")$Name
    #print(subset(getProbeInfo(IlluminaHumanMethylation450kmanifest,'SnpII'), Color=='Red'))
    #snpI.grn.names <- snpI.red.names <- snpII.names <- NULL
  }
  else snpI.red.names <- snpI.grn.names <- snpII.names <- NULL

  if(!idmr){
    data(iDMR, envir = environment())
    I.red.names <- I.red.names[!I.red.names %in% iDMR]
    I.grn.names <- I.grn.names[!I.grn.names %in% iDMR]
    II.names <- II.names[!II.names %in% iDMR]
  }

  if(!ch){
    II.names <- II.names[!grepl('^ch', II.names)]
  }

  if(is(object, "MethylSet")){
    if(allele == 'M') mat <- minfi::getMeth(object)
    else mat <- minfi::getUnmeth(object)

    if(type == 'I-red') mat <- mat[rownames(mat) %in% c(I.red.names, snpI.red.names),,drop=FALSE]
    else if(type == 'I-green') mat <- mat[rownames(mat) %in% c(I.grn.names, snpI.grn.names),,drop=FALSE]
    else mat <- mat[rownames(mat) %in% c(II.names, snpII.names),,drop=FALSE]
  }
  else if(is(object, "MethyLumiSet")){
    if(allele == 'M') mat <- methylumi::methylated(object)
    else mat <- methylumi::unmethylated(object)

    if(type == 'I-red') mat <- mat[rownames(mat) %in% c(I.red.names, snpI.red.names),,drop=FALSE]
    else if(type == 'I-green') mat <- mat[rownames(mat) %in% c(I.grn.names, snpI.grn.names),,drop=FALSE]
    else mat <- mat[rownames(mat) %in% c(II.names, snpII.names),,drop=FALSE]
  }

  mat
}


prepSamp <- function(sample){
  if(is.null(sample)) sample <- NULL
  else{
    sample[sample < 0] <- 0
    sample <- sample + 1
    sample <- log2(sample)
  }
  sample
}


