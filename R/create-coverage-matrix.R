#' Create normalized data frame from  input coverage bed files 
#' The function will read in bed files, coverage column will be specified by user
#' and collate them to make a data matrix and then normalize by the library size. 
#' If the library sizes are not provided by then use sum of reads across loci as proxy. 
#' Though this is not recommended.
#' 
#' @param inputPath Path to coverage or junction bed files. 
#' Assumptions for coverage bed file - 1st 3 columns are the loci of interest (chr start stop) and the last column is the coverage readcount
#' 
#' @param filePattern character string or regular expression to search inputPath to get all coverage files
#' 
#' @param junctionReads logical flag specifying whether junction reads are to be calculated OR not
#' 
#' @param junctionType character type of the junction bed file. Either 'tophat' or 'star'
#' 
#' @param librarySizeFile full path file name with library sizes. Should have two columns
#' 1st should be sample name - same as file name in @param inputPath and the 2nd as the library size. No of rows should be equal to
#' no of coverage files read from @param inputpath 
#' 
#' @param promoterFile Full path file name of promoter loci - needs to be a bed file if junctionReads flag is equal to TRUE
#' @param normalize logical flag specifying whether normalization by library size should be done or not
#' 
#' @return A data.frame object with normalized RNAseq reads in each sample (rows) 
#' @author Aditi Qamra
#' @export
#'
#' @examples
#' \dontrun{
#' norm_matrix <- createCoverageMatrix(inputPath="readcount/",
#'                                    filePattern="coverage_",
#'                                    junctionReads=FALSE,
#'                                    junctionType='',
#'                                    promoterFile='',
#'                                    normalize=TRUE,
#'                                    librarySizeFile="libsizes.txt")
#' 
#'norm_matrix <- createCoverageMatrix(inputPath="test/",
#'                                    filePattern=".junction.bed",
#'                                    junctionReads=TRUE,
#'                                    junctionType='tophat',
#'                                    promoterFile='test/aploci.bed',
#'                                    normalize=TRUE,
#'                                    librarySizeFile="libsizes.txt")
#' 
#' 
#' 
#' }
#'
#'


createCoverageMatrix=function(inputPath,
                                  filePattern,
                                  junctionReads,
                                  junctionType, 
                                  promoterFile,
                                  normalize,
                                  librarySizeFile,... ){
  
  if (missing(inputPath)) {
      stop(paste0('Error: inputPath of coverage or junction files is missing' ))
  }
  
  
  if (missing(filePattern)) {
    stop(paste0('Error: string to identify coverage or junction files in inputPath is missing' ))
  }
  
  
  if (missing(junctionReads)) {
    stop(paste0('Error: logical flag specifying whether junction reads are being read in or not is missing. ' ))
  }
  
  if (junctionReads==TRUE & !junctionType %in% c('tophat', 'star')) {
    stop(paste0('Error: Invalid junction type: ', junctionType, '! Possible values: "tophat" or "star"'))
  }  

  if (missing(normalize)) {
    stop(paste0('Error: logical flag to normalize count matrix is missing' ))
  }  
  
  if (missing(librarySizeFile)) {
    paste0('no library size file was provided. 
           sum of read counts across loci will be used as proxy' )
  }
  
  
  if (junctionReads!=TRUE){
    
    # read in files
    
    files <- dir(inputPath, pattern=filePattern, all.files=T, full.names=T)
    
    # combine
    
    promoterReadCounts <- data.table::fread(files[1],sep="\t", stringsAsFactors = F, header=F, data.table=F)
    colnames(promoterReadCounts)[ncol(promoterReadCounts)] <- gsub(filePattern, "",basename(files[1]))
    promoterReadCounts$name <- paste(promoterReadCounts[,1], promoterReadCounts[,2], promoterReadCounts[,3], sep="_")
    promoterReadCounts <- promoterReadCounts[,c("name",  gsub(filePattern, "",basename(files[1])))]
    
    for (f in 2:length(files)){
      dat <- fread(files[f],sep="\t", stringsAsFactors = F, header=F, data.table=F)
      promoterReadCounts <- cbind(promoterReadCounts, dat[,ncol(dat)])
      colnames(promoterReadCounts)[ncol(promoterReadCounts)] <- gsub(filePattern, "", basename(files[f]))
    }
    
    promoterReadCounts <- promoterReadCounts[!duplicated(promoterReadCounts$name), ]
    rownames(promoterReadCounts) <- promoterReadCounts$name
    promoterReadCounts$name <- NULL
    
  } else {
    
    files <- dir(inputPath, pattern=filePattern, all.files=T, full.names=T)
    
    promoterReadCounts <- as.data.frame(lapply(files,calculateJunctionReadCounts,
                                  junctionType=junctionType,
                                 promoterFile=promoterFile))
    
    colnames(promoterReadCounts) <- gsub(filePattern,"", colnames(promoterReadCounts))
  }
  

  
  if (normalize!=TRUE){
    
    return(as.data.frame(promoterReadCounts))
    
  } else {
    
    # Normalize by library size and region size
    
    dist_vec <- as.numeric(str_split_fixed(rownames(promoterReadCounts), "_",3)[,3])-as.numeric(str_split_fixed(rownames(promoterReadCounts), "_",3)[,2])
    
    # region size
    promoterReadCounts_norm <- apply(promoterReadCounts,2, function(e) e/dist_vec)
    
    # libsize
    
    if (!missing(librarySizeFile)) {
      libsize <- read.table(librarySizeFile,sep="\t", stringsAsFactors = F, header=F)
      rownames(libsize) <- gsub(filePattern,"", libsize[,1] )
      libsize <- libsize[colnames(promoterReadCounts_norm),]
      promoterReadCounts_norm <- apply(promoterReadCounts,2, function(e) e/libsize)
    } else {
      
      libsize <- apply(promoterReadCounts,2, sum)
      promoterReadCounts_norm <- as.data.frame(t(apply(promoterReadCounts_norm,1, function(e) e*1000000/libsize)))
    }
    
    # return data frame
    return(as.data.frame(promoterReadCounts_norm))
    
  }
 
} 

