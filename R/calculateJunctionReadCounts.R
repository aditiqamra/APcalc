#' Calculate sum of 1st exon-intron junction reads overlapping input promoter loci
#' 
#' @param junctionFile Path to junction files either from STAR or tophat
#' 
#' @param junctionType character type of the junction bed file. Either 'tophat' or 'star'
#' 
#' @param promoterFile full path file name of the promoter loci - should be a bed file. chr start stop as the 1st 3 columns in zero based format
#' 
#' @return A data.frame object with 1st exon junction reads for each promoter loci
#' @author Aditi Qamra
#' @export
#'
#' @examples
#' \dontrun{
#' junction_reads <- calculateJunctionReadCounts(junctionFile="test/abc.junction.bed",
#'                                    junctionType="tophat",
#'                                    promoterFile="test/aploci.bed")
#' 
#' 
#' 
#' 
#' }
#'
#'

calculateJunctionReadCounts <- function(junctionFile='', 
                                        junctionType,
                                        promoterFile,
                                        gencodeFile) {

  
  if (!junctionType %in% c('tophat', 'star')) {
    stop(paste0('Error: Invalid junction type: ', junctionType, '! Possible values: "tophat" or "star"'))
  }
  
  if (is.null(junctionFile)) {
    stop(paste0('Error: Please specify valid junction file path!'))
  }
  
  if (missing(promoterFile)) {
    stop(paste0('Error: Please specify valid promoter file path!'))
  }

  if(junctionType == 'tophat') {
    print(paste0('Reading Tophat junction files from: ', junctionFile))
    junctionTable <- readTopHatJunctions(junctionFile)
    seqlevelsStyle(junctionTable) <- 'UCSC'
    
    print('File loaded into memory')
    
  } else if(junctionType == 'star') {
    
    print(paste0('Reading STAR junction files from', junctionFile))
    junctionTable <- readSTARJunctions(junctionFile)
    seqlevelsStyle(junctionTable) <- 'UCSC'
    
    junctionTable$score <- junctionTable$um_reads  # to match the tophat style, uniquely mapped reads are used as score
    print('File loaded into memory')
  }
  
  print('Identifying 1st exon-intron junctions ')
  junctionTable.overlap <- GenomicAlignments::findOverlaps(junctionTable, gencodeFile, type = 'equal')
  firstintronjunctionTable <- junctionTable[queryHits(junctionTable.overlap)]
  
  print('Calculating junction counts')
  promoterloci <- read.table(promoterFile, header=F, sep="\t", stringsAsFactors = F)[,1:3]
  colnames(promoterloci) <- c("chromosome","start", "end")
  promoterloci <- makeGRangesFromDataFrame(promoterloci, starts.in.df.are.0based=TRUE)
  
  promoter.overlap <- GenomicAlignments::findOverlaps(promoterloci ,firstintronjunctionTable, ignore.strand=TRUE)
  promoterloci$junctionCounts <- rep(0, length(promoterloci))
  promoterloci$junctionCounts[queryHits(promoter.overlap)] <- firstintronjunctionTable$score[subjectHits(promoter.overlap)]
  promoterloci$name <- paste(seqnames(promoterloci), start(promoterloci)-1, end(promoterloci), sep="_")
  junctionCounts <- data.frame(stringsAsFactors=F, 
                               junctioncount=tapply(promoterloci$junctionCounts, as.factor(promoterloci$name), sum))
  colnames(junctionCounts) <- basename(junctionFile)
  return(as.data.frame(junctionCounts))
}

