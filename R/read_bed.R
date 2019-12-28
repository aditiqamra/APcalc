#' ReadBed function from package genomation - https://github.com/BIMSBbioinfo/genomation/blob/master/R/readData.R
#' 
#' @return a granges object
#' @export
#'
#'

read_bed<-function(file,track.line=FALSE,remove.unusual=FALSE,zero.based=TRUE)
{
  
  meta.cols=list(score=5,name=4,thickStart=7,  
                 thickEnd=8, 
                 itemRgb=9,  
                 blockCount=10, 
                 blockSizes=11,  
                 blockStarts=12 )
  
  file <- compressedAndUrl2temp(file)
  if(is.numeric(track.line)){
    skip = track.line
  }else if(track.line=="auto"){
    skip = detectUCSCheader(file)
  }else{
    skip = 0
  }
  
  df=read_delim(file,skip=skip,n_max=2,col_names=FALSE, delim="\t")
  numcol=ncol(df)
  
  if(numcol==3){
    df=readGeneric(file, chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = NULL,   zero.based = zero.based,
                   remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t")
  }else if(numcol==4){
    df=readGeneric(file, chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = meta.cols[1],   zero.based = zero.based,
                   remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t")    
  }else if(numcol==5){
    df=readGeneric(file, chr = 1, start = 2, end = 3, strand = NULL,
                   meta.cols = meta.cols[1:2],   zero.based = zero.based,
                   remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t")    
  }else if(numcol == 6){
    df=readGeneric(file, chr = 1, start = 2, end = 3, strand = 6,
                   meta.cols = meta.cols[1:2],   zero.based = zero.based,
                   remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t") 
  }else if(numcol > 6){
    df=readGeneric(file, chr = 1, start = 2, end = 3, strand = 6,
                   meta.cols = meta.cols[c(1:2,(3:8)[1:(numcol-6)])],   zero.based = zero.based,
                   remove.unusual =remove.unusual, header = FALSE, skip = skip, sep = "\t") 
  }
  df
  
} 
