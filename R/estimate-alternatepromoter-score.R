#' Calculate AP score given a readcount matrix and K4me3 annotation
#'
#' @param promoterReadCounts A data.frame object. The number of reads
#'   per promoter (rows) for each sample (cols). Rownames should be promoter names
#'
#' @param promoterAnnotation A list with elements. Number of elements
#'  equal to nrow (promoterReadCounts) and order same as promoterReadCounts
#'  indicating directionality i.e. 'gain' or 'loss'
#'  of promoter promoterReadCounts
#'
#' @param promoterMethod Character type.
#'   Identifying method of AP score calculation
#'   Can be either "medianbased" or "rankbased"
#'
#' @param medianThreshold A numeric value.
#'   Specify if promoterMethod is equal to "medianbased"
#'
#' @return A data.frame object. An AP score (column) calculated for each sample (rows) using  @param promoterMethod
#' @export
#'
#' @examples
#' \dontrun{
#' alternatePromoterScore <- calculateAlternatePromoterScore(promoterReadCounts,
#'                                                          promoterAnnotation,
#'                                                          promoterMethod = "medianbased",
#'                                                          medianThreshold = 4)
#' }


calculateAlternatePromoterScore=function(promoterReadCounts,
                                         promoterAnnotation,
                                         promoterMethod=c("medianbased", "rankbased"),
                                         medianThreshold,... ){

  if (!(promoterMethod == 'medianbased' | promoterMethod==  'rankbased')) {
    stop(paste0('Error: Invalid method type: ', promoterMethod, '! Possible values: "medianbased" or "rankbased"'))
  }

  if (!(is.data.frame(promoterReadCounts) | is.matrix(promoterReadCounts))) {
    stop(paste0('Error: Invalid data type: ', promoterReadCounts, '! Possible values: "data.frame" or "matrix" '))
  }

  if (missing(promoterAnnotation)) {
    promoterAnnotation <- c(rep("gain"), nrow(promoterReadCounts))
    paste0('Setting all promoter annotations to "gain"')
  } else if (!(is.vector(promoterAnnotation))){
    stop(paste0('Error: Invalid data type: ', promoterAnnotation, '! Possible values: "character vector" '))
  } else if (length(promoterAnnotation) != nrow(promoterReadCounts) ){
    stop(paste0('Error: ', promoterAnnotation, 'does not have elements equal to nrow(', promoterReadCounts,')'  ))
  }


  promoterReadCounts <- as.data.frame(promoterReadCounts)

  # Calculate APscore
  if (promoterMethod=="rankbased"){

    apgain <- promoterReadCounts[(tolower(promoterAnnotation)=="gain"),]
    gainscore <- as.data.frame(t(apply(apgain,1,function(e) rank(e))))

    if (unique(tolower(promoterAnnotation)=="loss")!=FALSE){

      aploss <- promoterReadCounts[(tolower(promoterAnnotation)=="loss"),]
      lossscore <- as.data.frame(t(apply(aploss,1,function(e) rank(e*-1))))

    }

  } else if (promoterMethod=="medianbased"){

    apgain <- promoterReadCounts[(tolower(promoterAnnotation)=="gain"),]
    apgainmedian <-  apply(apgain,1,median,na.rm = T)
    gainscore <- as.data.frame( apply(apgain,2,function(e) if_else(e>=(apgainmedian*medianThreshold),1,0) ))

    if (unique(tolower(promoterAnnotation)=="loss")!=FALSE){
      aploss <- promoterReadCounts[(tolower(promoterAnnotation)=="loss"),]
      aplossmedian <-  apply(aploss,1,median,na.rm = T)
      lossscore <- as.data.frame( apply(aploss,2,function(e) if_else(e<=(aplossmedian*(1/medianThreshold)),1,0) ))
    }
  }

  apscore <- data.frame(id=colnames(promoterReadCounts), stringsAsFactors=F)
  apscore$gainapscore <- apply(gainscore,2,sum, na.rm = T)

  if (unique(tolower(promoterAnnotation)=="loss")!=FALSE){
      apscore$lossapscore <- apply(lossscore,2,sum, na.rm = T)
      apscore$apscore <- apscore$gainapscore + apscore$lossapscore
      apscore <- dplyr::arrange(apscore, apscore)
      apscore <- apscore[,c("id", "apscore")]
  } else {

    apscore$apscore <- apscore$gainapscore
    apscore <- dplyr::arrange(apscore, apscore)
    apscore <- apscore[,c("id", "apscore")]
  }

  return(apscore)

}

