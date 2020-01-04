#' Given alternate promoter score for samples,
#' divide samples into groups and calculate survival differences
#'
#' @param alternatePromoterScore A data.frame object. Two columns containing sample ids and
#' numeric alternate promoter score
#'
#' @param clinicaldata A data.frame object. Should have an id column labelled "id"
#' and have survival time and event information
#'
#' @param alternatePromoterUpperQuantileThreshold Numeric value greater than 0 and less than 1.
#'   Identifying threshold to define samples into APhigh and APlow.
#'   APhigh samples will be alternatePromoterScore >= quantile(alternatePromoterScore,alternatePromoterUpperQuantileThreshold).
#'   APlow samples will be alternatePromoterScore <= quantile(alternatePromoterScore,(1-alternatePromoterUpperQuantileThreshold)).

#' @param survivalTime Column in clinicaldata specifying survival time
#'
#' @param survivalEvent Column in clinicaldata specifying event information
#'
#' @param title Character type. Label to plot on survival curve
#'
#' @return A data.frame object. Return alternatePromoterScore with AP group classification
#' @export
#'
#' @examples
#' \dontrun{
#' alternatePromoterGroups <- calculateSurvival(alternatePromoterScore,
#'                                                          clinicaldata,
#'                                                          alternatePromoterUpperQuantileThreshold=0.66,
#'                                                          survivalTime="pfs",
#'                                                          survivalEvent="progression",
#'                                                          title="AP Group -  66/33 Split")
#' }
#'
#' @seealso \code{\link{calculateAlternatePromoterScore}} for calculating
#' alternate promoter score for samples
#'


calculateSurvival=function(alternatePromoterScore,
                           clinicaldata,
                           alternatePromoterQuantileThreshold,
                           survivalTime,
                           survivalEvent,
                           title,...){

  if (!(is.numeric(alternatePromoterQuantileThreshold))) {
    stop(paste0('Error: Invalid type: ', alternatePromoterQuantileThreshold, '! Should be a numeric value'))
  }

  if (!(between(alternatePromoterQuantileThreshold,0,1))) {
    stop(paste0('Error: Invalid range for : ', alternatePromoterQuantileThreshold, '! Should be between 0 and 1'))
  }

  if (!(is.character(title))) {
    stop(paste0('Error: Invalid type: ', title, '! Should be character string '))
  }


  if (!(is.data.frame(alternatePromoterScore) | is.matrix(alternatePromoterScore))) {
    stop(paste0('Error: Invalid data type: ', promoterReadCounts, '! Possible values: "data.frame" or "matrix" '))
  }

  alternatePromoterScore <- as.data.frame(alternatePromoterScore)

  if (!(is.data.frame(clinicaldata) | is.matrix(clinicaldata))) {
    stop(paste0('Error: Invalid data type: ', clinicaldata, '! Possible values: "data.frame" or "matrix" '))
  }

  if(!("id" %in% colnames(clinicaldata))){
    stop(paste0("Error: Clinical data should have column named 'id' mapping to id column of alternatePromoterScore "))
  }


  clinicaldata <- as.data.frame(clinicaldata)
  meta <- dplyr::left_join(alternatePromoterScore, clinicaldata, by = "id")

  meta$apgroup <- ifelse(meta$apscore >= quantile(meta$apscore, alternatePromoterQuantileThreshold),"High",
                         ifelse(meta$apscore <= quantile(meta$apscore, (1-alternatePromoterQuantileThreshold)), "Low","Middle"))

  form <-  paste0( "Surv(time=", survivalTime,",event=",survivalEvent,") ~apgroup")
  print(head(meta))
  survfitfile <- survfit(as.formula(form), data = meta)
  
  p <- survfitfile %>% ggsurvplot(survfitfile, data=meta, pval = TRUE, risk.table = TRUE, risk.table.height = 0.34, surv.plot.height = 1, palette = "jco", main=title)
  print(p)
  return(survfitfile)
  
}
