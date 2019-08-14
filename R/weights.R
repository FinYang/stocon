
#' Finding weights that minimise variance of return using linear regression
#'
#' Finding weights that minimise variance of return using linear regression as described in
#' Fan, Zhang & Yu (2012)
#'
#'
#'
#' @param Rt Lists of return. e.g. output from \code{\link{sim_simple}}
#'
#' @return List of weights for each period
#' @author Yangzhuoran Yang
#' @importFrom magrittr %>%
#' @references Fan, J., Zhang, J., & Yu, K. (2012). Vast Portfolio Selection With Gross-Exposure Constraints. Journal of the American Statistical Association, 107(498), 592–606. https://doi.org/10.1080/01621459.2012.682825
#' @export
weights_lm <- function(Rt){
  reg_data <- cbind(Rt[,1], Rt[,1] - Rt[,2:NCOL(Rt)])
  model <- lm(reg_data[,1]~reg_data[,-1])
  weights <- c(1-sum(coef(model)[-1]), coef(model)[-1])
  names(weights) <- NULL
  return(weights)
}
