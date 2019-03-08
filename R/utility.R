# utility_function ----
# u <- function(x, alpha = 0.5){
#   -exp(-alpha*x) %>% return()
# }
#
#' Power utility
#'
#' @author Yangzhuoran Yang
#' @export
power_u <- function(x, lambda = 0.5){
  x^(1-lambda)/(1-lambda)
}
