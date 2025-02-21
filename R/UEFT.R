# Functions created for user
# Exact functional functions include EFT, UEFT
# Yiyi Li
# 2023-10-11

# Uniform exact function test
#' @name UEFT
#' @title Uniform Exact Functional Test on Two Discrete Random Variables
#' @aliases UEFT
#' @description
#' Perform the uniform exact functional test on a contingency table to determine if the
#' column variable is a function of the row variable.
#' @usage UEFT(input, correct, log.p)
#'
#' @param input A matrix of nonnegative integers representing a contingency table.
#' Column is the casual and row is the effect.
#' @param correct Logical; if implement the continuity correction.
#' The description is at details. The default is TRUE.
#' @param log.p Logical; if TRUE, the p-value is given as log(p). The default is FALSE. The default is FALSE.
#' @return
#' The exact p-value of the test.
#' @details
#' The uniform idea was implementated using uniform marginal distribution of a square table as null hypothesis.
#'
#' @author Yiyi Li, Joe Song
#' @useDynLib UniExactFunTest
#' @import Rcpp
#' @export
#' @note
#' The functions provide a direct entry into the C++ implementations of the exact functional test.
#'
#' @examples
#'  # Initial a table
#'  x = matrix(c(0,5,10,0,0,5), ncol=3)
#'  # With continuity correction
#'  UEFT(x)
#'  # Without continuity correction
#'  UEFT(x, FALSE)
#' @keywords nonparametric
#' @keywords exact test
#' @keywords functional
#' @keywords uniform
UEFT = function(input, correct = TRUE, log.p = FALSE){
  DNAME <- deparse(substitute(input))
  if(correct){
    method.text <- "Uniform Exact functional test with continuity correction"
  }else{
    method.text <- "Uniform Exact functional test"
  }

  theSum = sum(input)

  if( is.na(theSum) | is.null(theSum) | sum(input%%1!=0)>=1) {
    # Check whether numbers in x are all integers
    stop("ERROR: Exact test requires integer contingency tables!", call. = TRUE)
  }
  # Check sample size
  fun.chisq = getFunChisqStat(input)
  col.sums = colSums(input)
  colChisq = ifelse(theSum > 0, sum(col.sums^2) / theSum * ncol(input) - theSum, 0)

  maxFunChisq = theSum * ncol(input) - theSum - colChisq

  if(maxFunChisq == 0){
    estimate = 0
  }else{
    estimate = sqrt(abs(fun.chisq) / maxFunChisq)
  }


  names(fun.chisq) <- "statistic"

  if(theSum >= 150){
    warning("UEFT requires sample size < 150. Switch to Funchisq")
    df = (nrow(input)-1) * (ncol(input)-1)
    p.value = stats::pchisq( fun.chisq, df = df, lower.tail=FALSE, log.p=log.p )

    result = structure(list(statistic = fun.chisq, p.value = p.value, estimate = estimate,
                            data.name = DNAME, method = method.text),
                       class = "htest")
    return( result )
  }

  # Do test
  if(correct){
    p.value <- UniEFTC(input)
  }else{
    p.value <- UniEFT(input)
  }
  # Take log
  if(log.p) p.value <- log(p.value)

  names(fun.chisq) <- "statistic"
  return(structure(list(statistic = fun.chisq, p.value = p.value, estimate = estimate,
                        data.name = DNAME, method = method.text),
                   class = "htest"))
}
