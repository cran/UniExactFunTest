# Functions created for user
# Exact functional functions include EFT, UEFT
# Yiyi Li
# 2023-09-10

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
#' The continuity correction algorithm
#'
#' @author Yiyi Li, Joe Song
#' @useDynLib UniExactFunTest
#' @import Rcpp
#' @export
UEFT = function(input, correct = TRUE, log.p = FALSE){
  DNAME <- deparse(substitute(input))
  if(correct){
    method.text <- "Uniform Exact functional test with continuity correction"
  }else{
    method.text <- "Uniform Exact functional test"
  }

  if(sum(input%%1!=0)>=1) {
    # Check whether numbers in x are all integers
    stop("ERROR: Exact test requires integer contingency tables!", call. = TRUE)
  }
  # Check sample size
  fun.chisq = getFunChisqStat(input)
  names(fun.chisq) <- "statistic"

  if(sum(input) >= 150){
    warning("UEFT requires sample size < 150. Switch to Funchisq")
    df = (nrow(input)-1) * (ncol(input)-1)
    p.value = stats::pchisq( fun.chisq, df = df, lower.tail=FALSE, log.p=log.p )

    result = structure(list(statistic = fun.chisq, p.value = p.value, estimate = fun.chisq,
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
  return(structure(list(statistic = fun.chisq, p.value = p.value, estimate = fun.chisq,
                        data.name = DNAME, method = method.text),
                   class = "htest"))
}
