## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(UniExactFunTest)

## ----test, include=FALSE, eval=TRUE-------------------------------------------
xs <- 
  list(
    matrix(c(
      0, 5,
      5, 0
    ), nrow=2, byrow=TRUE),
    matrix(c(
      0, 5,
      3, 5
    ), nrow=2, byrow=TRUE),
    matrix(c(
      0, 5,
      12, 3
    ), nrow=2, byrow=TRUE)
  )

for(x in xs) {
  print(knitr::kable(x))
  UEFT(x)
  UEFT(t(x))

}

## ----ice cream----------------------------------------------------------------
x <- matrix(c(
  13, 32,
  17, 23
), nrow=2, byrow=TRUE,
dimnames = list(
  c("Ice cream eaten", "No ice cream"),
  c("Cases", "Controls")
))
knitr::kable(x)
UEFT(x)
UEFT(t(x))


## ----AIDS---------------------------------------------------------------------
x <- matrix(
  c(4,2,3,3,16,2), nrow=2, ncol=3, byrow=TRUE,
  dimnames=list(
              c("AIDS: yes", "AIDS: no"),
              c("Males", "Females", "Both")
))
knitr::kable(x)
UEFT(x)
UEFT(t(x))

## ----Mendenhall et al. 1984---------------------------------------------------
x <- matrix(c(
  21, 2,
  15, 3
), nrow=2, byrow=TRUE,
dimnames = list(
  c("Surgery", "Radiation therapy"),
  c("Cancer controlled", "Cancer not controlled")
))
knitr::kable(x)
UEFT(x)
UEFT(t(x))

## ----sex-handedness-----------------------------------------------------------
x <- matrix(c(
  43, 9,
  44, 4
), nrow=2, byrow=TRUE,
dimnames = list(
  c("Male", "Female"),
  c("Right-handed", "Left-handed")
))
knitr::kable(x)
UEFT(x)
UEFT(t(x))

