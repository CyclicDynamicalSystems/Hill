source("hill-core.R")
library(RUnit)

hill.test3d <- function() {
  h <- hill(3, 18, 3)
  
  checkEquals(h$n, 3)
  checkEquals(h$a, 18)
  checkEquals(h$b, 1)
  checkEquals(h$g, 3)
              
  checkEquals(h$equil$x, c(2, 2, 2))
  checkEquals(h$f(2), 2)
  checkEquals(h$df(2), -2.666666666666666666666)
}

hill.test3d()