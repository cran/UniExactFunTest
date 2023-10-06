# Code to test UEFT
library(testthat)
library(UniExactFunTest)

test_that("TestEFTAndUEFT", {
#
  mList = list()
  sList = list()
  tsList = list()
  UEFTP = list()
  UEFTRP = list()
  UEFTCP = list()
  UEFTCRP = list()

  # Test1
  mList[[1]] = matrix(c(12, 26, 18, 0, 8, 12), nrow=2, byrow = TRUE)
  sList[[1]] = 0.042556227
  tsList[[1]] = 0.027271581

  UEFTP[[1]] = 0.03825917
  UEFTRP[[1]] = 0.01112863
  UEFTCP[[1]] = 0.06294452 # should be 0.06294452
  UEFTCRP[[1]] = 0.02151119


  # Test2
  mList[[2]] = matrix(c(0,0,0,0,0,0,0,0,0), nrow=3, byrow = TRUE)
  sList[[2]] = 1
  tsList[[2]] = 1

  UEFTP[[2]] = 1
  UEFTRP[[2]] = 1
  UEFTCP[[2]] = 1
  UEFTCRP[[2]] = 1


  # Test3
  mList[[3]] = matrix(c(4,0,4,0,4,0,1,0,1), 3)
  sList[[3]] = 0.002997003
  tsList[[3]] = 0.0065490065

  UEFTP[[3]] = 0.004170433
  UEFTRP[[3]] = 0.008134722
  UEFTCP[[3]] = 0.003717371
  UEFTCRP[[3]] = 0.008729366


  # Test4
  mList[[4]] = matrix(rep(7,16), nrow=4)
  sList[[4]] = 1
  tsList[[4]] = 1

  UEFTP[[4]] = 1
  UEFTRP[[4]] = 1
  UEFTCP[[4]] = 1
  UEFTCRP[[4]] = 1


  # Test 5
  mList[[5]] = matrix(c(4,0,0,0,4,0,0,0,4), nrow=3, byrow = TRUE)
  sList[[5]] = 0.00017316017
  tsList[[5]] = 0.00017316017


  UEFTP[[5]] = 0.0001731602
  UEFTRP[[5]] = 0.0001731602
  UEFTCP[[5]] = 8.658009e-05
  UEFTCRP[[5]] = 8.658009e-05



  mList[[6]] = matrix(c(2,0,0,2), nrow=2, byrow = TRUE)
  sList[[6]] = 0.3333333333
  tsList[[6]] = 0.3333333333

  UEFTP[[6]] = 0.3333333
  UEFTRP[[6]] = 0.3333333
  UEFTCP[[6]] = 0.1666667
  UEFTCRP[[6]] = 0.1666667


  mList[[7]] = matrix(c(0,0,0,5,
                        5,1,0,0,
                        0,0,0,5,
                        5,0,1,0),
                      nrow=4, byrow = TRUE)
  sList[[7]] = 1.443345e-05
  tsList[[7]] = 1.443345e-05

  UEFTP[[7]] = 3.924468e-06
  UEFTRP[[7]] = 0.001089909
  UEFTCP[[7]] = 3.853827e-06
  UEFTCRP[[7]] = 0.001052053



  mList[[8]] = matrix(c(1,10,15,20,5,2,13,0,0), nrow=3, byrow = TRUE)
  sList[[8]] = 6.521275e-11
  tsList[[8]] = 1.384389e-09

  UEFTP[[8]] = 2.949504e-10
  UEFTRP[[8]] = 2.441368e-09
  UEFTCP[[8]] = 2.836923e-10
  UEFTCRP[[8]] = 2.495195e-09


  mList[[9]] = matrix(c(
    6,    0,    8,    9,    6,
    0,   11,    0,    0,    0,
    0,    0,    0,    0,    0
  ),
  nrow=3, byrow=TRUE
  )
  sList[[9]] = 4.325631e-10
  tsList[[9]] = 4.325631e-10

  UEFTP[[9]] = 1.791214e-06
  UEFTRP[[9]] = 2.998506e-10
  UEFTCP[[9]] = 1.778817e-06
  UEFTCRP[[9]] = 2.945847e-10



  mList[[10]] = matrix(c(
    0, 0, 0, 10, 9,
    0, 0, 0, 0, 0,
    8, 3, 10, 0, 0
  ),
  nrow=3, byrow=TRUE
  )
  sList[[10]] = 1.523433e-11  # Bug Different here, ori EFT got 1.523433e-11
  tsList[[10]] = 1.523433e-11

  UEFTP[[10]] = 1.040541e-05
  UEFTRP[[10]] = 1.449086e-13
  UEFTCP[[10]] = 1.04134e-05
  UEFTCRP[[10]] = 1.66459e-13



  mList[[11]] = matrix(c(
    8,    0,   11,    0,   11,
    0,    4,    0,    0,    0,
    0,    0,    0,    6,    0), nrow=3, byrow=TRUE
  )
  sList[[11]] =  5.617703e-12
  tsList[[11]] = 5.617703e-12

  UEFTP[[11]] = 2.642357e-07
  UEFTRP[[11]] = 2.264492e-10
  UEFTCP[[11]] = 2.670806e-07
  UEFTCRP[[11]] = 2.299357e-10


  # Test 13
  mList[[12]] = matrix( c(
    0, 2, 0,
    0, 2, 0,
    0, 2, 0 ), nrow = 3, byrow =T)
  sList[[12]] = 1
  tsList[[12]] = 1

  UEFTP[[12]] = 1
  UEFTRP[[12]] = 1
  UEFTCP[[12]] = 1
  UEFTCRP[[12]] = 1


  # Test 13
  mList[[13]] = matrix( c(
    1, 2, 1,
    2, 2, 3,
    1, 1, 0
  ), nrow = 3, byrow =T)
  sList[[13]] = 1
  tsList[[13]] = 0.8321678322

  UEFTP[[13]] = 0.7602398
  UEFTRP[[13]] = 0.7602398
  UEFTCP[[13]] = 0.8313115
  UEFTCRP[[13]] = 0.8509183


  mList[[14]] = matrix( c(
    1, 0, 0,
    8, 2, 0,
    0, 1, 1
  ), nrow = 3, byrow =T)
  sList[[14]] = 0.1818181818
  tsList[[14]] = 0.1818181818

  UEFTP[[14]] = 0.3526474
  UEFTRP[[14]] = 0.3526474
  UEFTCP[[14]] = 0.3370732
  UEFTCRP[[14]] = 0.3326104


  mList[[15]] =
    matrix( c(0,    0,    0,    7,
              7,    0,    0,    0,
              0,    0,    0,    20,
              0,    0,    0,    6
    ), nrow = 4, byrow =T)
  sList[[15]] = 1.072756491e-07
  tsList[[15]] = 1.072756491e-07 # Bug, EFT got 1.072756491e-07

  UEFTP[[15]] = 2.469333e-07
  UEFTRP[[15]] = 6.673822e-05
  UEFTCP[[15]] = 2.197891e-07
  UEFTCRP[[15]] = 6.010585e-05


  for(i in 1:length(mList)){
    T = mList[[i]]

    v1 = signif(UEFT(T,FALSE)$p.value, 6)
    v2 = signif(UEFT(t(T), FALSE)$p.value, 6)
    v3 = signif(UEFT(T)$p.value, 6)
    v4 = signif(UEFT(t(T))$p.value, 6)


    expect_equal(v1, signif(UEFTP[[i]], 6) )
    expect_equal(v2, signif(UEFTRP[[i]], 6) )
    expect_equal(v3, signif(UEFTCP[[i]], 6) )
    expect_equal(v4, signif(UEFTCRP[[i]], 6) )
  }

}# end test
)
