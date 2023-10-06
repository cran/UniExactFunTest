// The UEFT implementation code
#ifndef UEFT_H
#define UEFT_H

// use trimTable, getchisqInteger, preFactorial
#include "matrixTools.h" 
using namespace std;

////////////////////////
// UEFT section
// UniEFT
double getUniEFT(vector<vector<int>>, bool secondIm = true, bool iftrim = false, bool getCDF = false);

// The UEL not expand the sample size, but use approximatly P-value among all cases.
// getCDF: default is false. When true, it will output 1 - CDF
double getUniEFTC(const vector<vector<int>> &, bool getCDF = false);

#endif
