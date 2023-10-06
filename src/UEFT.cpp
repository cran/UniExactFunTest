#include "UEFT.h"
#include "fastEnu.h"

/////////////////////////////
// Uni exact functional test
double getUniEFT(vector<vector<int>> inputm, bool secondIm, bool iftrim, bool getCDF)
{
  if (iftrim == true)
  {
    inputm = trimTable(inputm);
  }
  unsigned r = inputm.size();

  if (r == 0 || r == 1)
  {
    return 1;
  }
  unsigned s = inputm[0].size();

  if (s == 0 || s == 1)
  {
    return 1;
  }

  if (getchisqInteger(inputm) <= 0)
  {
    return 1;
  }

  int N = 0;
  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < s; j++)
    {
      N += inputm[i][j];
    }
  }

  int valueForFac;
  valueForFac = (r > N ? r : N);
  valueForFac = (s > valueForFac ? s : valueForFac);

  vector<double> factorialNum = preFactorial(valueForFac);

  double oriChi;

  double UniEFT;
  vector<int> rowSum, colSum;

  if (r >= s)
  {
    if (secondIm)
    {
      colSum = disTable(N, s, s);
      rowSum = disTable(N, s, s);
      // Another implementation
      // colSum = disTable(N, s, s);
      // rowSum = disTable(N, s, r);
    }
    else
    {
      colSum = disTable(N, s, s);
      rowSum = disTable(N, r, r);
    }
    if (getCDF)
    {
      UniEFT = fastEnu::UEFTNetwork(rowSum, colSum, N, inputm, factorialNum, true);
    }
    else
    {
      UniEFT = fastEnu::UEFTNetwork(rowSum, colSum, N, inputm, factorialNum);
    }
  }
  else
  {
    colSum = disTable(N, r, s);
    rowSum = disTable(N, r, r);
    if (getCDF)
    {
      UniEFT = fastEnu::UEFTNetwork(rowSum, colSum, N, inputm, factorialNum, true);
    }
    else
    {
      UniEFT = fastEnu::UEFTNetwork(rowSum, colSum, N, inputm, factorialNum);
    }
  }

  return UniEFT;
}

// The UEL not expand the sample size, but use approximatly P-value among all cases.
double getUniEFTC(const vector<vector<int>> &inputm, bool getCDF)
{
  unsigned r = inputm.size();

  if (r == 0 || r == 1)
  {
    return 1;
  }
  unsigned s = inputm[0].size();

  if (s == 0 || s == 1)
  {
    return 1;
  }

  int N = 0;
  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < s; j++)
    {
      N += inputm[i][j];
    }
  }

  if (N == 0)
  {
    return 1;
  }

  vector<int> colSum, rowSum;
  if (r >= s)
  {
    colSum = disTable(N, s, s);
    rowSum = disTable(N, s, s);
    // Another implementation
    // colSum = disTable(N, s, s);
    // rowSum = disTable(N, s, r);
  }
  else
  {
    colSum = disTable(N, r, s);
    rowSum = disTable(N, r, r);
  }

  int valueForFac;
  valueForFac = (r > N ? r : N);
  valueForFac = (s > valueForFac ? s : valueForFac);
  vector<double> factorialNum = preFactorial(valueForFac);
  THREEP P3;

  P3 = fastEnu::UEFTCNetwork(rowSum, colSum, N, inputm, factorialNum, getCDF);
  double result;

  if (P3.stat <= 0)
  {
    return 1;
  }
  // Cases1: CDF, stat at enumerated table
  if (getCDF && P3.stat == P3.rightS)
  {
    result = P3.rightP;
  }
  // Cases2: CDF, stat at left enumerated table
  else if (getCDF && P3.stat == P3.leftS)
  {
    result = P3.midP;
  }
  // Case3: CDF, normal case
  else if (getCDF)
  {
    // No need leftP.
    // midP: left bottom, right top
    // rightP; right bottom.
    result = (P3.stat - P3.leftS) *
                 (P3.midP - P3.rightP) / (P3.rightS - P3.leftS) +
             P3.rightP;
  }
  // Pvalue case
  else
  {
    // midP: Leftbottom Pvalue and righttop Pvalue. leftP: leftTop pvalue. rightP: rightbottom Pbalue
    double pR, pL;
    pR = (P3.midP + P3.rightP) / 2;
    pL = (P3.midP + P3.leftP) / 2;
    // Case1: right S and left S same
    if (P3.rightS == P3.leftS)
    {
      result = pR;
    }
    // Case2: Normal case
    else
    {
      result = (P3.stat - P3.leftS) *
                   (pR - pL) / (P3.rightS - P3.leftS) +
               pL;
    }
  }
  // Return
  return result;
}
