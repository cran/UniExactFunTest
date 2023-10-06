

#include "matrixTools.h"

vector<vector<int>> trimTable(const vector<vector<int>> &table)
{ // DQP::trimTable() removes all-zero rows and all-zero columns
    //   from the input table
    // MS 9/21/2019
    vector<vector<int>> trimmed;

    size_t nrows = table.size();
    size_t ncols = nrows > 0 ? table[0].size() : 0;

    vector<int> rowsums(nrows, 0);
    vector<int> colsums(ncols, 0);

    for (auto i = 0u; i < nrows; ++i)
    {
        for (auto j = 0u; j < ncols; ++j)
        {
            rowsums[i] += table[i][j];
            colsums[j] += table[i][j];
        }
    }

    size_t nrows_trimmed(nrows), ncols_trimmed(ncols);

    for (auto i = 0u; i < nrows; ++i)
    {
        if (rowsums[i] == 0)
            nrows_trimmed--;
    }

    for (auto j = 0u; j < ncols; ++j)
    {
        if (colsums[j] == 0)
            ncols_trimmed--;
    }

    if (nrows_trimmed == nrows && ncols_trimmed == ncols)
    {
        trimmed = table;
    }
    else
    {
        trimmed = vector<vector<int>>(nrows_trimmed, vector<int>(ncols_trimmed));
        size_t i_t = 0u;
        for (auto i = 0u; i < nrows; ++i)
        {
            if (rowsums[i] == 0)
            {
                continue;
            }
            size_t j_t = 0u;
            for (auto j = 0u; j < ncols; ++j)
            {
                if (colsums[j] == 0)
                {
                    continue;
                }
                trimmed[i_t][j_t] = table[i][j];
                j_t++;
            }
            i_t++;
        }
    }

    return trimmed;
}

// Get Funchisq stat for X->Y
double getchisqInteger(const vector<vector<int>> &input)
{
  double sum1, sum2;
  // r, s, matrix are class variables
  sum1 = 0;
  sum2 = 0;

  double sumcol;
  vector<double> sumrow(input.size(), 0);

  double N = 0;
  double r = input.size();
  double s = input[0].size();

  for (unsigned i = 0; i < r; i++)
  {
    for (unsigned j = 0; j < s; j++)
    {
      sumrow[i] += input[i][j];
    }
    N += sumrow[i];
  }

  if (N == 0)
  {
    return 0;
  }

  for (int i = 0; i < r; i++)
  {
    if (sumrow[i] != 0)
    {
      for (int j = 0; j < s; j++)
      {
        {
          sum1 += (double)(input[i][j] * input[i][j]) / sumrow[i];
        }
      }
    }
  }

  for (int j = 0; j < s; j++)
  {
    sumcol = 0;
    for (int i = 0; i < r; i++)
    {
      sumcol += input[i][j];
    }
    sum2 += (double)(sumcol * sumcol) / (double)N;
  }

  return (sum1 - sum2) * (double)s;
}

vector<double> preFactorial(unsigned n)
{
  vector<double> result(n + 1);
  result[0] = 1;
  for (double x = 1; x <= (double)n; x++)
  {
    result[x] = x * result[x - 1];
  }

  return result;
}

vector<double> preLogFactorial(unsigned n)
{
  vector<double> result(n + 1);
  result[0] = 0.0;
  for (unsigned x = 1; x <= n; x++)
    result[x] = log(x) + result[x - 1];
  return result;
}


// n: sample size
// r: Number of cells to distribute values
// max: Total number of cells. If max < r, max - r cells will be 0
vector<int> disTable(int n, int r, int max)
{
    int a = n / r;
    vector<int> result(max, 0);
    for (int i = 0; i < r; i++)
    {
        result[i] = a;
    }

    int leftv = n - a * r;
    if (leftv > 0)
    {
        for (int i = 0; i < leftv; i++)
        {
            result[i]++;
        }
    }
    return result;
}
