//  EFT_DQP.cpp
//    Exact functional test implementaiton using dynamic programming
//
//  Created by Hien Nguyen on 7/24/2018.
//
//  Revision history:
//  2019-09-07 (Hien Nguyen): fix the error case when the branch and bound occurs at the root node.
//
//  2019-02-25 (Hien Nguyen): insulate code into the namespace DP.
//
//  2019-02-25 (Hien Nguyen): the file name changed from
//     EFTNetwork.cpp of 2.4.7 to EFT_DQP.cpp to distinguish from other
//     implementations of the exact functional test.
//

#include "fastEnu.h"

using namespace std;

// createKey: convert Rs to equivalent string
unsigned long int fastEnu::createKey(vector<int> Rs, int layer, int maxCSum)
{
  maxCSum++;
  Rs.push_back(layer);
  sort(Rs.begin(), Rs.end());
  unsigned long int eInt = 0;
  for (size_t x = 0; x < Rs.size(); x++)
  {
    eInt *= (unsigned long int)maxCSum;
    eInt += (unsigned long int)Rs[x];
  }
  return eInt;
}

// compute the weight between two nodes
// colchisq = colChisq(network[1][x].getCsum(), sumCol[0], squares, ROWMARGIN);
double fastEnu::colChisq(vector<int> &Rs1, vector<int> &Rs2, int sum, const vector<int> &squares, const double &ROWMARGIN)
{
  if (sum > 0)
  {
    double colchisq = 0.0;
    for (size_t x = 0; x < Rs2.size(); x++)
    {
      colchisq += squares[Rs1[x] - Rs2[x]];
    }
    colchisq = colchisq * ROWMARGIN / sum;
    return (colchisq);
  }
  else
    return 0;
}

double fastEnu::colChisq(vector<int> Rs1, const int &sum, const vector<int> &squares, const double &ROWMARGIN)
{
  if (sum > 0)
  {
    double colchisq = 0.0;
    for (size_t x = 0; x < Rs1.size(); x++)
    {
      colchisq += squares[Rs1[x]];
    }
    colchisq = colchisq * ROWMARGIN / sum;
    return (colchisq);
  }
  else
    return 0;
}

// compute the length from the current node to the end node
double fastEnu::length(const vector<int> &Rs1, const int &sum, int &layer, const vector<int> &Cs, const vector<double> &factorials)
{
  double length = factorials[sum];
  for (size_t x = 0; x < Rs1.size(); x++)
  {
    length /= factorials[Rs1[x]];
  }
  for (int x = 0; x < layer; x++)
  {
    length /= factorials[Cs[x]];
  }
  return (length);
}

// compute the length between two nodes
double fastEnu::length(const vector<int> &Rs1, const vector<int> &Rs2, const vector<double> &factorials)
{
  double length = 1.0;
  for (size_t x = 0; x < Rs2.size(); x++)
  {
    length /= factorials[Rs1[x] - Rs2[x]];
  }
  return (length);
}

// compute the funchisq without the fixed marginals
double fastEnu::funchisqByRow(const vector<vector<int>> &observedTable,
                              vector<int> &RSUM, const vector<int> &squares, double &ROWMARGIN)
{
  double rowchisq = 0.0;
  double funchisq = 0.0;
  for (size_t i = 0; i < observedTable.size(); i++)
  {
    rowchisq = 0;
    if (RSUM[i] > 0)
    {
      for (size_t j = 0; j < observedTable[0].size(); j++)
      {
        rowchisq += squares[observedTable[i][j]];
      }
      rowchisq = rowchisq * ROWMARGIN / RSUM[i];
    }
    funchisq += rowchisq;
  }
  return (funchisq);
}

// compute the funchisq without the fixed marginals
twoDouble fastEnu::funchisqForUni(const vector<vector<int>> &observedTable, int N, vector<int> &RSUM, vector<int> &uniCSum,
                                  const vector<int> &squares, double &ROWMARGIN)
{
  double rowchisq = 0.0;
  double funchisq = 0.0;

  vector<int> CSUM(observedTable[0].size(), 0);
  for (unsigned i = 0; i < observedTable.size(); i++)
  {
    for (unsigned j = 0; j < CSUM.size(); j++)
    {
      CSUM[j] += observedTable[i][j];
    }
  }

  for (size_t i = 0; i < observedTable.size(); i++)
  {
    rowchisq = 0;
    if (RSUM[i] > 0)
    {
      for (size_t j = 0; j < observedTable[0].size(); j++)
      {
        rowchisq += squares[observedTable[i][j]];
      }
      rowchisq = rowchisq * ROWMARGIN / RSUM[i];
    }
    funchisq += rowchisq;
  }
  // The term using uniform column sum
  double addTerm1 = 0;
  for (unsigned i = 0; i < uniCSum.size(); i++)
  {
    addTerm1 += squares[uniCSum[i]] * ROWMARGIN / N;
  }
  // The term use observed column sum
  double addTerm2 = 0;
  for (unsigned i = 0; i < CSUM.size(); i++)
  {
    addTerm2 += squares[CSUM[i]] * ROWMARGIN / N;
  }
  funchisq = funchisq + addTerm1 - addTerm2;

  // return the additional term and funchisq.
  twoDouble result;
  result.addTerm1 = addTerm1;
  result.funchisq = funchisq;
  return (result);
}

// compute the funchisq without the fixed marginals
double fastEnu::funchisqForAE(const vector<vector<int>> &observedTable,
                              vector<int> &RSUM, vector<int> &CSUM,
                              vector<int> &maxMCSUM,
                              const vector<int> &squares, double &ROWMARGIN)
{
  double N1 = 0;
  double N2 = 0;

  for (unsigned i = 0; i < RSUM.size(); i++)
  {
    N1 += RSUM[i];
  }
  for (unsigned i = 0; i < maxMCSUM.size(); i++)
  {
    N2 += maxMCSUM[i];
  }

  if (N1 == 0 || N2 == 0)
  {
    return 0;
  }

  double rowchisq = 0.0;
  double funchisq = 0.0;
  for (size_t i = 0; i < observedTable.size(); i++)
  {
    rowchisq = 0;
    if (RSUM[i] > 0)
    {
      for (size_t j = 0; j < observedTable[0].size(); j++)
      {
        rowchisq += squares[observedTable[i][j]];
      }
      rowchisq = rowchisq * ROWMARGIN / RSUM[i];
    }
    funchisq += rowchisq;
  }

  double addTerm1, addTerm2;

  addTerm1 = 0;
  // Addition evaluation for itself
  for (size_t i = 0; i < CSUM.size(); i++)
  {
    addTerm1 += squares[CSUM[i]] * ROWMARGIN / N1;
  }

  addTerm2 = 0;
  // Addition evaluation for oriM
  for (size_t i = 0; i < maxMCSUM.size(); i++)
  {
    addTerm2 += squares[maxMCSUM[i]] * ROWMARGIN / N2;
  }

  funchisq = funchisq + (addTerm1 - addTerm2);
  return (funchisq);
}

double fastEnu::lower_bound(int layer, const vector<int> &Rsum, const vector<int> &O_colsums, const double &ROWMARGIN, bool exact)
{
  // The case for UEFTC
  if (exact)
  {
    return 0;
  }
  double lower_bound = 0;

  vector<int> U(Rsum);

  size_t nrows = Rsum.size();

  vector<size_t> order(nrows);
  for (size_t q = 0; q < nrows; ++q)
  {
    order[q] = q;
  }

  // sort U in increasing order
  sort(order.begin(), order.end(),
       [&U](size_t i1, size_t i2)
       { return U[i1] < U[i2]; });

  for (int l = 0; l < layer; l++)
  {
    // find lower bound for row l
    int runsum = 0;

    for (size_t k = 0; k < nrows; ++k)
    {
      // accumulate the lower bound
      double xavg = (O_colsums[l] - runsum) / (double)(nrows - k);
      if (U[order[k]] < xavg)
      {
        if (O_colsums[l] > 0)
          lower_bound += U[order[k]] * U[order[k]] * ROWMARGIN / (double)O_colsums[l];
        runsum += U[order[k]];
      }
      else
      {
        if (O_colsums[l] > 0)
          lower_bound += (nrows - k) * xavg * xavg * ROWMARGIN / (double)O_colsums[l];
        break;
      }
    }
  }

  return lower_bound;
}

double fastEnu::upper_bound(int layer, const vector<int> &Rsum, const vector<int> &O_colsums, const double &ROWMARGIN, bool exact)
{
  if (exact)
  {
    return (double)DBL_MAX;
  }

  double upper_bound = 0;

  vector<int> U(Rsum);

  size_t nrows = Rsum.size();

  vector<size_t> order(nrows);
  for (size_t q = 0; q < nrows; ++q)
  {
    order[q] = q;
  }

  // sort U in decreasing order
  sort(order.begin(), order.end(),
       [&U](size_t i1, size_t i2)
       { return U[i1] > U[i2]; });

  for (int l = layer - 1; l >= 0; l--)
  {
    // find lower bound for row l
    if (O_colsums[l] > 0)
    {
      int runsum = 0;

      for (size_t k = 0; k < nrows; ++k)
      {
        // accumulate the lower bound
        int xmax = O_colsums[l] - runsum;
        if (U[order[k]] < xmax)
        {
          if (O_colsums[l] > 0)
            upper_bound += U[order[k]] * U[order[k]] * ROWMARGIN / O_colsums[l];
          runsum += U[order[k]];
        }
        else if (xmax != 0)
        {
          if (O_colsums[l] > 0)
            upper_bound += xmax * xmax * ROWMARGIN / O_colsums[l];
          runsum += xmax;
        }
        else
        {
          break;
        }
      }
    }
  }

  return upper_bound;
}

// enumerate the children nodes given the current node, the formula is given in Network Algorithm (Mehta and Patel) paper
void fastEnu::createNode(fastEnuNode &node, vector<int> Cs, const vector<int> &Rs, int layer, vector<int> &currCs, int &ncols, int sum1, int sum2,
                         const vector<int> &S, const int &i, const vector<int> &squares, const vector<double> &factorials, vector<fastEnuNode> &Layer,
                         const double &ROWMARGIN, unordered_map<unsigned long int, int> &hashTable, int maxCSum, bool exact)
{
  // When the current column sum completed.
  if (i == ncols)
  {
    // Obtain the length, weight of node
    double len = fastEnu::length(Cs, currCs, factorials);
    int colchisq = fastEnu::colChisq(Cs, currCs, Rs[layer], squares, ROWMARGIN);

    unsigned long int eKey = fastEnu::createKey(currCs, layer, maxCSum);
    unordered_map<unsigned long int, int>::const_iterator keyIt;
    keyIt = hashTable.find(eKey);

    // if the child node does not exist yet, insert it to the next layer as a new node
    if (keyIt == hashTable.end())
    {
      Layer.push_back(fastEnuNode(currCs, eKey));
      node.addChildLink((int)Layer.size() - 1, len, colchisq);
      // update hashTable by insertion
      hashTable.insert(make_pair(eKey, Layer.size() - 1));

      Layer[Layer.size() - 1].setMinPastChisq(node.getMinPastChisq() + colchisq);
      Layer[Layer.size() - 1].setMaxPastChisq(node.getMaxPastChisq() + colchisq);

      Layer[Layer.size() - 1].setLB(fastEnu::lower_bound(layer, currCs, Rs, ROWMARGIN, exact));
      Layer[Layer.size() - 1].setUB(fastEnu::upper_bound(layer, currCs, Rs, ROWMARGIN, exact));

      Layer[Layer.size() - 1].setLengthToEnd(length(currCs, S[layer - 1], layer, Rs, factorials));
    }
    // Changed to fix 0 rowSum problem
    else if (Layer.size() == 0 && keyIt != hashTable.end())
    {
    }
    else
    {
      // if the child node already exists, add a new link from the current node to that child node
      int index = keyIt->second;
      node.addChildLink(index, len, colchisq);

      Layer[index].setMinPastChisq(std::min(Layer[index].getMinPastChisq(), node.getMinPastChisq() + colchisq));
      Layer[index].setMaxPastChisq(std::max(Layer[index].getMaxPastChisq(), node.getMaxPastChisq() + colchisq));
    }
  } // end i == ncols
  else
  {
    // When current column sum not completed, get the next column sum
    int lowerbound, upperbound;

    sum1 += (i > 0 ? Cs[i - 1] : 0);
    sum2 += (i > 0 ? currCs[i - 1] : 0);
    // Get the bound for the value of node.
    lowerbound = std::max(0, Cs[i] - Rs[layer] + sum1 - sum2);
    upperbound = std::min(Cs[i], (layer - 1 >= 0 ? S[layer - 1] : 0) - sum2);

    for (int x = lowerbound; x <= upperbound; x++)
    {
      currCs[i] = x; // Get currecnt column sum at this layer.
      fastEnu::createNode(node, Cs, Rs, layer, currCs, ncols,
                          sum1, sum2, S, i + 1, squares, factorials, Layer, ROWMARGIN, hashTable, maxCSum, exact);
    }
  }
}

double fastEnu::EFTNetwork(const vector<vector<int>> &inputM, bool getCDF)
{
  vector<vector<int>> observedTable = trimTable(inputM);
  // vector<vector<int>> observedTable = inputM;
  int nrows = (int)observedTable.size();
  int ncols = nrows > 0 ? (int)observedTable[0].size() : 0;

  if (nrows == 0 || ncols == 0)
  {
    return 1.0;
  }

  int N = 0;
  vector<int> RowSums(nrows, 0);
  vector<int> ColSums(ncols, 0);

  for (int i = 0; i < nrows; i++)
  {
    for (int j = 0; j < ncols; j++)
    {
      N += observedTable[i][j];
      RowSums[i] += observedTable[i][j];
      ColSums[j] += observedTable[i][j];
    }
  }

  if (N == 0)
  {
    return 1;
  }

  int maxCSum = nrows;
  for (int j = 0; j < ncols; j++)
  {
    if (maxCSum < ColSums[j])
    {
      maxCSum = ColSums[j];
    }
  }

  vector<int> squares(N + 1);
  for (int x = 0; x < N + 1; x++)
    squares[x] = x * x;

  vector<double> factorials(N + 1);
  factorials[0] = 1.0;
  for (int x = 1; x <= N; x++)
    factorials[x] = x * factorials[x - 1];
  // for (int x = 1; x <= N; x++)  factorials[x] = factorial<double>(x);

  double marginal = factorials[N];
  for (int x = 0; x < nrows; x++)
  {
    marginal /= factorials[RowSums[x]];
  }
  for (int x = 0; x < ncols; x++)
  {
    marginal /= factorials[ColSums[x]];
  }

  std::vector<int> S(nrows);
  S[0] = RowSums[0];
  for (int x = 1; x < nrows; x++)
  {
    S[x] = S[x - 1] + RowSums[x];
  }

  double ROWMARGIN = 1;
  for (int i = 0; i < nrows; i++)
  {
    if (RowSums[i] > 0)
      ROWMARGIN *= RowSums[i];
  }

  // compute the adjusted funchisq
  double funchisq;
  funchisq = fastEnu::funchisqByRow(observedTable, RowSums, squares, ROWMARGIN);

  double funchisqRight = 0;
  for (unsigned i = 0; i < ColSums.size(); i++)
  {
    funchisqRight += squares[ColSums[i]] * ROWMARGIN / N;
  }

  if (funchisq - funchisqRight == 0)
  {
    return 1;
  }

  vector<vector<fastEnuNode>> network(nrows + 1);

  // create the sink node
  vector<int> currCs(ncols, 0);

  network[nrows].push_back(fastEnuNode(ColSums, 0));

  network[nrows][0].addPastLen(1.0, 0);
  network[nrows][0].setMaxPastChisq(0);
  network[nrows][0].setMinPastChisq(0);

  network[nrows][0].setLB(fastEnu::lower_bound(nrows, ColSums, RowSums, ROWMARGIN));
  network[nrows][0].setUB(fastEnu::upper_bound(nrows, ColSums, RowSums, ROWMARGIN));
  network[nrows][0].setLengthToEnd(fastEnu::length(ColSums, S[nrows - 1], nrows, RowSums, factorials));

  if (getCDF)
  {
    if (network[nrows][0].getLB() > funchisq)
    {
      return 1;
    }
  }
  else
  {
    if (network[nrows][0].getLB() >= funchisq)
    {
      return 1;
    }
  }

  // hash table
  unordered_map<unsigned long int, int> hashTable;

  // generate the network

  for (int layer = nrows; layer > 1; layer--)
  {
    for (size_t n = 0; n < network[layer].size(); n++)
    {
      if (getCDF)
      {
        if (network[layer][n].getLB() + network[layer][n].getMinPastChisq() <= funchisq)
        {
          if (network[layer][n].getUB() + network[layer][n].getMaxPastChisq() > funchisq)
          {
            fastEnu::createNode(network[layer][n], network[layer][n].getCsum(), RowSums, layer - 1, currCs, ncols, 0, 0, S, 0, squares, factorials, network[layer - 1], ROWMARGIN, hashTable, maxCSum);
          }
        }
      }
      else
      {
        if (network[layer][n].getLB() + network[layer][n].getMinPastChisq() < funchisq)
        {
          if (network[layer][n].getUB() + network[layer][n].getMaxPastChisq() >= funchisq)
          {
            fastEnu::createNode(network[layer][n], network[layer][n].getCsum(), RowSums, layer - 1, currCs, ncols, 0, 0, S, 0, squares, factorials, network[layer - 1], ROWMARGIN, hashTable, maxCSum);
          }
        }
      }
    }
  }

  double rowchisq;
  for (size_t x = 0; x < network[1].size(); x++)
  {
    rowchisq = colChisq(network[1][x].getCsum(), RowSums[0], squares, ROWMARGIN);
    network[1][x].setLB(rowchisq);
    network[1][x].setUB(rowchisq);
  }

  // compute upperbound and lowerbound for higher layer
  double minLB = 0;
  double maxUB = 0;
  double tempBound;

  for (int layer = 2; layer <= nrows; layer++)
  {
    for (size_t node = 0; node < network[layer].size(); node++)
    {

      // network[layer][node].setLengthToEnd(length(network[layer][node].getCsum(), S[layer - 1], layer, RowSums, factorials));

      minLB = DBL_MAX;
      maxUB = 0;

      if (network[layer][node].getSize() > 0)
      {

        for (int child = 0; child < network[layer][node].getSize(); child++)
        {

          rowchisq = network[layer][node].getColChisqToChildren(child);

          tempBound = rowchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getLB();
          if (minLB > tempBound)
            minLB = tempBound;

          tempBound = rowchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getUB();
          if (maxUB < tempBound)
            maxUB = tempBound;
        }

        network[layer][node].setLB(minLB);
        network[layer][node].setUB(maxUB);
      }
    }
  }

  double lengthSoFar = 1;
  double chisqSoFar = 0;
  double pvalue = 0;
  int pastSize = 0;

  // traverse the network in a breadth-first strategy
  for (int layer = nrows; layer >= 1; layer--)
  {
    for (size_t node = 0; node < network[layer].size(); node++)
    {
      pastSize = network[layer][node].getPastSize();

      for (int i = 0; i < pastSize; i++)
      {

        chisqSoFar = network[layer][node].getPastChisq(i);

        if (((!getCDF) && chisqSoFar + network[layer][node].getUB() < funchisq) ||
            (getCDF && chisqSoFar + network[layer][node].getUB() <= funchisq))
        {
        }
        else
        {

          lengthSoFar = network[layer][node].getPastLen(i);

          if ((!getCDF) && (chisqSoFar + network[layer][node].getLB() >= funchisq))
          {
            pvalue += lengthSoFar * network[layer][node].getLengthToEnd();
          }
          else if ((getCDF) && (chisqSoFar + network[layer][node].getLB() > funchisq))
          {
            pvalue += lengthSoFar * network[layer][node].getLengthToEnd();
          }
          else
          {
            for (int k = 0; k < network[layer][node].getSize(); k++)
            {
              network[layer - 1][network[layer][node].getChildrenIndex(k)].addPastLen(
                  network[layer][node].getLengthToChildren(k) * lengthSoFar,
                  network[layer][node].getColChisqToChildren(k) + chisqSoFar);
            } // end child
          }
        }
      } // end pastChisq

    } // end node
  }   // end layer
  hashTable.clear();
  pvalue /= marginal;

  return pvalue;
}

// funchisq is the statistic obtained as input
// RowSums, ColSums: The uniform marginal sum
// N: Sample size.
// oriM: the original matrix used to evaluat the Funchisq
// factorials: Const value to speed up.
// getCDF: True if compute CDF, not P-value. Default is False
double fastEnu::UEFTNetwork(vector<int> RowSums, vector<int> ColSums, int N, const vector<vector<int>> &oriM,
                            const vector<double> &factorials, bool getCDF)
{
  // To reverse the Row and Col
  int ncols = ColSums.size();
  int nrows = RowSums.size();
  if (oriM.size() == 0)
  {
    return 1.0;
  }

  if (ncols == 0 || nrows == 0)
  {
    return 1.0;
  }

  vector<int> squares(N + 1);
  for (int x = 0; x < N + 1; x++)
    squares[x] = x * x;

  double marginal = factorials[N];
  int maxCSum = nrows;
  for (int x = 0; x < ncols; x++)
  {
    marginal /= factorials[ColSums[x]];
    if (maxCSum < ColSums[x])
    {
      maxCSum = ColSums[x];
    }
  }
  for (int x = 0; x < nrows; x++)
  {
    marginal /= factorials[RowSums[x]];
  }

  std::vector<int> S(nrows);
  S[0] = RowSums[0];
  for (int x = 1; x < nrows; x++)
  {
    S[x] = S[x - 1] + RowSums[x];
  }

  double ROWMARGIN = 1;
  for (int i = 0; i < nrows; i++)
  {
    if (RowSums[i] > 0)
      ROWMARGIN *= RowSums[i];
  }

  vector<int> oriRowSums(oriM.size(), 0);

  for (unsigned i = 0; i < oriM.size(); i++)
  {
    for (unsigned j = 0; j < oriM[0].size(); j++)
    {
      oriRowSums[i] += oriM[i][j];
    }
  }

  twoDouble funchisqAdd = fastEnu::funchisqForUni(oriM, N, oriRowSums, ColSums, squares, ROWMARGIN);

  double funchisq = funchisqAdd.funchisq;

  if ((funchisq - funchisqAdd.addTerm1) / ROWMARGIN * nrows <= 0)
  {
    return 1.0;
  }

  // Start to build the network.
  vector<vector<fastEnuNode>> network(nrows + 1);

  // create the sink node
  network[nrows].push_back(fastEnuNode(ColSums, 0));
  network[nrows][0].addPastLen(1.0, 0);
  network[nrows][0].setMaxPastChisq(0);
  network[nrows][0].setMinPastChisq(0);

  network[nrows][0].setLB(fastEnu::lower_bound(nrows, ColSums, RowSums, ROWMARGIN));
  network[nrows][0].setUB(fastEnu::upper_bound(nrows, ColSums, RowSums, ROWMARGIN));
  network[nrows][0].setLengthToEnd(fastEnu::length(ColSums, S[nrows - 1], nrows, RowSums, factorials));

  if (getCDF)
  {
    if (network[nrows][0].getLB() > funchisq)
    {
      return 1;
    }
  }
  else
  {
    if (network[nrows][0].getLB() >= funchisq)
    {
      return 1;
    }
  }

  // hash table
  unordered_map<unsigned long int, int> hashTable;

  // generate the network
  vector<int> currCs(ncols);
  for (int layer = nrows; layer > 1; layer--)
  {
    for (size_t n = 0; n < network[layer].size(); n++)
    {
      if (getCDF)
      {
        if (network[layer][n].getLB() + network[layer][n].getMinPastChisq() <= funchisq)
        {
          if (network[layer][n].getUB() + network[layer][n].getMaxPastChisq() > funchisq)
          {
            fastEnu::createNode(network[layer][n], network[layer][n].getCsum(), RowSums, layer - 1, currCs, ncols, 0, 0, S, 0, squares, factorials, network[layer - 1], ROWMARGIN, hashTable, maxCSum);
          }
        }
      }
      else
      {
        if (network[layer][n].getLB() + network[layer][n].getMinPastChisq() < funchisq)
        {
          if (network[layer][n].getUB() + network[layer][n].getMaxPastChisq() >= funchisq)
          {
            fastEnu::createNode(network[layer][n], network[layer][n].getCsum(), RowSums, layer - 1, currCs, ncols, 0, 0, S, 0, squares, factorials, network[layer - 1], ROWMARGIN, hashTable, maxCSum);
          }
        }
      }
    }
  }

  double colchisq;
  for (size_t x = 0; x < network[1].size(); x++)
  {
    colchisq = colChisq(network[1][x].getCsum(), RowSums[0], squares, ROWMARGIN);
    network[1][x].setLB(colchisq);
    network[1][x].setUB(colchisq);
  }

  // compute upperbound and lowerbound for higher layer
  double minLB = 0;
  double maxUB = 0;
  double tempBound;

  for (int layer = 2; layer <= nrows; layer++)
  {
    for (size_t node = 0; node < network[layer].size(); node++)
    {

      // network[layer][node].setLengthToEnd(length(network[layer][node].getCsum(), S[layer - 1], layer, ColSums, factorials));

      minLB = DBL_MAX;
      maxUB = 0;

      if (network[layer][node].getSize() > 0)
      {

        for (int child = 0; child < network[layer][node].getSize(); child++)
        {

          colchisq = network[layer][node].getColChisqToChildren(child);

          tempBound = colchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getLB();
          if (minLB > tempBound)
            minLB = tempBound;

          tempBound = colchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getUB();
          if (maxUB < tempBound)
            maxUB = tempBound;
        }

        network[layer][node].setLB(minLB);
        network[layer][node].setUB(maxUB);
      }
    }
  }

  double lengthSoFar = 1;
  double chisqSoFar = 0;
  double pvalue = 0;
  int pastSize = 0;

  // traverse the network in a breadth-first strategy
  for (int layer = nrows; layer >= 1; layer--)
  {
    for (size_t node = 0; node < network[layer].size(); node++)
    {
      pastSize = network[layer][node].getPastSize();

      for (int i = 0; i < pastSize; i++)
      {

        chisqSoFar = network[layer][node].getPastChisq(i);

        if (
            ((!getCDF) && chisqSoFar + network[layer][node].getUB() < funchisq) ||
            (getCDF && chisqSoFar + network[layer][node].getUB() <= funchisq))
        {
        }
        else
        {

          lengthSoFar = network[layer][node].getPastLen(i);

          if ((!getCDF) && (chisqSoFar + network[layer][node].getLB() >= funchisq))
          {
            pvalue += lengthSoFar * network[layer][node].getLengthToEnd();
          }
          else if ((getCDF) && (chisqSoFar + network[layer][node].getLB() > funchisq))
          {
            pvalue += lengthSoFar * network[layer][node].getLengthToEnd();
          }
          else
          {
            for (int k = 0; k < network[layer][node].getSize(); k++)
            {
              network[layer - 1][network[layer][node].getChildrenIndex(k)].addPastLen(
                  network[layer][node].getLengthToChildren(k) * lengthSoFar,
                  network[layer][node].getColChisqToChildren(k) + chisqSoFar);
            } // end child
          }
        }
      } // end pastChisq

    } // end node
  }   // end layer

  hashTable.clear();
  pvalue /= marginal;

  return pvalue;
}

// The sub-function for UEFTCNetwork. Should used within UEFTCNetwork
netResult fastEnu::subUEFTCNetwork(const twoDouble funchisqAdd, const double funchisq,
                                   int maxCSum, const vector<int> &ColSums, const vector<int> &RowSums,
                                   const double &ROWMARGIN, const double marginal, const vector<int> &squares,
                                   const vector<int> &S, const vector<double> &factorials, const bool getCDF)
{
  netResult result;
  result.pvalue = 1;
  result.leftS = 0;
  result.rightS = DBL_MAX;
  result.ifL = false;
  result.ifR = false;

  // Set const values.
  int nrows = RowSums.size();
  int ncols = ColSums.size();
  // Initial network
  vector<vector<fastEnuNode>> network(nrows + 1);

  // create the start node
  network[nrows].push_back(fastEnuNode(ColSums, 0));
  network[nrows][0].addPastLen(1.0, 0);
  network[nrows][0].setMaxPastChisq(0);
  network[nrows][0].setMinPastChisq(0);

  network[nrows][0].setLB(fastEnu::lower_bound(nrows, ColSums, RowSums, ROWMARGIN, true));
  network[nrows][0].setUB(fastEnu::upper_bound(nrows, ColSums, RowSums, ROWMARGIN, true));
  network[nrows][0].setLengthToEnd(fastEnu::length(ColSums, S[nrows - 1], nrows, RowSums, factorials));

  if (getCDF)
  {
    if (network[nrows][0].getLB() > funchisq)
    {
      return result;
    }
  }
  else
  {
    if (network[nrows][0].getLB() >= funchisq)
    {
      return result;
    }
  }

  // hash table
  unordered_map<unsigned long int, int> hashTable;

  // generate the network
  vector<int> currCs(ncols, 0);
  for (int layer = nrows; layer > 1; layer--)
  {
    for (size_t n = 0; n < network[layer].size(); n++)
    {
      if (getCDF)
      {
        if (network[layer][n].getLB() + network[layer][n].getMinPastChisq() <= funchisq)
        {
          if (network[layer][n].getUB() + network[layer][n].getMaxPastChisq() > funchisq)
          {
            fastEnu::createNode(network[layer][n], network[layer][n].getCsum(), RowSums, layer - 1, currCs, ncols, 0, 0, S, 0, squares, factorials, network[layer - 1], ROWMARGIN, hashTable, maxCSum, true);
          }
        }
      }
      else
      {
        if (network[layer][n].getLB() + network[layer][n].getMinPastChisq() < funchisq)
        {
          if (network[layer][n].getUB() + network[layer][n].getMaxPastChisq() >= funchisq)
          {
            fastEnu::createNode(network[layer][n], network[layer][n].getCsum(), RowSums, layer - 1, currCs, ncols, 0, 0, S, 0, squares, factorials, network[layer - 1], ROWMARGIN, hashTable, maxCSum, true);
          }
        }

      } // end getCDF
    }   // end n loop
  }     // end layer loop

  double colchisq;
  for (size_t x = 0; x < network[1].size(); x++)
  {
    colchisq = colChisq(network[1][x].getCsum(), RowSums[0], squares, ROWMARGIN);
    network[1][x].setLB(colchisq);
    network[1][x].setUB(colchisq);
  }

  // compute upperbound and lowerbound for higher layer
  double minLB = 0;
  double maxUB = 0;
  double tempBound;

  for (int layer = 2; layer <= nrows; layer++)
  {
    for (size_t node = 0; node < network[layer].size(); node++)
    {

      // network[layer][node].setLengthToEnd(length(network[layer][node].getCsum(), S[layer - 1], layer, ColSums, factorials));

      minLB = DBL_MAX;
      maxUB = 0;

      if (network[layer][node].getSize() > 0)
      {

        for (int child = 0; child < network[layer][node].getSize(); child++)
        {

          colchisq = network[layer][node].getColChisqToChildren(child);

          tempBound = colchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getLB();
          if (minLB > tempBound)
            minLB = tempBound;

          tempBound = colchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getUB();
          if (maxUB < tempBound)
            maxUB = tempBound;
        }

        network[layer][node].setLB(minLB);
        network[layer][node].setUB(maxUB);
      }
    }
  }

  // Set initial values to get pvalue
  double lengthSoFar = 1;
  double chisqSoFar = 0;
  int pastSize = 0;
  double UB, LB;
  result.pvalue = 0;
  // traverse the network in a breadth-first strategy
  for (int layer = nrows; layer >= 1; layer--)
  {
    for (size_t node = 0; node < network[layer].size(); node++)
    {
      pastSize = network[layer][node].getPastSize();

      for (int i = 0; i < pastSize; i++)
      {

        chisqSoFar = network[layer][node].getPastChisq(i);

        UB = chisqSoFar + network[layer][node].getUB();
        // Update left S;
        if ((UB <= funchisq) && (UB > result.leftS))
        {
          result.leftS = UB;
          result.ifL = true;
        }

        if (((!getCDF) && UB < funchisq) || (getCDF && UB <= funchisq))
        {
        }
        else
        {
          lengthSoFar = network[layer][node].getPastLen(i);
          LB = chisqSoFar + network[layer][node].getLB();

          if (getCDF)
          {
            // Update rightS
            if (
                LB >= funchisq &&
                ((LB < result.rightS) || (!result.ifR)))
            {
              result.rightS = LB;
              result.ifR = true;
            }

            if (LB > funchisq)
            {
              result.pvalue += lengthSoFar * network[layer][node].getLengthToEnd();
            }
            else
            {
              for (int k = 0; k < network[layer][node].getSize(); k++)
              {
                network[layer - 1][network[layer][node].getChildrenIndex(k)].addPastLen(
                    network[layer][node].getLengthToChildren(k) * lengthSoFar,
                    network[layer][node].getColChisqToChildren(k) + chisqSoFar);
              } // end child
            }
          } // end getCDF
          else
          {
            if (LB >= funchisq)
            {
              result.pvalue += lengthSoFar * network[layer][node].getLengthToEnd();

              if (LB < result.rightS)
              {
                result.rightS = LB;
                result.ifR = true;
              }
            }
            else
            {
              for (int k = 0; k < network[layer][node].getSize(); k++)
              {
                network[layer - 1][network[layer][node].getChildrenIndex(k)].addPastLen(
                    network[layer][node].getLengthToChildren(k) * lengthSoFar,
                    network[layer][node].getColChisqToChildren(k) + chisqSoFar);
              } // end child
            }
          }
        }
      } // end pastChisq

    } // end node
  }   // end layer

  hashTable.clear();
  result.pvalue /= marginal;
  return result;
}

// Design for UEFTC with P-value approximatly among all cases.
// RowSums: the row sum of input matrix. For UEFT, use the uniform distribution
// ColSums: the column sum of input matrix. For UEFT, use the uniform distribution
// N: the sample size
// oriM: the input matrix
// factorials: the factorial values to improve the speed
// getCDF: Only for distribution.  Default is False, if True, the output will be 1 - CDF
THREEP fastEnu::UEFTCNetwork(vector<int> RowSums, vector<int> ColSums, int N, const vector<vector<int>> &oriM,
                             const vector<double> &factorials, bool getCDF)
{
  // Initial the stat and P-value
  THREEP result; // The result with stat and P
  result.rightS = 0;
  result.rightP = 1;
  result.leftS = 0;
  result.leftP = 1;
  result.midP = 1;
  result.stat = 0;
  bool ifR = false; // false if the maximum funchisq reach
  bool ifL = false; // false if the minimum funchisq reach

  int ncols = ColSums.size(); // row of table
  int nrows = RowSums.size(); // col of table

  // Case for empty input
  if (oriM.size() == 0)
  {
    return result;
  }
  if (ncols == 0 || nrows == 0)
  {
    return result;
  }

  // Get squares to speed up
  vector<int> squares(N + 1);
  for (int x = 0; x < N + 1; x++)
    squares[x] = x * x;

  // marginal: for accuracy
  double marginal = factorials[N];
  // maxCSum: for hashtable.
  int maxCSum = nrows;

  for (int x = 0; x < ncols; x++)
  {
    marginal /= factorials[ColSums[x]];
    if (maxCSum < ColSums[x])
    {
      maxCSum = ColSums[x];
    }
  }
  for (int x = 0; x < nrows; x++)
  {
    marginal /= factorials[RowSums[x]];
  }
  // S: for accuracy
  std::vector<int> S(nrows);
  S[0] = RowSums[0];
  for (int x = 1; x < nrows; x++)
  {
    S[x] = S[x - 1] + RowSums[x];
  }
  // ROWMARGIN: for accuracy
  double ROWMARGIN = 1;
  for (int i = 0; i < nrows; i++)
  {
    if (RowSums[i] > 0)
      ROWMARGIN *= RowSums[i];
  }

  // The the observed rowSum, different from the uniform row sum
  vector<int> oriRowSums(oriM.size(), 0);
  for (unsigned i = 0; i < oriM.size(); i++)
  {
    for (unsigned j = 0; j < oriM[0].size(); j++)
    {
      oriRowSums[i] += oriM[i][j];
    }
  }

  // Obtain the funchisq statistic and an addidional term
  twoDouble funchisqAdd = fastEnu::funchisqForUni(oriM, N, oriRowSums, ColSums, squares, ROWMARGIN);
  // The temporary stat used to get Pvalue.
  double funchisq = funchisqAdd.funchisq;
  // In case of 0 statistic
  if ((funchisq - funchisqAdd.addTerm1) <= 0)
  {
    return result;
  }

  // first pValue evaluation
  result.stat = (funchisq - funchisqAdd.addTerm1) / ROWMARGIN * ncols;
  // partial result of pvalue and statistic for each enumeration
  netResult partResult;
  partResult = subUEFTCNetwork(funchisqAdd, funchisq, maxCSum, ColSums, RowSums, ROWMARGIN, marginal,
                               squares, S, factorials, getCDF);

  // Update the result
  result.midP = partResult.pvalue;
  result.leftS = partResult.leftS;
  result.rightS = partResult.rightS;
  ifR = partResult.ifR;
  ifL = partResult.ifL;

  // Get left values. Only for P-value
  if (!getCDF)
  {
    funchisq = result.leftS;

    // Case1, normal case left stat exist and stat not zero.
    if ((funchisq - funchisqAdd.addTerm1) > 0 && ifL)
    {
      partResult = subUEFTCNetwork(funchisqAdd, funchisq, maxCSum, ColSums, RowSums, ROWMARGIN, marginal,
                                   squares, S, factorials, getCDF);
      // Update the result
      result.leftP = partResult.pvalue;
      result.leftS = (result.leftS - funchisqAdd.addTerm1) / ROWMARGIN * ncols;
    }
    // Case2: the left stat exist and the stat is 0
    else if (ifL)
    {
      result.leftS = 0;
    }
    // Case3: left stat not exist.
  }

  // Case1: if right stat never updated. Means the input stat is the maximum
  if (!ifR)
  {
    // Computer the maximum of funchisq
    // Right stat be the maximum
    double k = ncols < nrows ? (double)ncols : (double)nrows;
    result.rightS = ((double)N - (double)N / k) * (double)ncols;
    result.rightP = 0;
  }
  // Case2: strange case. if right stat is 0
  else if (result.rightS - funchisqAdd.addTerm1 <= 0)
  {
    result.rightP = result.midP;
    result.rightS = result.stat;
  }
  // Normal case3
  else
  {
    funchisq = result.rightS;
    // For right , use the CDF
    getCDF = true;
    partResult = subUEFTCNetwork(funchisqAdd, funchisq, maxCSum, ColSums, RowSums, ROWMARGIN, marginal,
                                 squares, S, factorials, getCDF);
    result.rightP = partResult.pvalue;
    result.rightS = (result.rightS - funchisqAdd.addTerm1) / ROWMARGIN * ncols;
  }

  return result;
}

// funchisq is the statistic obtained as input
double fastEnu::MFTNetwork(vector<int> RowSums, vector<int> ColSums, int N, double funchisq, const vector<double> &factorials, bool ifNull)
{
  double tempChisq;
  // To reverse the Row and Col
  int ncols = ColSums.size();
  int nrows = RowSums.size();

  if (ncols == 0 || nrows == 0)
  {
    return 1.0;
  }

  vector<int> squares(N + 1);
  for (int x = 0; x < N + 1; x++)
    squares[x] = x * x;

  double marginal = factorials[N];
  int maxCSum = nrows;
  for (int x = 0; x < ncols; x++)
  {
    marginal /= factorials[ColSums[x]];
    if (maxCSum < ColSums[x])
    {
      maxCSum = ColSums[x];
    }
  }
  for (int x = 0; x < nrows; x++)
  {
    marginal /= factorials[RowSums[x]];
  }

  std::vector<int> S(nrows);
  S[0] = RowSums[0];
  for (int x = 1; x < nrows; x++)
  {
    S[x] = S[x - 1] + RowSums[x];
  }

  double ROWMARGIN = 1;
  for (int i = 0; i < nrows; i++)
  {
    if (RowSums[i] > 0)
      ROWMARGIN *= RowSums[i];
  }

  // CHANGE
  // compute the adjusted funchisq
  if (ifNull)
  {
    funchisq = funchisq / ncols;

    for (unsigned i = 0; i < ncols; i++)
    {
      funchisq += (double)(ColSums[i] * ColSums[i]) / (double)N;
    }
  }

  funchisq = funchisq * ROWMARGIN;

  if (funchisq == 0)
  {
    return 1.0;
  }

  vector<vector<fastEnuNode>> network(nrows + 1);

  // create the start node
  network[nrows].push_back(fastEnuNode(ColSums, 0));
  network[nrows][0].addPastLen(1.0, 0);
  network[nrows][0].setMaxPastChisq(0);
  network[nrows][0].setMinPastChisq(0);

  network[nrows][0].setLB(fastEnu::lower_bound(nrows, ColSums, RowSums, ROWMARGIN));
  network[nrows][0].setUB(fastEnu::upper_bound(nrows, ColSums, RowSums, ROWMARGIN));
  network[nrows][0].setLengthToEnd(fastEnu::length(ColSums, S[nrows - 1], nrows, RowSums, factorials));

  tempChisq = network[nrows][0].getLB();
  if (tempChisq - funchisq > -0.0000001)
  {
    return 1;
  }

  // hash table
  unordered_map<unsigned long int, int> hashTable;

  // generate the network
  vector<int> currCs(ncols);
  for (int layer = nrows; layer > 1; layer--)
  {
    for (size_t n = 0; n < network[layer].size(); n++)
    {
      tempChisq = network[layer][n].getLB() + network[layer][n].getMinPastChisq();
      if (tempChisq < funchisq)
      {

        tempChisq = network[layer][n].getUB() + network[layer][n].getMaxPastChisq();
        if (tempChisq - funchisq > -0.0000001)
        {
          fastEnu::createNode(network[layer][n], network[layer][n].getCsum(), RowSums, layer - 1, currCs, ncols, 0, 0, S, 0, squares, factorials, network[layer - 1], ROWMARGIN, hashTable, maxCSum);
        }
      }
    }
  }

  double colchisq;
  for (size_t x = 0; x < network[1].size(); x++)
  {
    colchisq = colChisq(network[1][x].getCsum(), RowSums[0], squares, ROWMARGIN);
    network[1][x].setLB(colchisq);
    network[1][x].setUB(colchisq);
  }

  // compute upperbound and lowerbound for higher layer
  double minLB = 0;
  double maxUB = 0;
  double tempBound;

  for (int layer = 2; layer <= nrows; layer++)
  {
    for (size_t node = 0; node < network[layer].size(); node++)
    {

      // network[layer][node].setLengthToEnd(length(network[layer][node].getCsum(), S[layer - 1], layer, ColSums, factorials));

      minLB = DBL_MAX;
      maxUB = 0;

      if (network[layer][node].getSize() > 0)
      {

        for (int child = 0; child < network[layer][node].getSize(); child++)
        {

          colchisq = network[layer][node].getColChisqToChildren(child);

          tempBound = colchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getLB();
          if (minLB > tempBound)
            minLB = tempBound;

          tempBound = colchisq + network[layer - 1][network[layer][node].getChildrenIndex(child)].getUB();
          if (maxUB < tempBound)
            maxUB = tempBound;
        }

        network[layer][node].setLB(minLB);
        network[layer][node].setUB(maxUB);
      }
    }
  }

  double lengthSoFar = 1;
  double chisqSoFar = 0;
  double pvalue = 0;
  int pastSize = 0;

  // traverse the network in a breadth-first strategy
  for (int layer = nrows; layer >= 1; layer--)
  {
    for (size_t node = 0; node < network[layer].size(); node++)
    {
      pastSize = network[layer][node].getPastSize();

      for (int i = 0; i < pastSize; i++)
      {

        chisqSoFar = network[layer][node].getPastChisq(i);

        tempChisq = chisqSoFar + network[layer][node].getUB();
        if (tempChisq - funchisq < -0.0000001)
        {
        }
        else
        {

          lengthSoFar = network[layer][node].getPastLen(i);

          tempChisq = chisqSoFar + network[layer][node].getLB();
          // Change -0.0000001 to positive might keep correct
          if (tempChisq - funchisq > -0.0000001)
          {
            pvalue += lengthSoFar * network[layer][node].getLengthToEnd();
          }
          else
          {
            for (int k = 0; k < network[layer][node].getSize(); k++)
            {
              network[layer - 1][network[layer][node].getChildrenIndex(k)].addPastLen(
                  network[layer][node].getLengthToChildren(k) * lengthSoFar,
                  network[layer][node].getColChisqToChildren(k) + chisqSoFar);
            } // end child
          }
        }
      } // end pastChisq

    } // end node
  }   // end layer

  hashTable.clear();
  pvalue /= marginal;

  return pvalue;
}
