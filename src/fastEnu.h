// The fast enumeration algorithm code to obtain P-value
// Enumerate table fit the given row/col sum

#ifndef FASTEUN_H
#define FASTEUN_H

//  EFT_DQP.h
//    Exact functional test implementaiton using dynamic and
//    quadratic programming
//
//  Created by Hien Nguyen on 7/24/2018.
//
//  Revision history:
//  2019-02-25 (Hien Nguyen): created the namespace DQP.
//
//  2019-02-25 (Hien Nguyen): the file name changed from
//     EFTNetwork.h of 2.4.7 to EFT_DQP.h to distinguish from other
//     implementations of the exact functional test.
// 2021-01-06
// Add true hash table

#include "fastEnuNode.h"
#include <iostream>
#include <array>
#include <vector>
#include <time.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iterator>
#include <unordered_map>
#include "matrixTools.h"

// struct for UEFTC return
struct THREEP
{
    double leftS;
    double rightS;
    double midP;
    double leftP;
    double stat;
    double rightP;
};

// Structure save 2 double values
struct twoDouble
{
    double funchisq;
    double addTerm1;
};

struct netResult
{
    double leftS;
    double rightS;
    double pvalue;
    bool ifL;
    bool ifR;
};

namespace fastEnu
{
    // createKey: convert Rs to equivalent int
    unsigned long int createKey(vector<int> Rs, int layer, int maxCSum);

    // compute the weight between two nodes
    double colChisq(vector<int> & Rs1, vector<int> & Rs2, int sum, const vector<int> &squares, const double &ROWMARGIN);

    double colChisq(vector<int> Rs1, const int &sum, const vector<int> &squares, const double &ROWMARGIN);

    // compute the length from the current node to the end node
    double length(const vector<int> &Rs1, const int &sum, int &layer, const vector<int> &Cs, const vector<double> &factorials);

    // compute the length between two nodes
    double length(const vector<int> &Rs1, const vector<int> &Rs2, const vector<double> &factorials);

    // compute the funchisq without the fixed marginals
    double funchisqByRow(const vector<vector<int>> &observedTable, vector<int> &RSUM, const vector<int> &squares, double &ROWMARGIN);

    // Compute the funchisq designed for uniform EFT
    twoDouble funchisqForUni(const vector<vector<int>> &observedTable, int N, vector<int> &RSUM, vector<int> &uniCSum,
                             const vector<int> &squares, double &ROWMARGIN);

    double funchisqForAE(const vector<vector<int>> &observedTable,
                         vector<int> &RSUM, vector<int> &CSUM,
                         vector<int> &maxMCSUM,
                         const vector<int> &squares, double &ROWMARGIN);

    // enumerate the children nodes given the current node, the formula is given in Network Algorithm (Mehta and Patel) paper
    void createNode(fastEnuNode & node, vector<int> Cs, const vector<int> &Rs, int layer, vector<int> &currRs, int &ncols, int sum1, int sum2,
                    const vector<int> &S, const int &i, const vector<int> &squares, const vector<double> &factorials, vector<fastEnuNode> &Layer,
                    const double &ROWMARGIN, unordered_map<unsigned long int, int> &hashTable, int maxCSum, bool exact = false);

    // exact: if true, would only return 0 so that the bound would be evaluate only by DP
    double lower_bound(int layer, const vector<int> &Rsum, const vector<int> &O_colsums, const double &ROWMARGIN, bool exact = false);
    // exact: if ture, sould only return the maximum value. bound would be evaluate only by DP
    double upper_bound(int layer, const vector<int> &Rsum, const vector<int> &O_colsums, const double &ROWMARGIN, bool exact = false);

    // main program for EFT
    double EFTNetwork(const vector<vector<int>> &inputM, bool getCDF = false);

    // EFT used for MFT
    double MFTNetwork(vector<int> RowSums, vector<int> ColSums, int N, double funchisq, const vector<double> &factorials, bool ifNull = false);

    // EFT used for UEFT
    double UEFTNetwork(vector<int> RowSums, vector<int> ColSums, int N, const vector<vector<int>> &oriM,
                       const vector<double> &factorials, bool getCDF = false);

    // The subfunction for fast enumeration code
    // funchisqAdd: Funchisq statistic and additional term from input table
    // funchisq: the statistic use to build tree
    // result: The result of P-value, stat.
    // ColSums: The uniform colsum
    // RowSums: The uniform rowSum
    // squares: square values
    // S: sets of square values for accuracy
    // factorials: set of factorials values
    // getCDF: if you need to get CDF or P-values
    // ROWMARGIN and marginal: const values
    netResult subUEFTCNetwork(const twoDouble funchisqAdd, const double funchisq, int maxCSum, const vector<int> &ColSums,
                              const vector<int> &RowSums, const double &ROWMARGIN, const double marginal,
                              const vector<int> &squares, const vector<int> &S, const vector<double> &factorials, const bool getCDF);

    // Design for UEFTC with P-value approximatly among all cases.
    // RowSums: the row sum of input matrix. For UEFT, use the uniform distribution
    // ColSums: the column sum of input matrix. For UEFT, use the uniform distribution
    // N: the sample size
    // oriM: the input matrix
    // factorials: the factorial values to improve the speed
    // getCDF: Only for distribution.  Default is False, if True, the output will be 1 - CDF
    THREEP UEFTCNetwork(vector<int> RowSums, vector<int> ColSums, int N, const vector<vector<int>> &oriM,
                        const vector<double> &factorials, bool getCDF = false);
}

#endif
