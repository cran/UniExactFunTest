// Sets of tools used for sets of test: MFT, UEFT and others

#ifndef MATRIXTOOLS_H
#define MATRIXTOOLS_H

#include <numeric>
#include <cmath>
#include <vector>

using namespace std;

////////////////
// General code:

// trimTable Function. Delete the replicated 0 at same row/col
//  Created by Joe Song on 9/22/19.
//  Copyright © 2019 Joe Song. All rights reserved.
vector<vector<int>> trimTable(const vector<vector<int>> &table);


// Get Funchisq stat for X->Y
double getchisqInteger(const vector<vector<int>> &);

// Code to compute the Factorial value
vector<double> preFactorial(unsigned n);
// Code to compute the logFactorial value, seprate file
vector<double> preLogFactorial(unsigned n);


// Function to distribution n values to r cells
vector<int> disTable(int n, int r, int max);

#endif /* trimTable_h */
