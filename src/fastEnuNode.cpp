#include "fastEnuNode.h"
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

fastEnuNode::fastEnuNode()
{
}

// Cs as column sum.
fastEnuNode::fastEnuNode(vector<int> Cs, int key)
{
	Csum = Cs;
	equiKey = key;
}

void fastEnuNode::show()
{
	/*
	for (size_t i = 0; i < Csum.size(); i++) {
		cout << Csum[i] << " ";
	}
	cout << " :: ";
	for (size_t i = 0; i < ChildrenIndex.size(); i++) {
		cout << ChildrenIndex[i] << " ";
	}
	cout << "::" << ub << " :: " << lb << endl;
	*/
}

void fastEnuNode::addRsum(int Rs)
{
	Csum.push_back(Rs);
}

vector<int> fastEnuNode::getCsum()
{
	return Csum;
}

int fastEnuNode::getChildrenIndex(int ind)
{
	return ChildrenIndex[ind];
}

double fastEnuNode::getLengthToChildren(int ind)
{
	return lengthToChildren[ind];
}

double fastEnuNode::getColChisqToChildren(int ind)
{
	return colChisqToChildren[ind];
}

void fastEnuNode::setUB(double val)
{
	ub = val;
}
void fastEnuNode::setLB(double val)
{
	lb = val;
}

double fastEnuNode::getLB()
{
	return lb;
}

double fastEnuNode::getUB()
{
	return ub;
}

void fastEnuNode::setequiKey(int x)
{
	equiKey = x;
}

int fastEnuNode::getequiKey()
{
	return equiKey;
}

int fastEnuNode::getSize()
{
	return (int)ChildrenIndex.size();
}

void fastEnuNode::addPastLen(double len, double chisq)
{
	size_t j = 0;

	// Get iterator of hash
	unordered_map<double, int>::const_iterator nodeItr;
	nodeItr = nodeTable.find(chisq);

	// Case if chisq not exist in hash
	if (nodeItr == nodeTable.end())
	{
		pastLen.push_back(len);
		pastChisq.push_back(chisq);
		nodeTable.insert(std::make_pair(chisq, pastChisq.size() - 1));
	}
	else
	{
		// Case can find the key, update the value
		pastLen[nodeItr->second] += len;
	}
}

double fastEnuNode::getPastLen(int ind)
{
	return pastLen[ind];
}

int fastEnuNode::getPastSize()
{
	return (int)pastLen.size();
}

double fastEnuNode::getPastChisq(int ind)
{
	return pastChisq[ind];
}

void fastEnuNode::setLengthToEnd(double val)
{
	lengthToEnd = val;
}

double fastEnuNode::getLengthToEnd()
{
	return lengthToEnd;
}

int fastEnuNode::getPastChisqSize()
{
	return (int)pastChisq.size();
}

void fastEnuNode::quicksort(int left, int right)
{

	double pivot = pastChisq[(left + right) / 2];
	int i, j;
	// double temp;
	i = left;
	j = right;

	while (i <= j)
	{
		while (pastChisq[i] < pivot)
			i++;
		while (pastChisq[j] > pivot)
			j--;

		if (i <= j)
		{
			// swap arr[i] and arr[j]
			std::swap(pastChisq[i], pastChisq[j]);
			std::swap(pastLen[i], pastLen[j]);
			i++;
			j--;
		}
	}

	if (left < j)
		quicksort(left, j);
	if (i < right)
		quicksort(i, right);
}

int fastEnuNode::bSearch(int chisq)
{
	auto it = std::lower_bound(pastChisq.begin(), pastChisq.end(), chisq);
	std::size_t index = std::distance(pastChisq.begin(), it);
	return (int)index;
}

int fastEnuNode::isChildInList(int x)
{
	size_t i = 0;
	while (i < ChildrenIndex.size() && x != ChildrenIndex[i])
		i++;
	if (i < ChildrenIndex.size())
		return (int)i;
	return -1;
}

void fastEnuNode::addLength(double x, int index)
{
	lengthToChildren[index] += x;
}

void fastEnuNode::addChildLink(int index, double len, double colchisq)
{
	ChildrenIndex.push_back(index);
	lengthToChildren.push_back(len);
	colChisqToChildren.push_back(colchisq);
}

void fastEnuNode::setColChisqToChildren(int ind, double colchisq)
{
	colChisqToChildren[ind] = colchisq;
}

void fastEnuNode::setMinPastChisq(double val)
{
	minPastChisq = val;
}

double fastEnuNode::getMinPastChisq()
{
	return minPastChisq;
}

void fastEnuNode::setMaxPastChisq(double val)
{
	maxPastChisq = val;
}

double fastEnuNode::getMaxPastChisq()
{
	return maxPastChisq;
}
