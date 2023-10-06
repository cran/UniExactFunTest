// Node structure for fast enumaration algorithm
#ifndef FASTENUNODE_H
#define FASTENUNODE_H

#include <vector>
#include <string>
#include <unordered_map>

using namespace std;

class fastEnuNode {

public:
	fastEnuNode();
	fastEnuNode(vector<int> Rs, int key);

	void show();

	void addRsum (int Rs);
	vector<int> getCsum();

	int getChildrenIndex(int ind);
	int isChildInList(int ind);

	void addLength(double val, int ind);
	double getLengthToChildren(int ind);

	void setColChisqToChildren(int ind,  double colchisq);
	double getColChisqToChildren(int ind);

	void addPastLen(double len,  double chisq);
	double getPastLen(int ind);

	double getPastChisq(int ind);
	int bSearch(int chisq);
	int getPastChisqSize();

	int getPastSize();

	void setUB(double val);
	double getUB();

	void setLB(double val);
	double getLB();

	void setequiKey(int x);
	int getequiKey();

	int getSize();

	void setLengthToEnd(double val);
	double getLengthToEnd();

	void quicksort(int left, int right);

	void addChildLink(int index, double len, double colchisq);

	void setMinPastChisq(double val);
	double getMinPastChisq();

	void setMaxPastChisq(double val);
	double getMaxPastChisq();

private:
	vector<int> Csum;	// the remaining rowsums after the previous columns are enumerated
	int equiKey;		// the equivalent integer converted from Rsum to compare the nodes faster
	double lengthToEnd; // the length from this node to the end node, used when the whole branch is counted

	double ub;			// the upper bound for this node
	double lb;			// the lower bound for this node

	vector<int> ChildrenIndex;			// to store the indices of the children node
	vector<double> lengthToChildren;	// to store the lenghths to the children node
	vector<double> colChisqToChildren;	// to store the weights to the children node (weight = partial funchisq for the column enumerated)

	vector<double> pastLen;		// the cummulative lengths from the start node to this node, each entry of this list will be called lengthSoFar
	vector<double> pastChisq;	// the cummulative weights from the start node to this node, each entry of this list will be called chisqSoFar or weightSoFar


	unordered_map<double, int> nodeTable;
	// first double is chisq, second int is the index of pastLen

	double minPastChisq;
	double maxPastChisq;
};

#endif
