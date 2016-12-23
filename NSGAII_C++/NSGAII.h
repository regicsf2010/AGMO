/*
 * NSGAII.h
 *
 *  Created on: Jun 5, 2014
 *      Author: reginaldo
 */

#ifndef NSGAII_H_
#define NSGAII_H_

class IChromosome;
#include <vector>
using namespace std;

struct Domination {
	int dominationCount;
	vector<int> dominated;
};

class NSGAII {
private:
	int getParetoFrontSize(IChromosome * const chromosomes, const int &firstPos, const int &paretoFront, const int &popSize);

public:
	NSGAII();
	virtual ~NSGAII();

	/**
	 * This method can be overwritten.
	 * @param population The whole population
	 */
	virtual void fastNonDominatedSort(IChromosome * const chromosomes, const int &popSize);

	/**
	 * This method can be overwritten.
	 * @param population The whole population
	 */
	virtual void crowdingDistances(IChromosome * const chromosomes, const int &popSize);

};

#endif /* NSGAII_H_ */
