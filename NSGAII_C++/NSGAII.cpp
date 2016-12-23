/*
 * NSGAII.cpp
 *
 *  Created on: Jun 5, 2014
 *      Author: reginaldo
 */

#include "NSGAII.h"
#include "interfaces/IChromosome.h"
#include "Conf.h"
#include <limits>
using namespace std;

NSGAII::NSGAII() {/** Nothing to be done here **/}

NSGAII::~NSGAII() {/** Nothing to be done here **/}

void NSGAII::fastNonDominatedSort(IChromosome * const chromosomes, const int &popSize) {
	Domination S[popSize];
	vector<int> F;

	for (int p = 0; p < popSize; ++p) {
		for (int q = 0; q < popSize; ++q) {
			if(p == q) continue;
			if(chromosomes[p].dominates(&chromosomes[q])) // chromosomes+q
				S[p].dominated.push_back(q);
			else if(chromosomes[q].dominates(&chromosomes[p]))
				S[p].dominationCount++;
		}
		if (S[p].dominationCount == 0) {
			chromosomes[p].setParetoFront(1);
			F.push_back(p);
		}
	}

	/**
	 * 		***** OBSERVATION ******
	 * It needs C++11 for this kinf of loop.
	 * Otherwise, please change for old loop style:
	 * 		for(unsigned int p = 0; p < F.size(); ++p)
	 */

	int i = 1; // i is the front counter
	while(!F.empty()) {
		vector<int> Q;
		for (int p : F) {
			for (int q : S[p].dominated) {
				S[q].dominationCount--;
				if(S[q].dominationCount == 0){
					chromosomes[q].setParetoFront(i + 1);
					Q.push_back(q);
				}
			}
		}
		i++;
		F = Q;
	}
}

void NSGAII::crowdingDistances(IChromosome * const chromosomes, const int &popSize) {
	for (int i = 0; i < popSize; ++i)
		chromosomes[i].setCrowdingDistance(0);

	//TODO YOU HAVE TO SORT ALL THE ELEMENTS BY PARETO FRONT HERE : ASCENDENT ORDER

	int firstParetoFront = chromosomes[0].getParetoFront();
	int lastParetoFront = chromosomes[popSize - 1].getParetoFront();
	int firstPos = 0, size = 0;

	for (int i = firstParetoFront; i <= lastParetoFront; i++) {

		size = getParetoFrontSize(chromosomes, firstPos, i, popSize);

		if(size > 1) {
			for (int obj = 0; obj < NOBJECTIVES; obj++) {

				//TODO YOU HAVE TO SORT THE INDIVIDUAL BY OBJECTIVES HERE : ASCENDENT ORDER
				// EACH ITERATION OF THIS 'FOR' MEANS ONE DIFFERENT OBJECTIVE SORTING
				// DON'T FORGET TO SORT ONLY BY PARETO FRONT, IN OTHER WORDS, WINDOW [firstPos, size].

				chromosomes[firstPos].setCrowdingDistance(numeric_limits<double>::infinity());
				chromosomes[size - 1].setCrowdingDistance(numeric_limits<double>::infinity());
				for (int j = 1; j < size - 1; j++) {
					chromosomes[j].setCrowdingDistance(chromosomes[j].getCrowdingDistance() + chromosomes[j + 1].getObjective(obj) - chromosomes[j - 1].getObjective(obj));
				}
			}
		} else {
			chromosomes[firstPos].setCrowdingDistance(numeric_limits<double>::infinity());
		}
		firstPos += size;
	}
}

int NSGAII::getParetoFrontSize(IChromosome * const chromosomes, const int &firstPos, const int &paretoFront, const int &popSize) {
	int count = 0;
	for (int i = firstPos; i < popSize; i++) {
		if(chromosomes[i].getParetoFront() != paretoFront)
			break;
		count++;
	}
	return count;
}
