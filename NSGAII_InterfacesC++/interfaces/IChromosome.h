/*
 * IChromosome.h
 *
 *  Created on: Jun 5, 2014
 *      Author: reginaldo
 */

#ifndef ICHROMOSOME_H_
#define ICHROMOSOME_H_

class IChromosome {
private:
	int paretoFront;
	double crowdingDistance;

public:
	IChromosome() : paretoFront(0), crowdingDistance(0) {}
	virtual ~IChromosome() {}

	/**
	 * THIS METHOD MUST BE IMPLEMENTED BY THE CLIENT.
	 * @param chromosome The individual to evaluate
	 * @return True if (this) chromosome dominates argument chromosome
	 */
	virtual bool dominates(const IChromosome * const chromosome) = 0;

	/**
	 * THIS METHOD MUST BE IMPLEMENTED BY THE CLIENT.
	 * @param objective The objective that needs to be returned.
	 */
	virtual double getObjective(const int &objective) = 0;

	inline void setParetoFront(const int &paretoFront) {
		this->paretoFront = paretoFront;
	}

	inline int getParetoFront() const {
		return this->paretoFront;
	}

	inline void setCrowdingDistance(const double &crowdingDistance) {
		this->crowdingDistance = crowdingDistance;
	}

	inline int getCrowdingDistance() const {
		return this->crowdingDistance;
	}
};

#endif /* ICHROMOSOME_H_ */
