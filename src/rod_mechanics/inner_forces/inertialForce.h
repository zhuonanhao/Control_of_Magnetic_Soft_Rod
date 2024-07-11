#ifndef INERTIALFORCE_H
#define INERTIALFORCE_H

#include "eigenIncludes.h"
#include "rod_mechanics/elasticRod.h"
#include "time_steppers/timeStepper.h"

class inertialForce
{
public:
	inertialForce(elasticRod &m_rod, timeStepper &m_stepper);
	~inertialForce();
	void computeFi();
	void computeJi();

private:
	elasticRod *rod;
	timeStepper *stepper;
    			
    int ind1, ind2, mappedInd1, mappedInd2;	
    double f, jac;
};

#endif
