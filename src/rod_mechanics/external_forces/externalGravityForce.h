#ifndef EXTERNALGRAVITYFORCE_H
#define EXTERNALGRAVITYFORCE_H

#include "eigenIncludes.h"
#include "rod_mechanics/elasticRod.h"
#include "time_steppers/timeStepper.h"

class externalGravityForce
{
public:
	externalGravityForce(elasticRod &m_rod, timeStepper &m_stepper, Vector3d m_gVector);
	~externalGravityForce();
	void computeFg();
	void computeJg();
	
private:
	elasticRod *rod;
	timeStepper *stepper;
    Vector3d gVector;
    VectorXd massGravity;
    void setGravity();
};

#endif
