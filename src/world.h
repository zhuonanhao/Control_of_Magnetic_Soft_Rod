#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

#include <time.h>

// include elastic rod class
#include "rod_mechanics/elasticRod.h"

// include force classes
#include "rod_mechanics/inner_forces/elasticStretchingForce.h"
#include "rod_mechanics/inner_forces/elasticBendingForce.h"
#include "rod_mechanics/inner_forces/elasticTwistingForce.h"
#include "rod_mechanics/external_forces/externalGravityForce.h"
#include "rod_mechanics/control_forces/externalMagneticForce.h" 
#include "rod_mechanics/inner_forces/inertialForce.h"
#include "rod_mechanics/external_forces/externalContactForce.h"

// include external force
#include "rod_mechanics/external_forces/dampingForce.h"
// #include "dampingForce.cpp"

// include time stepper
#include "time_steppers/timeStepper.h"

// include input file and option
#include "initialization/setInput.h"

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	void setRodStepper();
	void updateTimeStep();
	int simulationRunning();
	int numPoints();
	double getScaledCoordinate(int i);
	double getCurrentTime();
	double getTotalTime();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);

	Vector3d getM1(int i);
	Vector3d getX(int i);

	double getTubeRadius();
	Vector3d getTubeNode(int i);
	int getTubeNv();
		
private:

	// Physical parameters
	double RodLength;
	double helixradius;
	double helixpitch;
	double rodRadius;
	int numVertices;
	double youngM;
	double Poisson;
	double shearM;
	double deltaTime;
	double totalTime;
	double density;
	Vector3d gVector;
	double viscosity;
	double scaleRender;
	double speed;

	double dBar;
    double stiffness;

	Vector3d baVector;
    Vector3d brVector;
    double muZero;

	// Viscous drag coefficients
    double eta_per, eta_par;
    
	double tol, stol;
	int maxIter; // maximum number of iterations
	double characteristicForce;
	double forceTol;
	
	// Geometry
	MatrixXd vertices;
	
	// Rod
	elasticRod *rod;
	
	// set up the time stepper
	timeStepper *stepper;
	double *totalForce;
	double currentTime;
	
	// declare the forces
	elasticStretchingForce *m_stretchForce;
	elasticBendingForce *m_bendingForce;
	elasticTwistingForce *m_twistingForce;
	inertialForce *m_inertialForce;	
	externalGravityForce *m_gravityForce;
	externalMagneticForce *m_magneticForce;
	dampingForce *m_dampingForce;
	externalContactForce *m_externalContactForce;
	
	int Nstep;
	int timeStep;
	int iter;

	void rodGeometry();
	void rodBoundaryCondition();
    
	bool render; // should the OpenGL rendering be included?
	bool saveData; // should data be written to a file?

	Vector3d xInitial;
};

#endif
