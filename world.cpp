#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean
	
	// Physical parameters
	RodLength = m_inputData.GetScalarOpt("RodLength");      // meter
    helixradius = m_inputData.GetScalarOpt("helixradius");  // meter
    helixpitch = m_inputData.GetScalarOpt("helixpitch");    // meter
    rodRadius = m_inputData.GetScalarOpt("rodRadius");      // meter

    gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
    maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
	numVertices = m_inputData.GetIntOpt("numVertices");     // int_num
	youngM = m_inputData.GetScalarOpt("youngM");            // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");          // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");          // kg/m^3
	viscosity = m_inputData.GetScalarOpt("viscosity");      // viscosity in Pa-s

	baVector = m_inputData.GetVecOpt("baVector");           // magnetic field
	brVector = m_inputData.GetVecOpt("brVector");           // magnetic field
	muZero = m_inputData.GetScalarOpt("muZero");            // magnetic field

	dBar = m_inputData.GetScalarOpt("dBar");
	stiffness = m_inputData.GetScalarOpt("stiffness");

	scaleRender = m_inputData.GetScalarOpt("scaleRender");

	speed = m_inputData.GetScalarOpt("speed");
	
	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus
	
	// Viscous drag coefficients using Resistive Force Theory
	eta_per = 4.0*M_PI*viscosity/( log(2*helixpitch/rodRadius) + 0.5);
    eta_par = 2.0*M_PI*viscosity/( log(2*helixpitch/rodRadius) - 0.5);
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
	
	int systemRet = system("mkdir datafiles"); //make the directory

	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	ostringstream name;
	name.precision(6);
	name << fixed;
    name << "datafiles/simDER";
    name << "_rodLength_" << RodLength;
    name << "_Ba_" << baVector(0);
    name << ".txt";

    outfile.open(name.str().c_str());
    outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	//if (timeStep == Nstep)
	{
		outfile << currentTime << " " << m_externalContactForce->headContactForce << endl;
	}
}

void world::setRodStepper()
{
	// Set up geometry
	rodGeometry();	

	// Create the rod 
	rod = new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
		youngM, shearM, RodLength);

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = M_PI * pow(rodRadius ,4)/4.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;
	
	// Set up boundary condition
	rodBoundaryCondition();
	
	// setup the rod so that all the relevant variables are populated
	rod->setup();
	// End of rod setup
	
	// set up the time stepper
	stepper = new timeStepper(*rod);
	totalForce = stepper->getForce();

	// declare the forces
	m_stretchForce = new elasticStretchingForce(*rod, *stepper);
	m_bendingForce = new elasticBendingForce(*rod, *stepper);
	m_twistingForce = new elasticTwistingForce(*rod, *stepper);
	m_inertialForce = new inertialForce(*rod, *stepper);
	m_gravityForce = new externalGravityForce(*rod, *stepper, gVector);
	m_magneticForce = new externalMagneticForce(*rod, *stepper, baVector, brVector, muZero);
	m_dampingForce = new dampingForce(*rod, *stepper, viscosity, eta_per, eta_par);
	m_externalContactForce = new externalContactForce(*rod, *stepper, dBar, stiffness);

	Nstep = totalTime/deltaTime;

	// Allocate every thing to prepare for the first iteration
	rod->updateTimeStep();
	
	timeStep = 0;
	currentTime = 0.0;

	xInitial = rod->getVertex(rod->nv - 1);
}

// Setup geometry
void world::rodGeometry()
{

    vertices = MatrixXd(numVertices, 3);

    double deltaL = RodLength / (numVertices - 1);

    for (int i = 0; i < numVertices; i++)
    {
    	vertices(i, 0) = deltaL * i - RodLength;
        vertices(i, 1) = 0.0;
        vertices(i, 2) = 0.0;
    }

}

void world::rodBoundaryCondition()
{
	// Apply boundary condition
	rod->setVertexBoundaryCondition(rod->getVertex(0), 0);
	rod->setVertexBoundaryCondition(rod->getVertex(1), 1);
	rod->setThetaBoundaryCondition(0.0, 0);
	
	for (int i = 0; i < rod->nv; i++)
	{
		rod->setOneVertexBoundaryCondition(0.0, i, 2);
	}

	for (int i = 0; i < rod->ne; i++)
	{
		rod->setThetaBoundaryCondition(0.0, i);
	}
}
	

void world::updateTimeStep()
{
	Vector3d x1 = rod->getVertex(0);
	x1(0) = x1(0) + speed * deltaTime;
	rod->setVertexBoundaryCondition(x1, 0);

	Vector3d x2 = rod->getVertex(1);
	x2(0) = x2(0) + speed * deltaTime;
	rod->setVertexBoundaryCondition(x2, 1);

	double normf = forceTol * 10.0;	
	double normf0 = 0;
	
	bool solved = false;
	
	iter = 0;

	// Start with a trial solution for our solution x
	rod->updateGuess(); // x = x0 + u * dt
		
	while (solved == false)
	{
		rod->prepareForIteration();
		
		stepper->setZero();

		// Compute the forces and the jacobians
		m_inertialForce->computeFi();
		m_inertialForce->computeJi();
			
		m_stretchForce->computeFs();
		m_stretchForce->computeJs();
			
		m_bendingForce->computeFb();
		m_bendingForce->computeJb();
		
		m_twistingForce->computeFt();
		m_twistingForce->computeJt();

		m_gravityForce->computeFg();
		m_gravityForce->computeJg();

		m_magneticForce->computeFm(currentTime);
		m_magneticForce->computeJm();
		
		m_dampingForce->computeFd();
		m_dampingForce->computeJd();

		m_externalContactForce->computeFc();
		
		// Compute norm of the force equations.
		normf = 0;
		for (int i = 0; i < rod->uncons; i++)
		{
			normf += totalForce[i] * totalForce[i];
		}

		normf = sqrt(normf);
		if (iter == 0) normf0 = normf;
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}
		
		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			rod->updateNewtonX(totalForce); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			//cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}
	
	rod->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << " iter=" << iter << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;
	
	if (solved == false)
	{
		//timeStep = Nstep; // we are exiting
	}
}

int world::simulationRunning()
{
	if (timeStep < Nstep) 
	{
		return 1;
	}
	else 
	{
		return -1;
	}
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int i)
{
	return rod->x[i] * scaleRender;
}

Vector3d world::getX(int i)
{
	Vector3d xCurrent1 = rod->getVertex(i);
	Vector3d xCurrent2 = rod->getVertex(i);

	return ((xCurrent1 + xCurrent2) / 2) * scaleRender;
}

Vector3d world::getM1(int i)
{
	Vector3d xCurrent1 = rod->getVertex(i);
	Vector3d xCurrent2 = rod->getVertex(i);

	Vector3d m1 = rod->m1.row(i);

	return ((xCurrent1 + xCurrent2) / 2 + m1 * RodLength/50) * scaleRender;
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}

double world::getTubeRadius()
{
	return rod->tubeRadius * scaleRender;
}

int world::getTubeNv()
{
	return rod->tubeNv;
}

Vector3d world::getTubeNode(int i)
{
	Vector3d xCurrent = rod->tubeNode.row(i);

	return xCurrent * scaleRender;
}