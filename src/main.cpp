/**
 * simDER
 * simDER stands for "[sim]plified [D]iscrete [E]lastic [R]ods"
 * Dec 2017
 * This code is based on previous iterations. 
 * */

//This line is for mac
//#include <GLUT/glut.h>

//This is for linux
#include <GL/glut.h>

#include <iostream>
#include <fstream>
#include <string>
#include "eigenIncludes.h"

// Rod and stepper are included in the world
#include "world.h"
#include "initialization/setInput.h"

world myWorld;
int NPTS;
ofstream outfile;

static void Key(unsigned char key, int x, int y)
{
  switch (key) // ESCAPE to quit
  {
	case 27:
		exit(0);
  }
}

/* Initialize OpenGL Graphics */
void initGL() 
{
	glClearColor(0.7f, 0.7f, 0.7f, 0.0f); // Set background color to black and opaque
	glClearDepth(10.0f);                   // Set background depth to farthest
	//glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	//glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);   // Enable smooth shading
	//glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections

	glLoadIdentity();
	//gluLookAt(0.00, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	gluLookAt(0.05, 0.05, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
	glPushMatrix();

	//glMatrixMode(GL_MODELVIEW);
}

void display(void)
{
	while ( myWorld.simulationRunning() > 0)
	{
		//  Clear screen and Z-buffer
		glClear(GL_COLOR_BUFFER_BIT);

		// draw axis
		double axisLen = 1;
		glLineWidth(0.5);
		
		glBegin(GL_LINES);
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(-axisLen, 0.0, 0.0);
			glVertex3f(axisLen, 0.0, 0.0);

			glColor3f(0.0, 1.0, 0.0);
			glVertex3f(0.0, -axisLen, 0.0);
			glVertex3f(0.0, axisLen, 0.0);

			glColor3f(0.0, 0.0, 1.0);
			glVertex3f(0.0, 0.0, -axisLen);
			glVertex3f(0.0, 0.0, axisLen);
		glEnd();
		
		//draw a line
		glColor3f(1.0,0.0,0.0);
		glLineWidth(3.0);
		
		glBegin(GL_LINES);
		for (int i=0; i < NPTS-1; i++)
		{
			glVertex3f( myWorld.getScaledCoordinate(4*i), myWorld.getScaledCoordinate(4*i+1), myWorld.getScaledCoordinate(4*i+2));
			glVertex3f( myWorld.getScaledCoordinate(4*(i+1)), myWorld.getScaledCoordinate(4*(i+1)+1), myWorld.getScaledCoordinate(4*(i+1)+2));
		}
		glEnd();

		glColor3f(0.1, 0.1, 0.1);
		glPointSize(1.0);
		for (int i=0; i < myWorld.getTubeNv()-2; i = i + 1)
		{
			Vector3d xCurrent1 = myWorld.getTubeNode(i);
			Vector3d xCurrent2 = myWorld.getTubeNode(i+1);

			Vector3d xAverage = (xCurrent1 + xCurrent2) / 2;
			Vector3d xTangent = (xCurrent2 - xCurrent1) / (xCurrent2 - xCurrent1).norm();

			Vector3d xRand;

			xRand(0) = 0.1;
			xRand(1) = 0.1;
			xRand(2) = 0.1;

			xRand = xRand / xRand.norm();

			Vector3d xM1 = xTangent.cross(xRand);
			xM1 = xM1 / xM1.norm();
			Vector3d xM2 = xM1.cross(xTangent);
			xM2 = xM2 / xM2.norm();

			double deltaKK = 2 * M_PI / 100;

			glBegin(GL_POINTS);

			for (int kk = 0; kk < 100; kk++)
			{
				Vector3d xPlot = xAverage + myWorld.getTubeRadius() * ( xM1 * cos(kk * deltaKK) + xM2 * sin(kk * deltaKK) );
				glVertex3f( xPlot(0), xPlot(1), xPlot(2));
			}

			glEnd();
		}
		
		glFlush();
		
		// Update step
		myWorld.updateTimeStep();
		myWorld.CoutData(outfile);
	}
	exit(1);
}

int main(int argc,char *argv[])
{
	setInput inputData;
	inputData = setInput();
	inputData.LoadOptions(argv[1]);
	inputData.LoadOptions(argc,argv);
	//read input parameters from txt file and cmd

	myWorld = world(inputData);
	myWorld.setRodStepper();

	myWorld.OpenFile(outfile);

	bool render = myWorld.isRender();
	if (render) // if OpenGL visualization is on
	{
		NPTS = myWorld.numPoints();
	
		glutInit(&argc,argv);
		glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize (1000, 1000);
		glutInitWindowPosition (100, 100);
		glutCreateWindow ("simDER");
		initGL();
		glutKeyboardFunc(Key);
		glutDisplayFunc(display);
		glutMainLoop();
	}	
	else
	{
		while ( myWorld.simulationRunning() > 0)
		{
			myWorld.updateTimeStep(); // update time step
			myWorld.CoutData(outfile); // write data to file
		}
	}

	// Close (if necessary) the data file
	myWorld.CloseFile(outfile);
	
	return 0;
}

