// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <experimental/filesystem> 	///< numeric_limits
#include <iostream>					///< strings and c-strings
#include "MpsInflowOutflow.h"

using namespace std;

// Constructor declaration
MpsInflowOutflow::MpsInflowOutflow()
{
}
// Destructor declaration
MpsInflowOutflow::~MpsInflowOutflow()
{
}

/// Constructs a plane in 3D space by specifying a single point on the plane and the surface normal pointing toward the fluid domain 
void MpsInflowOutflow::setPlaneFromPointNormal(MpsParticleSystem *PSystem, const int ii)
{
	// Normalize the normal vector
	double normN = sqrt(PSystem->inOutflowNormal[ii*3  ]*PSystem->inOutflowNormal[ii*3  ] + 
						PSystem->inOutflowNormal[ii*3+1]*PSystem->inOutflowNormal[ii*3+1] + 
						PSystem->inOutflowNormal[ii*3+2]*PSystem->inOutflowNormal[ii*3+2]);

	if(normN > PSystem->epsilonZero) {
		Pio.normal[0] = PSystem->inOutflowNormal[ii*3  ]/normN;
		Pio.normal[1] = PSystem->inOutflowNormal[ii*3+1]/normN;
		Pio.normal[2] = PSystem->inOutflowNormal[ii*3+2]/normN;
	}
	else {
		for(unsigned int i = 0; i < 3; i++)
			Pio.normal[i] = 0.0;
	}
	Pio.a = Pio.normal[0];
	Pio.b = Pio.normal[1];
	Pio.c = Pio.normal[2];
	//Pio.d = - (Pt[0]*Pio.normal[0] + Pt[1]*Pio.normal[1] + Pt[2]*Pio.normal[2]);
	Pio.d = (	PSystem->inOutflowPt[ii*3  ]*Pio.normal[0] + 
				PSystem->inOutflowPt[ii*3+1]*Pio.normal[1] + 
				PSystem->inOutflowPt[ii*3+2]*Pio.normal[2]);

	printf("Inflow/Outflow id: %2d", PSystem->inOutflowPlaneID[ii]);
	printf(", Seed Point: %3.4f %3.4f %3.4f", PSystem->inOutflowPt[ii*3], PSystem->inOutflowPt[ii*3+1], PSystem->inOutflowPt[ii*3+2]);
	printf(", Normal vector: %3.4f %3.4f %3.4f", PSystem->inOutflowNormal[ii*3], PSystem->inOutflowNormal[ii*3+1], PSystem->inOutflowNormal[ii*3+2]);
	
	if (PSystem->inOutflowTypeBC[ii] == 0) {
		printf(", Uniform inlet velocity: %3.4f %3.4f %3.4f\n", PSystem->inOutflowVel[ii*3], PSystem->inOutflowVel[ii*3+1], PSystem->inOutflowVel[ii*3+2]);
	}
	else 
		if (PSystem->inOutflowTypeBC[ii] == 1) {
		printf(", Constant Pressure: %3.4f\n", PSystem->inOutflowPress[ii]);
	}				
}

// Calculates the euclidean distance between point P1 and the plane
double MpsInflowOutflow::calcSignedDistance(const double *P1)
{
	return (Pio.a * P1[0] + Pio.b * P1[1] + Pio.c * P1[2] + Pio.d);
}