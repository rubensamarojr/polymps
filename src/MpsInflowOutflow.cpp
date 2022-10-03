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
void MpsInflowOutflow::setPlaneFromPointNormal(MpsParticleSystem *PSystem, const int ii) {
	
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

	// Set plane variables
	Pio.ID = PSystem->inOutflowPlaneID[ii];
	Pio.press = PSystem->inOutflowPress[ii];
	Pio.velocity[0] = PSystem->inOutflowVel[ii*3  ];
	Pio.velocity[1] = PSystem->inOutflowVel[ii*3+1];
	Pio.velocity[2] = PSystem->inOutflowVel[ii*3+2];

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

// Calculates the euclidean signed distance between point P1 and the plane
double MpsInflowOutflow::calcSignedDistance(const double Px, const double Py, const double Pz) {
	return (Pio.a * Px + Pio.b * Py + Pio.c * Pz + Pio.d);
}

// Impose motion to particles
void MpsInflowOutflow::imposeMotionParticles(MpsParticleSystem *PSystem, MpsParticle *Particles) {
	
	double Amp[3], AccX, AccY, AccZ;
	Amp[0]= 0.1; Amp[1]= 0.0; Amp[2]= 0.0;	///< Amplitude motion (m)
	
	// Cte motion
	AccX = Amp[0];	AccY = Amp[1];	AccZ = Amp[2];

	// Oscillatory motion
	// double Per = 0.1;	///< Period (s)
	// double wo = 2*3.1416/Per;
	// AccX = -Amp[0] * wo * wo * sin(wo * PSystem->timeCurrent);
	// AccY = -Amp[1] * wo * wo * sin(wo * PSystem->timeCurrent);
	// AccZ = -Amp[2] * wo * wo * sin(wo * PSystem->timeCurrent);

#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
	if(Particles->particleType[i] == PSystem->fluid) {
		Particles->vel[i*3  ] += AccX * PSystem->timeStep;
		Particles->vel[i*3+1] += AccY * PSystem->timeStep;
		Particles->vel[i*3+2] += AccZ * PSystem->timeStep;
	}
	}
}

// Check particles in the Inflow/Outflow region
void MpsInflowOutflow::checkInOutflowParticles(MpsParticleSystem *PSystem, MpsParticle *Particles) {

	// Set number of inOutflow (IO) particles to Zero
	Particles->numIOParticles = 0;
	// Update number of particles including IO particles
	Particles->numRealAndIOParticles = Particles->numParticles + Particles->numIOParticles;

	// Auxiliar variables
	int numIOParticlesAux = 0;
	int iLastParticle = Particles->numParticles;

//#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
	if(Particles->particleType[i] == PSystem->fluid) {
		
		Particles->signDist[i] = calcSignedDistance(Particles->pos[i*3], Particles->pos[i*3+1], Particles->pos[i*3+2]);

		// Verify if the signed distance of particle i from the IOplan is positive
		if(Particles->signDist[i] > PSystem->epsilonZero) {
			// Particle inside the fluid domain
			
			// Verify if the signed distance of particle i from the IOplan is shorter than partDist
			if(Particles->signDist[i] < PSystem->partDist) {
				// Particle near the IOplan and its position shorter than partDist

				// Create two IO particles along the line normal to the IOplan and passing from the i-th position
				// IO particles are placed at distances partDist and 2*partdDist from i-th position
				double posIO1x, posIO1y, posIO1z, posIO2x, posIO2y, posIO2z;
				posIO1x = Particles->pos[i*3  ] - (PSystem->partDist)*Pio.a;
				posIO1y = Particles->pos[i*3+1] - (PSystem->partDist)*Pio.b;
				posIO1z = Particles->pos[i*3+2] - (PSystem->partDist)*Pio.c;
				posIO2x = Particles->pos[i*3  ] - (2.0*PSystem->partDist)*Pio.a;
				posIO2y = Particles->pos[i*3+1] - (2.0*PSystem->partDist)*Pio.b;
				posIO2z = Particles->pos[i*3+2] - (2.0*PSystem->partDist)*Pio.c;

				// Assign velocities equal to those of i-th particle
				double velIO1x, velIO1y, velIO1z, velIO2x, velIO2y, velIO2z;
				velIO1x = Particles->vel[i*3  ];
				velIO1y = Particles->vel[i*3+1];
				velIO1z = Particles->vel[i*3+2];
				velIO2x = Particles->vel[i*3  ];
				velIO2y = Particles->vel[i*3+1];
				velIO2z = Particles->vel[i*3+2];
				
				Particles->pos[iLastParticle*3  ] = posIO1x;
				Particles->pos[iLastParticle*3+1] = posIO1y;
				Particles->pos[iLastParticle*3+2] = posIO1z;
				Particles->vel[iLastParticle*3  ] = velIO1x;
				Particles->vel[iLastParticle*3+1] = velIO1y;
				Particles->vel[iLastParticle*3+2] = velIO1z;

				// Set some variables
				Particles->particleBC[iLastParticle] = PSystem->other;
				Particles->particleType[iLastParticle] = PSystem->inOutflowPartID;
				Particles->press[iLastParticle] = Pio.press;

				// Update lastParticle ID
				iLastParticle++;

				Particles->pos[iLastParticle*3  ] = posIO2x;
				Particles->pos[iLastParticle*3+1] = posIO2y;
				Particles->pos[iLastParticle*3+2] = posIO2z;
				Particles->vel[iLastParticle*3  ] = velIO2x;
				Particles->vel[iLastParticle*3+1] = velIO2y;
				Particles->vel[iLastParticle*3+2] = velIO2z;

				// Set some variables
				Particles->particleBC[iLastParticle] = PSystem->other;
				Particles->particleType[iLastParticle] = PSystem->inOutflowPartID;
				Particles->press[iLastParticle] = Pio.press;
				
				// Update lastParticle ID
				iLastParticle++;

				numIOParticlesAux += 2;
			}
			// Verify if the signed distance of particle i from the IOplan is shorter than reS
			else if(Particles->signDist[i] < PSystem->reS) {
				// Particle near the IOplan and its position ranges between partDist and reS

				// Computes ui . normal
				double ui_n = Particles->vel[i*3] * Pio.normal[0] + Particles->vel[i*3+1] * Pio.normal[1] + Particles->vel[i*3+2] * Pio.normal[2];

				// Verify the velocity component normal to the IOplan
				if(ui_n > PSystem->epsilonZero) {
					// The velocity component normal to the IOplan is positive, then apply Inflow condition

					// Check if the release of a new particle is required
					// Conical scan region is applied
					
					// The direction of the conical region axis is the sum of the vectors rij pointing towards i-th particle

					// Effective neighboring j-th particles are verified if is outside the scan region
				}
				else{
					// The velocity component normal to the IOplan is negative, then apply Outflow condition
					// No new effective particle has to be generated

				}
			}
		}
		else {
			// The signed distance is negative, then the Particle is outside the fluid domain
			// Set particle as ghost
		}
		
	}
	}


	// Updates number of inOutflow (IO) particles
	Particles->numIOParticles = numIOParticlesAux;
	// Updates number of particles including IO particles
	Particles->numRealAndIOParticles = Particles->numParticles + Particles->numIOParticles;
}

// Clear variables only of the inOutflow (IO) particles
void MpsInflowOutflow::clearVariablesInOutflowParticles(MpsParticleSystem *PSystem, MpsParticle *Particles) {
//#pragma omp parallel for schedule(dynamic,64)
	for(int i=Particles->numParticles-1; i<Particles->numRealAndIOParticles; i++) {
	
	// Update some data of IO Particle
	Particles->particleType[i]=PSystem->ghost;
	Particles->particleBC[i]=PSystem->other;
	Particles->particleNearWall[i]=false;
	Particles->nearMeshType[i]=meshType::FIXED;
	Particles->distParticleWall2[i]=10e8*PSystem->partDist;
	// Set zero to velocity and press of lastParticle
	for (int j = 0; j < 3; j++){
		//pos[i*3+j]=0.0;
		//Posk[i*3+j]=0.0;
		Particles->vel[i*3+j] = 0.0;
	}
	Particles->press[i] = 0.0;
	// Set maximum position to lastParticle
	Particles->pos[i*3  ] = PSystem->domainMaxX - PSystem->partDist;
	Particles->pos[i*3+1] = PSystem->domainMaxY - PSystem->partDist;
	if (PSystem->dim == 2)
		Particles->pos[i*3+2] = 0.0;
	else
		Particles->pos[i*3+2] = PSystem->domainMaxZ - PSystem->partDist;
	
	}
}