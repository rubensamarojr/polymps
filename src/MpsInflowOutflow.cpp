// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <experimental/filesystem> 	///< numeric_limits
#include <iostream>					///< strings and c-strings
#include "MpsInflowOutflow.h"

#include <sys/time.h>	///< time
#include <iomanip>		///< std::setprecision

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

// Check real particles in the Inflow/Outflow region and create real and IO particles
void MpsInflowOutflow::checkRealParticlesInOutflow(MpsParticleSystem *PSystem, MpsParticle *Particles) {

	//double minDistIOcreation = 1.1*PSystem->partDist;
	// Auxiliar variables
	int numCreatedParticlesAux, idLastRealParticle, idLastRealandIOParticle;
//#pragma omp parallel
//{
	// Set auxiliar variables
	numCreatedParticlesAux = 0;						///< Number of created particles
	idLastRealParticle = Particles->numParticles;	///< Array ID of the last real particle + 1
	idLastRealandIOParticle = Particles->numParticles + numCreatedParticlesAux;	///< Array ID of the last real+created particle + 1

	// Set number of created particles to zero
	Particles->numCreatedParticles = 0;

	int sumAux, sumAux1, sumAux2;
	sumAux = sumAux1 = sumAux2 = 0;

	struct timeval time_now, time_end;
	double time_taken;

// 	// start timer.
// 	gettimeofday(&time_now, NULL);
// //#pragma omp parallel for reduction(+:sumAux, sumAux1, sumAux2)
// // #pragma omp parallel for
// 	for(int i=0; i<Particles->numParticles; i++) {
// 	if(Particles->particleType[i] == PSystem->fluid) {

// 		if(i%2==0||i%3==0)
// 		{
// // #pragma omp critical
// {
// // 			int tid = omp_get_thread_num();
// // 			printf("\nomp thread %d", tid);
//  			// int i1, i2;
//  			// i1 = sumAux;
//  			// i2 = sumAux+1;
// // #pragma omp atomic
// 			sumAux+=2;
// // #pragma omp atomic
// 			sumAux1+=4;
// // #pragma omp atomic
// 			sumAux2+=3;
// } // end omp critical
// 		}
// 	}
// }
// 	// stop timer.
// 	gettimeofday(&time_end, NULL);
// 	// Calculating total time taken by the program.
// 	time_taken = (time_end.tv_sec - time_now.tv_sec) * 1e6;
// 	time_taken = (time_taken + (time_end.tv_usec - time_now.tv_usec)) * 1e-6;
// 	printf("\n Time: %.5e sec - sumA: %d - sumA1: %d - sumA2: %d", time_taken, sumAux, sumAux1, sumAux2);


	// start timer.
	gettimeofday(&time_now, NULL);

// Be aware with "race condition" in OpenMP
//#pragma omp parallel for schedule(dynamic,64) shared (idLastRealandIOParticle) reduction(+:numCreatedParticlesAux)
//#pragma omp for nowait
//#pragma omp parallel for shared (idLastRealandIOParticle) reduction(+:numCreatedParticlesAux)
//#pragma omp parallel for reduction(+:numCreatedParticlesAux)
//#pragma omp parallel for schedule(dynamic,64) reduction(+:numCreatedParticlesAux)
//#pragma omp parallel for reduction(+:sumAux)
//#pragma omp parallel for reduction(+:numCreatedParticlesAux, sumAux)
#pragma omp parallel for reduction(+:numCreatedParticlesAux)
	for(int i=0; i<Particles->numParticles; i++) {
	if(Particles->particleType[i] == PSystem->fluid) {
		
		Particles->signDist[i] = calcSignedDistance(Particles->pos[i*3], Particles->pos[i*3+1], Particles->pos[i*3+2]);

// To prevent race condition, the access to the shared (global) variable "idLastRealandIOParticle" must be synchronized.
//#pragma omp critical (before_if)
//{
		// Verify if the signed distance of particle i from the IOplan is positive
		if(Particles->signDist[i] > PSystem->epsilonZero) {
			// Particle inside the fluid domain
			
			// Verify if the signed distance of particle i from the IOplan is shorter than partDist
			if(Particles->signDist[i] <= PSystem->partDist) {
				// Particle near the IOplan and its position shorter than partDist

				// Assign true to particle is in the IOregion
				Particles->isInIORegion[i] = true;

// Specifies that code is only executed on one thread at a time
//#pragma omp critical (after_if)
//{
				// Create two IO particles along the line normal to the IOplan and passing from the i-th position
				// Each thread of the "parallel for" has private (local) variables of two IO particles index
				// Be aware that they are defined using the shared (global) variable idLastRealandIOParticle
				int idIO1, idIO2;

// Specifies that code is only executed on one thread at a time
// #pragma omp critical (after_if)
// {
// 				idIO1 = idLastRealandIOParticle;
// 				idIO2 = idIO1+1;
// 				// Update shared (global) variables lastParticle ID and number of IO particles
// 				idLastRealandIOParticle += 2;
// 				// numCreatedParticlesAux += 2;
// } //end of omp critical (after_if)

// idIO1 atomically captures the original value of idLastRealandIOParticle. Then, atomically updates idLastRealandIOParticle
#pragma omp atomic capture
				{
					idIO1 = idLastRealandIOParticle;
					idLastRealandIOParticle += 2;
				}
				idIO2 = idIO1+1;

				// IO particles are placed at distances partDist and 2*partdDist from i-th position
				Particles->pos[idIO1*3  ] = Particles->pos[i*3  ] - (PSystem->partDist)*Pio.a;
				Particles->pos[idIO1*3+1] = Particles->pos[i*3+1] - (PSystem->partDist)*Pio.b;
				Particles->pos[idIO1*3+2] = Particles->pos[i*3+2] - (PSystem->partDist)*Pio.c;
				Particles->pos[idIO2*3  ] = Particles->pos[i*3  ] - (2.0*PSystem->partDist)*Pio.a;
				Particles->pos[idIO2*3+1] = Particles->pos[i*3+1] - (2.0*PSystem->partDist)*Pio.b;
				Particles->pos[idIO2*3+2] = Particles->pos[i*3+2] - (2.0*PSystem->partDist)*Pio.c;
				
				// Assign velocities equal to those of i-th particle
				Particles->vel[idIO1*3  ] = Particles->vel[i*3  ];
				Particles->vel[idIO1*3+1] = Particles->vel[i*3+1];
				Particles->vel[idIO1*3+2] = Particles->vel[i*3+2];
				Particles->vel[idIO2*3  ] = Particles->vel[i*3  ];
				Particles->vel[idIO2*3+1] = Particles->vel[i*3+1];
				Particles->vel[idIO2*3+2] = Particles->vel[i*3+2];

				// Set some variables
				Particles->particleBC[idIO1] = PSystem->other;
				Particles->particleType[idIO1] = PSystem->inOutflowPartID;
				Particles->press[idIO1] = Pio.press;
				Particles->particleBC[idIO2] = PSystem->other;
				Particles->particleType[idIO2] = PSystem->inOutflowPartID;
				Particles->press[idIO2] = Pio.press;

				// Update lastParticle ID
// The atomic construct allows multiple threads to safely update a shared (global) variable
//#pragma omp atomic
				//idLastRealandIOParticle++;

				// Update shared (global) variable number of IO particles
				numCreatedParticlesAux += 2;

//} //end of omp critical (after_if)
			}
			// Verify if the signed distance of particle i from the IOplan is shorter than reS
			else if(Particles->signDist[i] < PSystem->reS) {
				// Particle near the IOplan and its position ranges between partDist and reS

				// Verify if the particle was in the IOregion in the previous step
				if(Particles->isInIORegion[i] == true) {

					// Assign false to particle is in the IOregion
					Particles->isInIORegion[i] = false;

					// Create one Real particle and two IO particles along the line normal to the IOplan and passing from the idReal-th position
					// Each thread of the "parallel for" has private (local) variables of two IO particles index
					// Be aware that they are defined using the shared (global) variable idLastRealandIOParticle
					int idReal, idIO1, idIO2;

// idReal atomically captures the original value of idLastRealandIOParticle. Then, atomically updates idLastRealandIOParticle
#pragma omp atomic capture
					{
						idReal = idLastRealandIOParticle;
						idLastRealandIOParticle += 3;
					}
					idIO1 = idReal+1;
					idIO2 = idIO1+1;
					
					// Real particle is placed at distance partDist from i-th position
					Particles->pos[idReal*3  ] = Particles->pos[i*3  ] - (PSystem->partDist)*Pio.a;
					Particles->pos[idReal*3+1] = Particles->pos[i*3+1] - (PSystem->partDist)*Pio.b;
					Particles->pos[idReal*3+2] = Particles->pos[i*3+2] - (PSystem->partDist)*Pio.c;
					// Assign velocities equal to those of i-th particle
					Particles->vel[idReal*3  ] = Particles->vel[i*3  ];
					Particles->vel[idReal*3+1] = Particles->vel[i*3+1];
					Particles->vel[idReal*3+2] = Particles->vel[i*3+2];

					// IO particles are placed at distances partDist and 2*partdDist from idReal-th position
					Particles->pos[idIO1*3  ] = Particles->pos[idReal*3  ] - (PSystem->partDist)*Pio.a;
					Particles->pos[idIO1*3+1] = Particles->pos[idReal*3+1] - (PSystem->partDist)*Pio.b;
					Particles->pos[idIO1*3+2] = Particles->pos[idReal*3+2] - (PSystem->partDist)*Pio.c;
					Particles->pos[idIO2*3  ] = Particles->pos[idReal*3  ] - (2.0*PSystem->partDist)*Pio.a;
					Particles->pos[idIO2*3+1] = Particles->pos[idReal*3+1] - (2.0*PSystem->partDist)*Pio.b;
					Particles->pos[idIO2*3+2] = Particles->pos[idReal*3+2] - (2.0*PSystem->partDist)*Pio.c;
					// Assign velocities equal to those of idReal-th particle
					Particles->vel[idIO1*3  ] = Particles->vel[idReal*3  ];
					Particles->vel[idIO1*3+1] = Particles->vel[idReal*3+1];
					Particles->vel[idIO1*3+2] = Particles->vel[idReal*3+2];
					Particles->vel[idIO2*3  ] = Particles->vel[idReal*3  ];
					Particles->vel[idIO2*3+1] = Particles->vel[idReal*3+1];
					Particles->vel[idIO2*3+2] = Particles->vel[idReal*3+2];

					// Set some variables
					Particles->particleBC[idReal] = PSystem->inner;
					Particles->particleBC[idIO1] = PSystem->other;
					Particles->particleBC[idIO2] = PSystem->other;
					
					Particles->particleType[idReal] = PSystem->fluid;
					Particles->particleType[idIO1] = PSystem->inOutflowPartID;
					Particles->particleType[idIO2] = PSystem->inOutflowPartID;
					
					Particles->press[idReal] = Pio.press;
					Particles->press[idIO1] = Pio.press;
					Particles->press[idIO2] = Pio.press;
					
					Particles->isInIORegion[idReal] = true;
					Particles->isInIORegion[i] = false;
					Particles->isInIORegion[i] = false;

					// Update lastParticle ID
	// The atomic construct allows multiple threads to safely update a shared (global) variable
	//#pragma omp atomic
					//idLastRealandIOParticle++;

					// Update shared (global) variable number of created particles
					numCreatedParticlesAux += 3;

				}
				else {
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
		}
		else {
			// The signed distance is negative, then the Particle is outside the fluid domain
			// Set particle as ghost
		}
	
//} //end of omp critical (before_if)

	} // end of if(Particles->particleType[i]...
	} // end of for(int i=0;...

//} //end of omp parallel

	// if(numCreatedParticlesAux != 0 && numCreatedParticlesAux != 200)
		// printf("\n NumParticles: %d", numCreatedParticlesAux);

	// printf("\n NumParticles: %d", sumAux);
	
	// Updates number of Real particles
	// Particles->numParticles = numRealParticlesAux;
	// Updates number of inOutflow (IO) particles
	Particles->numCreatedParticles = numCreatedParticlesAux;
	// Updates number of particles including IO particles
	Particles->numRealAndIOParticles = Particles->numParticles + Particles->numCreatedParticles;

	// stop timer.
	gettimeofday(&time_end, NULL);
	// Calculating total time taken by the program.
	time_taken = (time_end.tv_sec - time_now.tv_sec) * 1e6;
	time_taken = (time_taken + (time_end.tv_usec - time_now.tv_usec)) * 1e-6;
	printf("\n Time: %.5e sec - NumPart: %d - numCreatedPart: %d - numRealAndIOPart: %d - idLastRealandIOParticle: %d", time_taken, Particles->numParticles, Particles->numCreatedParticles, Particles->numRealAndIOParticles, idLastRealandIOParticle);
}


// Swap ID of real particles in the array Particles->numRealAndIOParticles to the array Particles->numParticles
void MpsInflowOutflow::swapIdRealAndIOParticlesInOutflow(MpsParticleSystem *PSystem, MpsParticle *Particles) {

	// Auxiliar variables
	int numCreatedRealParticlesAux, idLastRealParticle;

	// Set auxiliar variables considering that all created particles are IO particles, i.e.,
	// the number of created real particles is zero
	Particles->numIOParticles = Particles->numCreatedParticles;	///< Number of inOutflow (IO) particles
	numCreatedRealParticlesAux = 0;								///< Number of created Real particles
	idLastRealParticle = Particles->numParticles;				///< Array ID of the last real particle + 1

	struct timeval time_now, time_end;
	double time_taken;

	// start timer.
	gettimeofday(&time_now, NULL);

// Be aware with "race condition" in OpenMP
// Loop only for the created particles
// #pragma omp parallel for reduction(+:numCreatedRealParticlesAux)
	for(int i=0; i<Particles->numCreatedParticles; i++) {

		int idIO = i + Particles->numParticles;
		// Verify if the created particle is real
		if(Particles->particleBC[idIO] == PSystem->inner) {
		
		// Array ID to be filled with a Real particle
		int idReal;

// idReal atomically captures the original value of idLastRealParticle. Then, atomically updates idLastRealParticle
// #pragma omp atomic capture
			idReal = idLastRealParticle++;
			
			// if(Particles->particleBC[idReal] != PSystem->inner)
			// {
				// Update the variables in the array
				// Swapping elements of the new real particle and IO
				swap(Particles->pos[idReal*3  ], Particles->pos[idIO*3  ]);
				swap(Particles->pos[idReal*3+1], Particles->pos[idIO*3+1]);
				swap(Particles->pos[idReal*3+2], Particles->pos[idIO*3+2]);
				swap(Particles->vel[idReal*3  ], Particles->vel[idIO*3  ]);
				swap(Particles->vel[idReal*3+1], Particles->vel[idIO*3+1]);
				swap(Particles->vel[idReal*3+2], Particles->vel[idIO*3+2]);
				swap(Particles->particleBC[idReal], Particles->particleBC[idIO]);
				swap(Particles->particleType[idReal], Particles->particleType[idIO]);
				swap(Particles->press[idReal], Particles->press[idIO]);
			// }
			
			// Update shared (global) variable number of Real, IO particles
			numCreatedRealParticlesAux++;
		}
	
	} // end for(int idIO=0;...
	
	// Updates number of Real particles
	Particles->numParticles += numCreatedRealParticlesAux;
	// Updates number of inOutflow (IO) particles
	Particles->numIOParticles -= numCreatedRealParticlesAux;
	// Updates number of particles including IO particles
	// Particles->numRealAndIOParticles = Particles->numParticles + Particles->numIOParticles;

	// stop timer.
	gettimeofday(&time_end, NULL);
	// Calculating total time taken by the program.
	time_taken = (time_end.tv_sec - time_now.tv_sec) * 1e6;
	time_taken = (time_taken + (time_end.tv_usec - time_now.tv_usec)) * 1e-6;
	printf("\n Time: %.5e sec - NumPart: %d - NumIOPart: %d - NumRealIoPart: %d - idLastRealParticle: %d", time_taken, Particles->numParticles, Particles->numIOParticles, Particles->numRealAndIOParticles, idLastRealParticle);

}


/////////////////////////////////////////////////
// Check IO particles in the Inflow/Outflow region
void MpsInflowOutflow::checkIOParticlesInOutflow(MpsParticleSystem *PSystem, MpsParticle *Particles) {

	// Set number of inOutflow (IO) particles to Zero
	//Particles->numIOParticles = 0;
	// Update number of particles including IO particles
	//Particles->numRealAndIOParticles = Particles->numParticles + Particles->numIOParticles;

	// Auxiliar variables
	int numRealParticlesAux, numIOParticlesAux, idLastRealParticle, idLastRealandIOParticle;

	// Set auxiliar variables
	numRealParticlesAux = Particles->numParticles;	///< Number of Real particles
	numIOParticlesAux = Particles->numIOParticles;	///< Number of inOutflow (IO) particles
	idLastRealParticle = Particles->numParticles;	///< Array ID of the last real particle + 1
	idLastRealandIOParticle = Particles->numRealAndIOParticles;	///< Array ID of the last real+IO particle + 1

	struct timeval time_now, time_end;
	double time_taken;

	// start timer.
	gettimeofday(&time_now, NULL);

// Be aware with "race condition" in OpenMP
// Loop only for IO particles
//#pragma omp parallel for reduction(+:numIOParticlesAux)
	for(int idIO=Particles->numParticles; idIO<Particles->numRealAndIOParticles; idIO++) {
	// if(Particles->particleType[idIO] == PSystem->fluid) {
		
		Particles->signDist[idIO] = calcSignedDistance(Particles->pos[idIO*3], Particles->pos[idIO*3+1], Particles->pos[idIO*3+2]);

// To prevent race condition, the access to the shared (global) variable "idLastRealParticle" must be synchronized.
//#pragma omp critical (before_if)
//{
		// Verify if the signed distance of IO particle from the IOplan is positive
		// if(Particles->signDist[idIO] >= PSystem->epsilonZero) {
		if(Particles->signDist[idIO] >= 0) {
			// IO particle inside the fluid domain

			// Array ID to be filled with a Real particle
			int idReal;

// idReal atomically captures the original value of idLastRealParticle. Then, atomically updates idLastRealParticle
// #pragma omp atomic capture
			{
				idReal = idLastRealParticle;
				idLastRealParticle += 1;
			}
			
			// Update the variables in the array
			// Swapping elements of the new real particle and IO
			swap(Particles->pos[idReal*3  ], Particles->pos[idIO*3  ]);
			swap(Particles->pos[idReal*3+1], Particles->pos[idIO*3+1]);
			swap(Particles->pos[idReal*3+2], Particles->pos[idIO*3+2]);
			swap(Particles->vel[idReal*3  ], Particles->vel[idIO*3  ]);
			swap(Particles->vel[idReal*3+1], Particles->vel[idIO*3+1]);
			swap(Particles->vel[idReal*3+2], Particles->vel[idIO*3+2]);
			Particles->particleBC[idReal] = PSystem->inner;
			Particles->particleType[idReal] = PSystem->fluid;
			swap(Particles->press[idReal], Particles->press[idIO]);

			// Create one IO particle along the line normal to the IOplan and passing from the new idReal-th particle's position
			// Each thread of the "parallel for" has private (local) variables of IO particle index
			// Be aware that it is defined using the shared (global) variable idLastRealandIOParticle
			int idIO2;

// idIO2 atomically captures the original value of idLastRealandIOParticle. Then, atomically updates idLastRealandIOParticle
// #pragma omp atomic capture
			{
				idIO2 = idLastRealandIOParticle;
				idLastRealandIOParticle += 1;
			}

			// IO particle is placed at distance partDist from idReal-th position
			Particles->pos[idIO2*3  ] = Particles->pos[idReal*3  ] - (2.0*PSystem->partDist)*Pio.a;
			Particles->pos[idIO2*3+1] = Particles->pos[idReal*3+1] - (2.0*PSystem->partDist)*Pio.b;
			Particles->pos[idIO2*3+2] = Particles->pos[idReal*3+2] - (2.0*PSystem->partDist)*Pio.c;
			
			// Assign velocities equal to those of idReal-th particle
			Particles->vel[idIO2*3  ] = Particles->vel[idReal*3  ];
			Particles->vel[idIO2*3+1] = Particles->vel[idReal*3+1];
			Particles->vel[idIO2*3+2] = Particles->vel[idReal*3+2];

			// Set some variables
			Particles->particleBC[idIO2] = PSystem->other;
			Particles->particleType[idIO2] = PSystem->inOutflowPartID;
			Particles->press[idIO2] = Pio.press;

			// Update shared (global) variable number of Real, IO particles
			numRealParticlesAux += 1;
			// numIOParticlesAux += 1;
			
		}
	
//} //end of omp critical (before_if)

	// }
	} // end for(int i=0;...

//} //end of omp parallel
	
	// stop timer.
	gettimeofday(&time_end, NULL);
	// Calculating total time taken by the program.
	time_taken = (time_end.tv_sec - time_now.tv_sec) * 1e6;
	time_taken = (time_taken + (time_end.tv_usec - time_now.tv_usec)) * 1e-6;
	printf("\n Time: %.5e sec - NumParticles: %d - idLastRealParticle: %d", time_taken, numIOParticlesAux, idLastRealParticle);

	// Updates number of Real particles
	Particles->numParticles = numRealParticlesAux;
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