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
	Pio.d = - (	PSystem->inOutflowPt[ii*3  ]*Pio.normal[0] + 
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
	
	int motType = 1;

	double Amp[3], AccX, AccY, AccZ;
	
	// Cte motion
	if (motType == 0) {
		Amp[0] = 0.1;	Amp[1] = 0.0; Amp[2] = 0.0;	///< Amplitude motion (m)
		AccX = Amp[0];	AccY = Amp[1];	AccZ = Amp[2];
	}
	// Oscillatory motion
	else if (motType == 1)
	{
		Amp[0] = 0.25;	Amp[1] = 0.0; Amp[2] = 0.0;	///< Amplitude motion (m)
		double Per = 0.25;	///< Period (s)
		double wo = 2*3.1416/Per;
		AccX = Amp[0] * wo * wo * (-0.1 + sin(wo * PSystem->timeCurrent));
		AccY = Amp[1] * wo * wo * sin(wo * PSystem->timeCurrent);
		AccZ = Amp[2] * wo * wo * sin(wo * PSystem->timeCurrent);
	}

#pragma omp parallel for
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		if(Particles->particleType[i] == PSystem->fluid) {
			if (motType == 0) {
				Particles->acc[i*3  ] += AccX;
				// Particles->vel[i*3  ] += AccX * PSystem->timeStep;
			}
			else {
				Particles->acc[i*3  ] += AccX * (Particles->pos[i*3+2] - 0.25);
				// Particles->vel[i*3  ] += AccX * (Particles->pos[i*3+2] - 0.25)/0.045 * PSystem->timeStep;
			}

			Particles->acc[i*3+1] += AccY;
			Particles->acc[i*3+2] += AccZ;

			// Particles->vel[i*3+1] += AccY * PSystem->timeStep;
			// Particles->vel[i*3+2] += AccZ * PSystem->timeStep;

			// Particles->pos[i*3  ] += Particles->vel[i*3  ]*PSystem->timeStep;
			// Particles->pos[i*3+1] += Particles->vel[i*3+1]*PSystem->timeStep;
			// Particles->pos[i*3+2] += Particles->vel[i*3+2]*PSystem->timeStep;
		}
	}
}

// Set initial values to inOutflow variables
void MpsInflowOutflow::setInOutflowVariables(MpsParticleSystem *PSystem, MpsParticle *Particles) {
	// Set number of created and deleted particles to zero
	Particles->numCreatedParticles = 0;
	Particles->numDeletedParticles = 0;
	// Set number of Real + IO particles to number of Real particles, i.e.,
	// the number of IO particles is zero
	Particles->numRealAndIOParticles = Particles->numParticles;
	Particles->numIOParticles = 0;

	// Set value of euclidean signed distance to infinity
#pragma omp parallel for
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		if(Particles->particleType[i] == PSystem->fluid) {
			Particles->signDist[i] = PSystem->nearInfinity;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Check particles in the Inflow/Outflow region, create real and IO particles, or delete real particles
void MpsInflowOutflow::checkCreateDeleteParticlesInOutflow(MpsParticleSystem *PSystem, MpsParticle *Particles) {
	
	//double minDistIOcreation = 1.1*PSystem->partDist;
	// Auxiliar variables
	int numCreatedParticlesAux, numDeletedParticlesAux, idLastRealandIOParticle, atLeastOneRealParticleCreated;
	// int idLastRealParticle;
//#pragma omp parallel
//{
	// Set auxiliar variables
	numCreatedParticlesAux = 0;						///< Number of created particles
	numDeletedParticlesAux = 0;						///< Number of deleted particles
	atLeastOneRealParticleCreated = 0;				///< If > 0, at least one real particle was created
	// idLastRealParticle = Particles->numParticles;	///< Array ID of the last real particle + 1
	idLastRealandIOParticle = Particles->numRealAndIOParticles;	///< Array ID of the last real+created particle + 1

	// printf("\n Ioid: %2d - Time: %.5e sec - NumPart: %4d - idLaRePar: %4d - numReIOPart: %4d - idLaReIOPar: %4d"
	// 	" - numCrPart: %3d - numDelPart: %3d ", 
	// 	Pio.ID, 0.0, Particles->numParticles, idLastRealParticle, Particles->numRealAndIOParticles, idLastRealandIOParticle,
	// 	numCreatedParticlesAux, numDeletedParticlesAux);

	// Set number of created particles to zero
	// Particles->numCreatedParticles = 0;

	// int sumAux, sumAux1, sumAux2;
	// sumAux = sumAux1 = sumAux2 = 0;

	struct timeval time_now, time_end;
	double time_taken;

// 	// start timer.
// 	gettimeofday(&time_now, NULL);
// //#pragma omp parallel for reduction(+:sumAux, sumAux1, sumAux2)
// // #pragma omp parallel for
// 	for(int ip=0; ip<Particles->numParticles; ip++) {
// 	int i = Particles->particleID[ip];
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

#pragma omp parallel for reduction(+:numCreatedParticlesAux, numDeletedParticlesAux, atLeastOneRealParticleCreated)
	for(int ip=0; ip<Particles->numParticles; ip++) {
	int i = Particles->particleID[ip];
	if(Particles->particleType[i] == PSystem->fluid) {
		
		// Computes euclidean signed distance between particle i and the current plane
		double signDistAux = calcSignedDistance(Particles->pos[i*3], Particles->pos[i*3+1], Particles->pos[i*3+2]);

// To prevent race condition, the access to the shared (global) variable "idLastRealandIOParticle" must be synchronized.
//#pragma omp critical (before_if)
//{
		// Verify if the signed distance of particle i from the IOplan is positive
		if(signDistAux > 0.0) {
			// Particle inside the fluid domain
			
			// Verify if the signed distance of particle i from the IOplan is shorter than partDist
			if(signDistAux < PSystem->partDist) {
			// if(signDistAux < PSystem->distCollisionLimit) {
				// Particle near the IOplan and its position shorter than partDist

				// Assign true to particle is in the IOregion
				Particles->isInIORegion[i] = true;

// Specifies that code is only executed on one thread at a time
//#pragma omp critical (after_if)
//{
				// Create two IO particles along the line normal to the IOplan and passing from the i-th position
				// Each thread of the "parallel for" has private (local) variables of two IO particles index
				// Be aware that they are defined using the shared (global) variable idLastRealandIOParticle
				int ipIO1, ipIO2, idIO1, idIO2;

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
					ipIO1 = idLastRealandIOParticle;
					idLastRealandIOParticle += 2;
				}
				ipIO2 = ipIO1+1;
				idIO1 = Particles->particleID[ipIO1];
				idIO2 = Particles->particleID[ipIO2];

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
				Particles->particleBC[idIO2] = PSystem->other;
				Particles->particleType[idIO1] = PSystem->inOutflowParticle;
				Particles->particleType[idIO2] = PSystem->inOutflowParticle;
				Particles->signDist[idIO1] = -signDistAux; // signDistAux is positive
				Particles->signDist[idIO2] = -signDistAux - PSystem->partDist;

				Particles->motherID[idIO1] = i;
				Particles->motherID[idIO2] = i;
				// Particles->press[idIO1] = Pio.press;
				// Particles->press[idIO2] = Pio.press;
				// w: ioPlan wall
				// Pio = Pw * di-io / di-w - Pi * (di-io - di-w) / di-w
				// Particles->press[idIO1] = Pio.press * 2.0 - Particles->press[i];
				// Particles->press[idIO2] = (Pio.press * (2.0 * signDistAux + PSystem->partDist) - Particles->press[i] * (signDistAux + PSystem->partDist) ) / (signDistAux + 0.01 * PSystem->reS2);
				
				Particles->ioflowID[idIO1] = Pio.ID;
				Particles->ioflowID[idIO2] = Pio.ID;

				Particles->pndi[idIO1] = Particles->pndi[i];
				Particles->pndi[idIO2] = Particles->pndi[i];
				Particles->numNeigh[idIO1] = Particles->numNeigh[i];
				Particles->numNeigh[idIO2] = Particles->numNeigh[i];
				Particles->pressAverage[idIO1] = 0.0;
				Particles->pressAverage[idIO2] = 0.0;
				Particles->pndki[idIO1] = Particles->pndki[i];
				Particles->pndki[idIO2] = Particles->pndki[i];
				Particles->pndski[idIO1] = Particles->pndski[i];
				Particles->pndski[idIO2] = Particles->pndski[i];
				Particles->pndSmall[idIO1] = Particles->pndSmall[i];
				Particles->pndSmall[idIO2] = Particles->pndSmall[i];

				// Update lastParticle ID
// The atomic construct allows multiple threads to safely update a shared (global) variable
//#pragma omp atomic
				//idLastRealandIOParticle++;

				// Update shared (global) variable number of IO particles
				numCreatedParticlesAux += 2;

//} //end of omp critical (after_if)
			}
			// Verify if the signed distance of particle i from the IOplan is shorter than reS
			else if(signDistAux < PSystem->reS) {
			// else if(signDistAux < 0.9 * 2.0 * PSystem->partDist) {
			// else if(signDistAux < 2.0 * PSystem->distCollisionLimit) {
				// Particle near the IOplan and its position ranges between partDist and reS

				// Verify if the particle was in the IOregion in the previous step
				if(Particles->isInIORegion[i] == true) {

					// Assign false to particle is in the IOregion
					Particles->isInIORegion[i] = false;

					// Create one Real particle and two IO particles along the line normal to the IOplan and passing from the idReal-th position
					// Each thread of the "parallel for" has private (local) variables of two IO particles index
					// Be aware that they are defined using the shared (global) variable idLastRealandIOParticle
					int ipReal, ipIO1, ipIO2, idReal, idIO1, idIO2;

// idReal atomically captures the original value of idLastRealandIOParticle. Then, atomically updates idLastRealandIOParticle
#pragma omp atomic capture
					{
						ipReal = idLastRealandIOParticle;
						idLastRealandIOParticle += 3;
					}
					ipIO1 = ipReal+1;
					ipIO2 = ipIO1+1;

					idReal = Particles->particleID[ipReal];
					idIO1 = Particles->particleID[ipIO1];
					idIO2 = Particles->particleID[ipIO2];
					
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
					Particles->particleBC[idReal] = Particles->particleBC[i];
					Particles->particleBC[idIO1] = PSystem->other;
					Particles->particleBC[idIO2] = PSystem->other;
					
					Particles->particleType[idReal] = PSystem->fluid;
					Particles->particleType[idIO1] = PSystem->inOutflowParticle;
					Particles->particleType[idIO2] = PSystem->inOutflowParticle;
					
					Particles->motherID[idIO1] = idReal;
					Particles->motherID[idIO2] = idReal;

					Particles->press[idReal] = Particles->press[i];
					// Particles->press[idIO1] = Pio.press;
					// Particles->press[idIO2] = Pio.press;

					Particles->ioflowID[idReal] = Pio.ID;
					Particles->ioflowID[idIO1] = Pio.ID;
					Particles->ioflowID[idIO2] = Pio.ID;
					
					Particles->isInIORegion[idReal] = true;
					Particles->isInIORegion[idIO1] = false;
					Particles->isInIORegion[idIO2] = false;

					Particles->signDist[idReal] = signDistAux - PSystem->partDist; // signDistAux is positive
					Particles->signDist[idIO1] = -signDistAux;
					Particles->signDist[idIO2] = -signDistAux - PSystem->partDist;

					Particles->pndi[idReal] = Particles->pndi[i];
					Particles->pndi[idIO1] = Particles->pndi[i];
					Particles->pndi[idIO2] = Particles->pndi[i];
					
					Particles->numNeigh[idReal] = Particles->numNeigh[i];
					Particles->pressAverage[idReal] = Particles->pressAverage[i];
					Particles->pndki[idReal] = Particles->pndki[i];
					Particles->pndski[idReal] = Particles->pndski[i];
					Particles->pndSmall[idReal] = Particles->pndSmall[i];
					Particles->numNeigh[idIO1] = Particles->numNeigh[i];
					Particles->numNeigh[idIO2] = Particles->numNeigh[i];
					Particles->pressAverage[idIO1] = 0.0;
					Particles->pressAverage[idIO2] = 0.0;
					Particles->pndki[idIO1] = Particles->pndki[i];
					Particles->pndki[idIO2] = Particles->pndki[i];
					Particles->pndski[idIO1] = Particles->pndski[i];
					Particles->pndski[idIO2] = Particles->pndski[i];
					Particles->pndSmall[idIO1] = Particles->pndSmall[i];
					Particles->pndSmall[idIO2] = Particles->pndSmall[i];

					// Update lastParticle ID
// The atomic construct allows multiple threads to safely update a shared (global) variable
//#pragma omp atomic
					//idLastRealandIOParticle++;

					// Update shared (global) variable number of created particles
					numCreatedParticlesAux += 3;

					// Update shared (global) variable that indicates that at least one real particle was created
					atLeastOneRealParticleCreated++;
				}
				else {
					// Computes ui . normal
					double ui_n = Particles->vel[i*3] * Pio.normal[0] + Particles->vel[i*3+1] * Pio.normal[1] + Particles->vel[i*3+2] * Pio.normal[2];

					// Verify the velocity component normal to the IOplan
					if (ui_n > 0.0) {
						// The velocity component normal to the IOplan is positive, then apply Inflow condition

						// Check if the release of a new particle is required
						// Conical scan region is applied
						
						// The direction of the conical region axis is the sum of the vectors rij pointing towards i-th particle

						// Effective neighboring j-th particles are verified if is outside the scan region
					}
					else {
						// The velocity component normal to the IOplan is negative, then apply Outflow condition
						// No new effective particle has to be generated

					}
				}
			}

			// Gets minimum absolute value
			if(signDistAux*signDistAux < Particles->signDist[i]*Particles->signDist[i])
				Particles->signDist[i] = signDistAux;

		}
		else {
			// The signed distance is negative, then the Particle is outside the fluid domain
			// Set particle data as Ghost data
			Particles->setParticleDataToGhost(i, PSystem);
			// printf("\nSAINDO Ioid: %d - PartID: %d", Pio.ID, i);

// 			// ID of "Last Real" and "Last Real+IO" particles
// 			int idReal, idRealIO;

// // Specifies that code is only executed on one thread at a time
// // #pragma omp critical (after_else)
// {
// 				idLastRealParticle--;
// 				idReal = idLastRealParticle;
// 				idLastRealandIOParticle--;
// 				idRealIO = idLastRealandIOParticle;
// } //end of omp critical (after_else)
			
// // idReal atomically captures the updated value of idLastRealandIOParticle.
// // #pragma omp atomic capture
// // 			{
// // 				idLastRealParticle -= 1;
// // 				idReal = idLastRealParticle;
// // 			}

// 			// Swap the data between "Last Real Particle" and i-th Ghost particle
// 			Particles->swapDataLastRealPartAndGhostPart(i, idReal, PSystem);

// // idRealIO atomically captures the updated value of idLastRealandIOParticle.
// // #pragma omp atomic capture
// // 			{
// // 				idLastRealandIOParticle -= 1;
// // 				idRealIO = idLastRealandIOParticle;
// // 			}

// 			// Swap the data between "Last Real+IO Particle" and "Last Real Particle" previously assigned as ghost
// 			Particles->swapDataLastRealIOPartAndLastRealPart(idReal, idRealIO, PSystem);

// // } //end of omp critical (after_else)
			
			// Update shared (global) variable number of deleted particles
			numDeletedParticlesAux++;
		}
	
//} //end of omp critical (before_if)

	} // end of if(Particles->particleType[i]...
	} // end of for(int i=0;...

//} //end of omp parallel

	// if(numCreatedParticlesAux != 0 && numCreatedParticlesAux != 200)
		// printf("\n NumParticles: %d", numCreatedParticlesAux);
	
	
	// Updates number of created particles
	Particles->numCreatedParticles += numCreatedParticlesAux;
	// Updates number of deleted particles
	Particles->numDeletedParticles += numDeletedParticlesAux;
	// Updates number of Real particles
	// Particles->numParticles -= numDeletedParticlesAux;
	// Updates number of particles including created particles
	// Particles->numRealAndIOParticles += numCreatedParticlesAux - numDeletedParticlesAux;
	Particles->numRealAndIOParticles += numCreatedParticlesAux;
	// Updates number of inOutflow (IO) particles considering no Real particles were created
	Particles->numIOParticles += numCreatedParticlesAux;

	// Assign realParticleCreated to false, i.e., no real particles was created
	Particles->realParticleCreated = false;
	if(atLeastOneRealParticleCreated > 0) {
		// Assign realParticleCreated to true, i.e., the number of created real particles is > 0
		Particles->realParticleCreated = true;
	}

	// stop timer.
	gettimeofday(&time_end, NULL);
	// Calculating total time taken by the program.
	time_taken = (time_end.tv_sec - time_now.tv_sec) * 1e6;
	time_taken = (time_taken + (time_end.tv_usec - time_now.tv_usec)) * 1e-6;
	// printf("\n Ioid: %2d - Time: %.5e sec - NumPart: %4d - idLaRePar: %4d - numReIOPart: %4d - idLaReIOPar: %4d"
	// 	" - numCrPart: %3d - numDelPart: %3d ", 
	// 	Pio.ID, time_taken, Particles->numParticles, idLastRealParticle, Particles->numRealAndIOParticles, idLastRealandIOParticle,
	// 	numCreatedParticlesAux, numDeletedParticlesAux);

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}


// Swap the data between Real particles in the array Particles->numRealAndIOParticles and the array Particles->numParticles
void MpsInflowOutflow::swapIdRealAndIOParticlesInOutflow(MpsParticleSystem *PSystem, MpsParticle *Particles) {

	// Auxiliar variables
	int numCreatedRealParticlesAux, idLastRealParticle;

	// Set auxiliar variables considering that all created particles are IO particles, i.e.,
	// the number of created real particles is zero
	numCreatedRealParticlesAux = 0;								///< Number of created Real particles
	idLastRealParticle = Particles->numParticles;				///< Array ID of the last real particle + 1

	struct timeval time_now, time_end;
	double time_taken;

	// start timer.
	gettimeofday(&time_now, NULL);

// Be aware with "race condition" in OpenMP
// Loop only for the created particles
// #pragma omp parallel for reduction(+:numCreatedRealParticlesAux)
	// for(int i=0; i<Particles->numCreatedParticles; i++) {
	// 	int idIO = i + Particles->numParticles;
	for(int idIOp=Particles->numParticles; idIOp<Particles->numRealAndIOParticles; idIOp++) {
		int idIO = Particles->particleID[idIOp];
	// Verify if the created particle is real
		// if(Particles->particleBC[idIO] == PSystem->inner) {
		if(Particles->particleType[idIO] == PSystem->fluid) {
		
			// Array ID to be filled with a Real particle
			int ipReal;
			// int idReal;

// idReal atomically captures the original value of idLastRealParticle. Then, atomically updates idLastRealParticle
// #pragma omp atomic capture
			ipReal = idLastRealParticle++;
			// idReal = Particles->particleID[ipReal];

			// if(Particles->particleBC[idReal] != PSystem->inner)
			// if(Particles->particleType[idReal] != PSystem->fluid)
			// {
				// Swapping particle ID in the array
				swap(Particles->particleID[ipReal], Particles->particleID[idIOp]);

				// Update the variables in the array
				// Swapping elements of the new real particle and IO
				// swap(Particles->pos[idReal*3  ], Particles->pos[idIO*3  ]);
				// swap(Particles->pos[idReal*3+1], Particles->pos[idIO*3+1]);
				// swap(Particles->pos[idReal*3+2], Particles->pos[idIO*3+2]);
				// swap(Particles->vel[idReal*3  ], Particles->vel[idIO*3  ]);
				// swap(Particles->vel[idReal*3+1], Particles->vel[idIO*3+1]);
				// swap(Particles->vel[idReal*3+2], Particles->vel[idIO*3+2]);
				// swap(Particles->particleBC[idReal], Particles->particleBC[idIO]);
				// swap(Particles->particleType[idReal], Particles->particleType[idIO]);
				// swap(Particles->press[idReal], Particles->press[idIO]);

				// swap(Particles->numNeigh[idReal], Particles->numNeigh[idIO]);
				// swap(Particles->pressAverage[idReal], Particles->pressAverage[idIO]);
				// swap(Particles->pndi[idReal], Particles->pndi[idIO]);
				// swap(Particles->pndki[idReal], Particles->pndki[idIO]);
				// swap(Particles->pndski[idReal], Particles->pndski[idIO]);
				// swap(Particles->pndSmall[idReal], Particles->pndSmall[idIO]);
			// }
			
			// Update shared (global) variable number of Real, IO particles
			numCreatedRealParticlesAux++;
		}
	} // end for(int idIO=0;...


	// Updates number of Real particles
	Particles->numParticles += numCreatedRealParticlesAux;
	// Updates number of inOutflow (IO) particles
	// Particles->numIOParticles = Particles->numCreatedParticles - numCreatedRealParticlesAux;
	Particles->numIOParticles -= numCreatedRealParticlesAux;
	// Updates number of particles including IO particles
	// Particles->numRealAndIOParticles = Particles->numParticles + Particles->numIOParticles;
	

// 	if(PSystem->mpsType == calcPressType::IMPLICIT_PND || PSystem->mpsType == calcPressType::IMPLICIT_PND_DIVU) {
// 		// Resize vector of pressure
// 		Particles->pressurePPE.resize(Particles->numParticles); // Resizing a dynamic-size vector
// #pragma omp parallel for
// 		for(int ip=0; ip<Particles->numParticles; ip++) {
// 			int i = Particles->particleID[ip];
// 			Particles->pressurePPE(i) = Particles->press[i];
// 		}
// 	}

	
	// // Auxliar variables
	// int newNumRealParticles = Particles->numParticles;
	// int newNumRealIOParticles = Particles->numRealAndIOParticles;
	
	// // Move ghost particles to the last positons of the array
	// // Reverse for loop
	// // https://stackoverflow.com/questions/275994/whats-the-best-way-to-do-a-backwards-loop-in-c-c-c
	// for (int ip = Particles->numParticles; ip --> 0; )
	// {
	// 	int i = Particles->particleID[ip];
	// 	if (Particles->particleType[i] == PSystem->ghost)
	// 	{
	// 		// ID of last Real particle
	// 		int ipLastRealParticle = newNumRealParticles - 1;
	// 		int iLastRealParticle = Particles->particleID[ipLastRealParticle];

	// 		// Swap the data between "Last Real Particle" and i-th Ghost particle
	// 		Particles->swapDataLastRealPartAndGhostPart(i, iLastRealParticle, PSystem);

	// 		// Decrease number of Real particles
	// 		newNumRealParticles--;

	// 		// ID of last Real+IO particle
	// 		int ipLastRealIOParticle = newNumRealIOParticles - 1;
	// 		int iLastRealIOParticle = Particles->particleID[ipLastRealIOParticle];

	// 		// Swap the data between "Last Real+IO Particle" and "Last Real Particle" previously assigned as ghost
	// 		Particles->swapDataLastRealIOPartAndLastRealPart(iLastRealParticle, iLastRealIOParticle, PSystem);
			
	// 		// Decrease number of Real+IO particles
	// 		newNumRealIOParticles--;	
	// 	}
	// }

	// // Update number of Real and Real+IO particles
	// Particles->numParticles = newNumRealParticles;
	// Particles->numRealAndIOParticles = newNumRealIOParticles;




	// stop timer.
	gettimeofday(&time_end, NULL);
	// Calculating total time taken by the program.
	time_taken = (time_end.tv_sec - time_now.tv_sec) * 1e6;
	time_taken = (time_taken + (time_end.tv_usec - time_now.tv_usec)) * 1e-6;
	// printf("\n IOid: %2d - Time: %.5e sec - NumPart: %4d - NumIOPart: %4d - NumReIoPart: %4d - idLastRePar: %3d", 
	// 	Pio.ID, time_taken, Particles->numParticles, Particles->numIOParticles, Particles->numRealAndIOParticles, idLastRealParticle);

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Verify overlaped IO-IO or IO-Real particles, assign overlaped IO particles as Ghost and updates the current number of ghost particles
void MpsInflowOutflow::checkOverlapedIOParticles(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	
	// Auxliar variables
	// int newNumRealParticles = numParticles;
	// int newNumRealIOParticles = numRealAndIOParticles;

	// int numIOParticlesAux = Particles->numIOParticles;
	double limitDist2 = 0.6*PSystem->partDist;
	limitDist2 *= limitDist2;
	int numGhostParticlesAux = 0;

	// Set outside particles as ghost
#pragma omp parallel for reduction(+:numGhostParticlesAux)
	for(int ioID=0; ioID<Particles->numIOParticles; ioID++) {

		int ip = ioID + Particles->numParticles;	// Global ID in the array Real+IO particles
		int i = Particles->particleID[ip];
	
		if(Particles->particleType[i] == PSystem->inOutflowParticle) {


			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			// double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			// double ni = Particles->pndSmall[i];
			// double sum = 0.0;

			int ix, iy, iz;
			Buckets->bucketCoordinates(ix, iy, iz, posXi, posYi, posZi, PSystem);
			int minZ = (iz-1)*((int)(PSystem->dim-2.0)); int maxZ = (iz+1)*((int)(PSystem->dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*PSystem->numBucketsXY + jy*PSystem->numBucketsX + jx;
				int j = Particles->firstParticleInBucket[jb];
				if(j == -1) continue;
				double plx, ply, plz;
				Particles->getPeriodicLengths(jb, plx, ply, plz, PSystem);
				while(true) {
					double v0ij, v1ij, v2ij, dstij2;
					// double v0imj, v1imj, v2imj, dstimj2;
					
					// Particle square distance r_ij^2 = (Xj - Xi_temporary_position)^2
					Particles->sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2, plx, ply, plz);
					// Mirror particle square distance r_imj^2 = (Xj - Xim_temporary_position)^2
					// Particles->sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2, plx, ply, plz);

					// If j and i is overlaped 
					if(dstij2 < limitDist2) {
					if(j != i) {

						// if(Particles->particleType[j] == PSystem->fluid) {
						// 	// Set particle data as Ghost data
						// 	Particles->setParticleDataToGhost(j, PSystem);
						// }
						// else {
						// 	// Set particle data as Ghost data
						// 	Particles->setParticleDataToGhost(i, PSystem);
						// }
						
						// Set particle data as Ghost data
						Particles->setParticleDataToGhost(i, PSystem);

						numGhostParticlesAux++;
					}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
		}
	}

	// Updates the number of ghost particles in the current step
	Particles->numGhostParticles += numGhostParticlesAux;

}
/////////////////////////////////////////////////
// Check IO particles in the Inflow/Outflow region
// void MpsInflowOutflow::checkIOParticlesInOutflow(MpsParticleSystem *PSystem, MpsParticle *Particles) {

// 	// Set number of inOutflow (IO) particles to Zero
// 	//Particles->numIOParticles = 0;
// 	// Update number of particles including IO particles
// 	//Particles->numRealAndIOParticles = Particles->numParticles + Particles->numIOParticles;

// 	// Auxiliar variables
// 	// int numRealParticlesAux, numIOParticlesAux, idLastRealParticle, idLastRealandIOParticle;
// 	int numRealParticlesAux, idLastRealParticle, idLastRealandIOParticle;

// 	// Set auxiliar variables
// 	numRealParticlesAux = 0;						///< Number of Real particles
// 	// numIOParticlesAux = Particles->numIOParticles;	///< Number of inOutflow (IO) particles
// 	idLastRealParticle = Particles->numParticles;	///< Array ID of the last real particle + 1
// 	idLastRealandIOParticle = Particles->numRealAndIOParticles;	///< Array ID of the last real+IO particle + 1

// 	struct timeval time_now, time_end;
// 	double time_taken;

// 	// start timer.
// 	gettimeofday(&time_now, NULL);

// // Be aware with "race condition" in OpenMP
// // Loop only for IO particles
// //#pragma omp parallel for reduction(+:numIOParticlesAux)
// 	for(int idIOp=Particles->numParticles; idIOp<Particles->numRealAndIOParticles; idIOp++) {
// 		int idIO = Particles->particleID[idIOp];
// 	// if(Particles->particleType[idIO] == PSystem->fluid) {
		
// 		// Computes euclidean signed distance between particle i and the current plane
// 		Particles->signDist[idIO] = calcSignedDistance(Particles->pos[idIO*3], Particles->pos[idIO*3+1], Particles->pos[idIO*3+2]);

// // To prevent race condition, the access to the shared (global) variable "idLastRealParticle" must be synchronized.
// //#pragma omp critical (before_if)
// //{
// 		// Verify if the signed distance of IO particle from the IOplan is positive
// 		// if(Particles->signDist[idIO] >= PSystem->epsilonZero) {
// 		if(Particles->signDist[idIO] >= 0) {
// 			// IO particle inside the fluid domain

// 			// Array ID to be filled with a Real particle
// 			int idReal;

// // idReal atomically captures the original value of idLastRealParticle. Then, atomically updates idLastRealParticle
// // #pragma omp atomic capture
// 			{
// 				idReal = idLastRealParticle;
// 				idLastRealParticle += 1;
// 			}
			
// 			// Update the variables in the array
// 			// Swapping elements of the new real particle and IO
// 			swap(Particles->pos[idReal*3  ], Particles->pos[idIO*3  ]);
// 			swap(Particles->pos[idReal*3+1], Particles->pos[idIO*3+1]);
// 			swap(Particles->pos[idReal*3+2], Particles->pos[idIO*3+2]);
// 			swap(Particles->vel[idReal*3  ], Particles->vel[idIO*3  ]);
// 			swap(Particles->vel[idReal*3+1], Particles->vel[idIO*3+1]);
// 			swap(Particles->vel[idReal*3+2], Particles->vel[idIO*3+2]);
// 			Particles->particleBC[idReal] = PSystem->inner;
// 			Particles->particleType[idReal] = PSystem->fluid;
// 			swap(Particles->press[idReal], Particles->press[idIO]);

// 			// Create one IO particle along the line normal to the IOplan and passing from the new idReal-th particle's position
// 			// Each thread of the "parallel for" has private (local) variables of IO particle index
// 			// Be aware that it is defined using the shared (global) variable idLastRealandIOParticle
// 			int idIO2;

// // idIO2 atomically captures the original value of idLastRealandIOParticle. Then, atomically updates idLastRealandIOParticle
// // #pragma omp atomic capture
// 			{
// 				idIO2 = idLastRealandIOParticle;
// 				idLastRealandIOParticle += 1;
// 			}

// 			// IO particle is placed at distance partDist from idReal-th position
// 			Particles->pos[idIO2*3  ] = Particles->pos[idReal*3  ] - (2.0*PSystem->partDist)*Pio.a;
// 			Particles->pos[idIO2*3+1] = Particles->pos[idReal*3+1] - (2.0*PSystem->partDist)*Pio.b;
// 			Particles->pos[idIO2*3+2] = Particles->pos[idReal*3+2] - (2.0*PSystem->partDist)*Pio.c;
			
// 			// Assign velocities equal to those of idReal-th particle
// 			Particles->vel[idIO2*3  ] = Particles->vel[idReal*3  ];
// 			Particles->vel[idIO2*3+1] = Particles->vel[idReal*3+1];
// 			Particles->vel[idIO2*3+2] = Particles->vel[idReal*3+2];

// 			// Set some variables
// 			Particles->particleBC[idIO2] = PSystem->other;
// 			Particles->particleType[idIO2] = PSystem->inOutflowParticle;
// 			Particles->press[idIO2] = Pio.press;

// 			// Update shared (global) variable number of Real, IO particles
// 			numRealParticlesAux += 1;
// 			// numIOParticlesAux += 1;
			
// 		}
	
// //} //end of omp critical (before_if)

// 	// }
// 	} // end for(int i=0;...

// //} //end of omp parallel

// 	// Updates number of Real particles
// 	Particles->numParticles += numRealParticlesAux;
// 	// Updates number of inOutflow (IO) particles
// 	// Particles->numIOParticles = numIOParticlesAux;
// 	// Updates number of particles including IO particles
// 	Particles->numRealAndIOParticles = Particles->numParticles + Particles->numIOParticles;

// 	// stop timer.
// 	gettimeofday(&time_end, NULL);
// 	// Calculating total time taken by the program.
// 	time_taken = (time_end.tv_sec - time_now.tv_sec) * 1e6;
// 	time_taken = (time_taken + (time_end.tv_usec - time_now.tv_usec)) * 1e-6;
// 	// printf("\n Time: %.5e sec - NumParticles: %d - idLastRealParticle: %d", time_taken, numIOParticlesAux, idLastRealParticle);
// 	printf("\n Time: %.5e sec - NumPart: %d - numRealAndIOPart: %d - idLastRealParticle: %d", time_taken, Particles->numParticles, Particles->numRealAndIOParticles, idLastRealParticle);
// }


// // Clear variables only of the inOutflow (IO) particles
// void MpsInflowOutflow::clearVariablesInOutflowParticles(MpsParticleSystem *PSystem, MpsParticle *Particles) {
// //#pragma omp parallel for schedule(dynamic,64)
// 	for(int ip=Particles->numParticles-1; ip<Particles->numRealAndIOParticles; ip++) {
// 	int i = Particles->particleID[ip];
// 	// Update some data of IO Particle
// 	Particles->particleType[i]=PSystem->ghost;
// 	Particles->particleBC[i]=PSystem->other;
// 	Particles->particleNearWall[i]=false;
// 	Particles->nearMeshType[i]=meshType::FIXED;
// 	Particles->distParticleWall2[i]=PSystem->nearInfinity;
// 	// Set zero to velocity and press of lastParticle
// 	for (int j = 0; j < 3; j++){
// 		//pos[i*3+j]=0.0;
// 		//Posk[i*3+j]=0.0;
// 		Particles->vel[i*3+j] = 0.0;
// 	}
// 	Particles->press[i] = 0.0;
// 	// Set maximum position to lastParticle
// 	Particles->pos[i*3  ] = PSystem->domainMaxX - PSystem->partDist;
// 	Particles->pos[i*3+1] = PSystem->domainMaxY - PSystem->partDist;
// 	if (PSystem->dim == 2)
// 		Particles->pos[i*3+2] = 0.0;
// 	else
// 		Particles->pos[i*3+2] = PSystem->domainMaxZ - PSystem->partDist;
	
// 	}
// }