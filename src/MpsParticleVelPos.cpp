// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <experimental/filesystem> 	///< numeric_limits
#include <iostream>					///< cout
#include "MpsParticleVelPos.h"

using namespace std;

// Constructor declaration
MpsParticleVelPos::MpsParticleVelPos()
{
}
// Destructor declaration
MpsParticleVelPos::~MpsParticleVelPos()
{
}

// Prediction of particle velocity and position
void MpsParticleVelPos::predictionVelocityPosition(MpsParticleSystem *PSystem, MpsParticle *Particles) {
#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
			Particles->vel[i*3  ] += Particles->acc[i*3  ]*PSystem->timeStep;	Particles->vel[i*3+1] += Particles->acc[i*3+1]*PSystem->timeStep;	Particles->vel[i*3+2] += Particles->acc[i*3+2]*PSystem->timeStep;
			//if(Particles->particleType[i] == PSystem->fluid) {
			Particles->pos[i*3  ] += Particles->vel[i*3  ]*PSystem->timeStep;	Particles->pos[i*3+1] += Particles->vel[i*3+1]*PSystem->timeStep;	Particles->pos[i*3+2] += Particles->vel[i*3+2]*PSystem->timeStep;
			//}
		}
		Particles->acc[i*3]=Particles->acc[i*3+1]=Particles->acc[i*3+2]=0.0;
		Particles->dvelCollision[i*3]=Particles->dvelCollision[i*3+1]=Particles->dvelCollision[i*3+2]=0.0;
		Particles->wallParticleForce1[i*3]=Particles->wallParticleForce1[i*3+1]=Particles->wallParticleForce1[i*3+2]=0.0;
		Particles->wallParticleForce2[i*3]=Particles->wallParticleForce2[i*3+1]=Particles->wallParticleForce2[i*3+2]=0.0;
		Particles->npcdDeviation[i*3]=Particles->npcdDeviation[i*3+1]=Particles->npcdDeviation[i*3+2]=0.0;
		Particles->numNeighWallContribution[i]=0;
		Particles->particleNearWall[i]=false;
		// Set squared distance of particle to triangle mesh to ~infinite
		Particles->distParticleWall2[i] = 10e8*PSystem->partDist;
		Particles->numNeighborsSurfaceParticles[i]=0.0;

		if(PSystem->wallType == boundaryWallType::POLYGON) {
			// Set mirrored particle to ~infinite if wall particles are used
			Particles->mirrorParticlePos[i*3  ] = 10e8*PSystem->partDist; Particles->mirrorParticlePos[i*3+1] = 10e8*PSystem->partDist; Particles->mirrorParticlePos[i*3+2] = 10e8*PSystem->partDist;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Correction of velocities and positions
void MpsParticleVelPos::correctionVelocityPosition(MpsParticleSystem *PSystem, MpsParticle *Particles) {
	PSystem->velMax = 0.0;						// Maximum flow velocity
	// double auxiliar[5] = {1.2, -3.3, 4.3, -0.3, 5.6};

	// https://stackoverflow.com/questions/39989473/use-openmp-in-c11-to-find-the-maximum-of-the-calculated-values
#pragma omp parallel
{
	double local_vMax = 0.0;
#pragma omp for
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
			Particles->vel[i*3  ]+=Particles->acc[i*3  ]*PSystem->timeStep;	Particles->vel[i*3+1]+=Particles->acc[i*3+1]*PSystem->timeStep;	Particles->vel[i*3+2]+=Particles->acc[i*3+2]*PSystem->timeStep;
			//if(Particles->particleType[i] == PSystem->fluid) {
			Particles->pos[i*3  ]+=Particles->acc[i*3  ]*PSystem->timeStep*PSystem->timeStep;	Particles->pos[i*3+1]+=Particles->acc[i*3+1]*PSystem->timeStep*PSystem->timeStep;	Particles->pos[i*3+2]+=Particles->acc[i*3+2]*PSystem->timeStep*PSystem->timeStep;
			//}
			Particles->acc[i*3]=Particles->acc[i*3+1]=Particles->acc[i*3+2]=0.0;

			//Particles->pos[i*3  ]=Particles->Posk[i*3 ]+Particles->vel[i*3  ]*PSystem->timeStep;	Particles->pos[i*3+1]=Particles->Posk[i*3+1]+Particles->vel[i*3+1]*PSystem->timeStep;	Particles->pos[i*3+2]=Particles->Posk[i*3+2]+Particles->vel[i*3+2]*PSystem->timeStep;
			//Particles->Posk[i*3  ]=Particles->pos[i*3  ];	Particles->Posk[i*3+1]=Particles->pos[i*3+1];	Particles->Posk[i*3+2]=Particles->pos[i*3+2];
			//Particles->Velk[i*3  ]=Particles->vel[i*3  ];	Particles->Velk[i*3+1]=Particles->vel[i*3+1];	Particles->Velk[i*3+2]=Particles->vel[i*3+2];

			//Particles->wallParticleForce1[i*3  ]=Particles->acc[i*3  ];	Particles->wallParticleForce1[i*3+1]=Particles->acc[i*3+1];	Particles->wallParticleForce1[i*3+2]=Particles->acc[i*3+2];
			//Particles->wallParticleForce2[i*3  ]=Particles->Acv[i*3  ];	Particles->wallParticleForce2[i*3+1]=Particles->Acv[i*3+1];	Particles->wallParticleForce2[i*3+2]=Particles->Acv[i*3+2];
			//Particles->acc[i*3]=Particles->acc[i*3+1]=Particles->acc[i*3+2]=0.0;
			//Particles->Acv[i*3]=Particles->Acv[i*3+1]=Particles->Acv[i*3+2]=0.0;

			double vMod2 = Particles->vel[i*3  ]*Particles->vel[i*3  ] + Particles->vel[i*3+1]*Particles->vel[i*3+1] + Particles->vel[i*3+2]*Particles->vel[i*3+2];
			if(vMod2 > local_vMax*local_vMax)
				local_vMax = sqrt(vMod2);
		}
	}

#pragma omp critical
	{
		if (local_vMax > PSystem->velMax)
			PSystem->velMax = local_vMax;
	}
}
	PSystem->CFLcurrent = PSystem->timeStep*PSystem->velMax/PSystem->partDist;

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Update velocity at wall and dummy particles
void MpsParticleVelPos::updateVelocityParticlesWallDummy(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	double velWallx = 0.0;
	double velWally = 0.0;
	double velWallz = 0.0;
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
	if(Particles->particleType[i] == PSystem->dummyWall) {
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double duXi = 0.0;	double duYi = 0.0;	double duZi = 0.0;
		double ni = 0.0;
		
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
				
				// Particle square distance r_ij^2 = (Xj - Xi_temporary_position)^2
				Particles->sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2, plx, ply, plz);

				if(dstij2 < PSystem->reL2) {
//				if(j != i) {
				if(j != i && Particles->particleType[j] == PSystem->fluid) {
					double dst = sqrt(dstij2);
					double wL = Particles->weight(dst, PSystem->reL, PSystem->weightType);
					ni += wL;
					duXi += Particles->vel[j*3  ]*wL;
					duYi += Particles->vel[j*3+1]*wL;
					duZi += Particles->vel[j*3+2]*wL;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		if(ni > PSystem->epsilonZero) {
			Particles->vel[i*3  ] = 2.0*velWallx - duXi/ni;
			Particles->vel[i*3+1] = 2.0*velWally - duYi/ni;
			Particles->vel[i*3+2] = 2.0*velWallz - duZi/ni;
		}
		else {
			Particles->vel[i*3  ] = 2.0*velWallx - duXi;
			Particles->vel[i*3+1] = 2.0*velWally - duYi;
			Particles->vel[i*3+2] = 2.0*velWallz - duZi;
		}
	}}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}
