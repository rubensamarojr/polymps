// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <experimental/filesystem> 	///< numeric_limits
#include <iostream>					///< cout
#include "MpsParticleCollision.h"

using namespace std;

// Constructor declaration
MpsParticleCollision::MpsParticleCollision()
{
}
// Destructor declaration
MpsParticleCollision::~MpsParticleCollision()
{
}

// Check collisions between particles
// Step-by-step improvement of MPS method in simulating violent free-surface motions and impact-loads
// https://doi.org/10.1016/j.cma.2010.12.001
void MpsParticleCollision::checkParticleCollisions(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
			// double mi = Particles->Dns[partType::PSystem->fluid];
			double mi;
			if(Particles->PTYPE[i] == 1) mi = PSystem->DNS_FL1;
			else mi = PSystem->DNS_FL2;
			double posXi = Particles->pos[i*3  ];double posYi = Particles->pos[i*3+1];double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];double velYi = Particles->vel[i*3+1];double velZi = Particles->vel[i*3+2];
			double dVelXi = 0.0;double dVelYi = 0.0;double dVelZi = 0.0;
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			
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
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle square distance r_ij^2 = (Xj - Xi_temporary_position)^2
					Particles->sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2, plx, ply, plz);
					// Mirror particle square distance r_imj^2 = (Xj - Xim_temporary_position)^2
					Particles->sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2, plx, ply, plz);

					// If j is inside the neighborhood of i and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < PSystem->distCollisionLimit2 && (dstij2 < dstimj2 || PSystem->wallType == boundaryWallType::PARTICLE)) {
						if(j != i) {
							double fDT = (velXi-Particles->vel[j*3  ])*v0ij+(velYi-Particles->vel[j*3+1])*v1ij+(velZi-Particles->vel[j*3+2])*v2ij;
							if(fDT > 0.0) {
								double mj;
								if(Particles->particleType[j]==PSystem->fluid)
								{
									if(Particles->PTYPE[j] == 1) mj = PSystem->DNS_FL1;
									else mj = PSystem->DNS_FL2;
								}
								else
								{
									mj = Particles->Dns[partType::WALL];
								}

								fDT *= PSystem->restitutionCollision*mj/(mi+mj)/dstij2;
								if(Particles->particleType[j]==PSystem->fluid)
								{
									dVelXi -= v0ij*fDT;		dVelYi -= v1ij*fDT;		dVelZi -= v2ij*fDT;
								}
								else
								{
									dVelXi -= 2.0*v0ij*fDT;	dVelYi -= 2.0*v1ij*fDT;	dVelZi -= 2.0*v2ij*fDT;
								}
							}
							/*
							double fDT = (Particles->vel[j*3  ]-velXi)*v0+(Particles->vel[j*3+1]-velYi)*v1+(Particles->vel[j*3+2]-velZi)*v2;
							double mj;
							if(Particles->particleType[j]==PSystem->fluid)
								mj = Particles->Dns[partType::PSystem->fluid];
							else
								mj = Particles->Dns[partType::WALL];
							fDT *= PSystem->restitutionCollision*mj/(mi+mj)/dst2;
							velXi2 += v0*fDT;		vecYi2 += v1*fDT;		vecZi2 += v2*fDT;
							*/
						}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}

			Particles->dvelCollision[i*3  ]=dVelXi;	Particles->dvelCollision[i*3+1]=dVelYi;	Particles->dvelCollision[i*3+2]=dVelZi;
			//accStar[i*3  ]=Particles->vel[i*3  ]+dVelXi;	accStar[i*3+1]=Particles->vel[i*3+1]+dVelYi;	accStar[i*3+2]=Particles->vel[i*3+2]+dVelZi;
		}
	}
#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
			// CHANGED !!!
			//Particles->pos[i*3  ]+=(acc[i*3  ]-Particles->vel[i*3  ])*PSystem->timeStep; Particles->pos[i*3+1]+=(acc[i*3+1]-Particles->vel[i*3+1])*PSystem->timeStep; Particles->pos[i*3+2]+=(acc[i*3+2]-Particles->vel[i*3+2])*PSystem->timeStep;
			Particles->vel[i*3  ]+=Particles->dvelCollision[i*3  ];	Particles->vel[i*3+1]+=Particles->dvelCollision[i*3+1];	Particles->vel[i*3+2]+=Particles->dvelCollision[i*3+2];
			
			//Velk[i*3  ]=Particles->vel[i*3  ];	Velk[i*3+1]=Particles->vel[i*3+1];	Velk[i*3+2]=Particles->vel[i*3+2];
			//Particles->pos[i*3  ]=Posk[i*3  ]+Particles->vel[i*3  ]*PSystem->timeStep; Particles->pos[i*3+1]=Posk[i*3+1]+Particles->vel[i*3+1]*PSystem->timeStep; Particles->pos[i*3+2]=Posk[i*3+2]+Particles->vel[i*3+2]*PSystem->timeStep;
		}
		Particles->dvelCollision[i*3  ]=0.0;	Particles->dvelCollision[i*3+1]=0.0;	Particles->dvelCollision[i*3+2]=0.0;
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Check collisions between particles (Dynamic Particle Collision)
// Enhanced weakly-compressible MPS method for violent free-surface flows: Role of particle regularization techniques
// https://doi.org/10.1016/j.jcp.2021.110202
void MpsParticleCollision::checkDynamicParticleCollisions(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	
	double Wij5 = 0.5*0.5*0.5*0.5*3.0;
	//double Wij5 = 0.5*0.5;
	double global_pmax = 0.0;
	double gravityMod = sqrt(PSystem->gravityX*PSystem->gravityX + PSystem->gravityY*PSystem->gravityY + PSystem->gravityZ*PSystem->gravityZ);
	
	// Compute maximum pressure on the walls
#pragma omp parallel
{
	double local_pmax = 0.0;
#pragma omp for
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->wall) {
			local_pmax = max(local_pmax, Particles->press[i]);
		}
	}
#pragma omp critical
	{
		if (local_pmax > global_pmax)
			global_pmax = local_pmax;
	}
}
	// Compute collision and repulsive terms and the dynamic coefficients
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
	//		double mi = Particles->Dns[partType::PSystem->fluid];
			double mi;
			if(Particles->PTYPE[i] == 1) mi = PSystem->DNS_FL1;
			else mi = PSystem->DNS_FL2;

			double posXi = Particles->pos[i*3  ];double posYi = Particles->pos[i*3+1];double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];double velYi = Particles->vel[i*3+1];double velZi = Particles->vel[i*3+2];
			double dVelXi = 0.0;double dVelYi = 0.0;double dVelZi = 0.0;
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];

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
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle square distance r_ij^2 = (Xj - Xi_temporary_position)^2
					Particles->sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2, plx, ply, plz);
					// Mirror particle square distance r_imj^2 = (Xj - Xim_temporary_position)^2
					Particles->sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2, plx, ply, plz);

					// If j is inside the neighborhood of i and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < PSystem->partDist*PSystem->partDist && (dstij2 < dstimj2 || PSystem->wallType == boundaryWallType::PARTICLE)) {
						if(j != i) {
							double mj;
							if(Particles->particleType[j]==PSystem->fluid)
							{
								if(Particles->PTYPE[j] == 1) mj = PSystem->DNS_FL1;
								else mj = PSystem->DNS_FL2;
							}
							else
							{
								mj = Particles->Dns[partType::WALL];
							}
							double dst = sqrt(dstij2);
							
							// inter-particle distance
							double Wij = 0.0;
							if(dst > PSystem->epsilonZero)
							{
								double w1 = 1.0 - dst*PSystem->invPartDist;
								double w2 = 4.0*dst*PSystem->invPartDist + 1.0;
								Wij = w1*w1*w1*w1*w2;
								//Wij = w1*w1;
							}
							double chi = sqrt(Wij/Wij5);
							double kappa = 0.0;
							if(dst < 0.5*PSystem->partDist)
							{
								kappa = 1.0;
							}
							else
							{
								kappa = chi;
							}

							double fDT = (velXi-Particles->vel[j*3  ])*v0ij+(velYi-Particles->vel[j*3+1])*v1ij+(velZi-Particles->vel[j*3+2])*v2ij;
							if(fDT > 0.0) {
								fDT *= kappa*2.0*mj/(mi+mj)/dstij2;
								//fDT *= PSystem->restitutionCollision*mj/(mi+mj)/dstij2;
								if(Particles->particleType[j]==PSystem->fluid)
								{
									dVelXi -= v0ij*fDT;		dVelYi -= v1ij*fDT;		dVelZi -= v2ij*fDT;
								}
								else //if(particleBC[i] == PSystem->surface)
								{
									dVelXi -= 2.0*v0ij*fDT;		dVelYi -= 2.0*v1ij*fDT;		dVelZi -= 2.0*v2ij*fDT;
								}
							}
							else
							{
								// Dynamic background pressure
								double pmax = 2.0/3.0*Particles->RHO[i]*gravityMod*0.2;
								double pmin = Particles->RHO[i]*gravityMod*PSystem->partDist;

								double ptil = max(min(PSystem->lambdaCollision*fabs(Particles->press[i]+Particles->press[j]), PSystem->lambdaCollision*pmax), pmin);
								double pb = ptil*chi;
								double rep = PSystem->timeStep/mi*chi*pb/dstij2;
								if(Particles->particleType[j]==PSystem->fluid)
								{
									dVelXi -= v0ij*rep;		dVelYi -= v1ij*rep;		dVelZi -= v2ij*rep;
								}
								else //if(particleBC[i] == PSystem->surface)
								{
									dVelXi -= 2.0*v0ij*rep;		dVelYi -= 2.0*v1ij*rep;		dVelZi -= 2.0*v2ij*rep;
								}
							}
						}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}

			Particles->dvelCollision[i*3  ]=dVelXi;	Particles->dvelCollision[i*3+1]=dVelYi;	Particles->dvelCollision[i*3+2]=dVelZi;
		}
	}
	// Update velocity and position
#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
			
			Particles->vel[i*3  ]+=Particles->dvelCollision[i*3  ];	Particles->vel[i*3+1]+=Particles->dvelCollision[i*3+1];	Particles->vel[i*3+2]+=Particles->dvelCollision[i*3+2];
			Particles->pos[i*3  ]+=Particles->dvelCollision[i*3  ]*PSystem->timeStep; Particles->pos[i*3+1]+=Particles->dvelCollision[i*3+1]*PSystem->timeStep; Particles->pos[i*3+2]+=Particles->dvelCollision[i*3+2]*PSystem->timeStep;
			/*
			double drNew[3], drMod, drMin, duNew[3];
			duNew[0] = Particles->dvelCollision[i*3  ];
			duNew[1] = Particles->dvelCollision[i*3+1];
			duNew[2] = Particles->dvelCollision[i*3+2];
			drNew[0] = Particles->dvelCollision[i*3  ]*PSystem->timeStep;
			drNew[1] = Particles->dvelCollision[i*3+1]*PSystem->timeStep;
			drNew[2] = Particles->dvelCollision[i*3+2]*PSystem->timeStep;
			drMod = (drNew[0]*drNew[0] + drNew[1]*drNew[1] + drNew[2]*drNew[2]);
			drMin = min(0.1*PSystem->partDist, drMod);
			if(drMin > PSystem->epsilonZero)
			{
				Particles->pos[i*3  ]+=drMin*drNew[0]/drMod;
				Particles->pos[i*3+1]+=drMin*drNew[1]/drMod;
				Particles->pos[i*3+2]+=drMin*drNew[2]/drMod;
				Particles->vel[i*3  ]+=drMin*duNew[0]/drMod;
				Particles->vel[i*3+1]+=drMin*duNew[1]/drMod;
				Particles->vel[i*3+2]+=drMin*duNew[2]/drMod;
			}
			*/
		}
		Particles->dvelCollision[i*3  ]=0.0;	Particles->dvelCollision[i*3+1]=0.0;	Particles->dvelCollision[i*3+2]=0.0;
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}