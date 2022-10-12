// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <experimental/filesystem> 	///< numeric_limits
#include <iostream>					///< strings and c-strings
#include "MpsBucket.h"

using namespace std;

// Constructor declaration
MpsBucket::MpsBucket()
{
}
// Destructor declaration
MpsBucket::~MpsBucket()
{
}

// Allocation of memory for buckets data
void MpsBucket::allocateMemory(MpsParticleSystem *PSystem, MpsParticle *Particles) {
	
	// First guess of buckets values
	PSystem->bucketSide = PSystem->reL*(1.0+PSystem->cflNumber);				// Length of one bucket side
	PSystem->invBucketSide = 1.0/PSystem->bucketSide;

	PSystem->numBucketsX = (int)((PSystem->domainMaxX - PSystem->domainMinX)*PSystem->invBucketSide) + 3;		// Number of buckets in the x direction in the analysis domain
	PSystem->numBucketsY = (int)((PSystem->domainMaxY - PSystem->domainMinY)*PSystem->invBucketSide) + 3;		// Number of buckets in the y direction in the analysis domain
	PSystem->numBucketsZ = (int)((PSystem->domainMaxZ - PSystem->domainMinZ)*PSystem->invBucketSide) + 3;		// Number of buckets in the z direction in the analysis domain
	if((int)PSystem->dim == 2) {	PSystem->numBucketsZ = 1; }

	Particles->bucketTypeBC = 	 		(int*)malloc(sizeof(int) * PSystem->numBC);			// Type of Domain Boundary Condition in Bucket
	Particles->periodicDirection =	 	(bool*)malloc(sizeof(bool) * PSystem->numBC*3);		// Periodic direction in domain (x, y, or z)
	//Particles->periodicLength = 	 	(double*)malloc(sizeof(double) * PSystem->numBC*3);	// Periodic length in x, y and z direction
	for(int b=0; b<PSystem->numBC; b++){
		// domainTypeBC == 0: None
		// domainTypeBC == 1: Periodic
		Particles->bucketTypeBC[b] = PSystem->domainTypeBC;
		Particles->periodicDirection[b*3  ] = PSystem->periodicDirectionX;
		Particles->periodicDirection[b*3+1] = PSystem->periodicDirectionY;
		Particles->periodicDirection[b*3+2] = PSystem->periodicDirectionZ;
	}
	PSystem->periodicLength[0] = PSystem->periodicLength[1] = PSystem->periodicLength[2] = 0.0;		// Periodic length in x, y and z direction

	// Compute domain limits and adjust the buckets values
	if(PSystem->wallType == boundaryWallType::PARTICLE || PSystem->domainTypeBC == 1) {
		calcDomainLimits(PSystem, Particles);
		PSystem->invBucketSide = 1.0/PSystem->bucketSide;
	}
	PSystem->numBucketsXY = PSystem->numBucketsX*PSystem->numBucketsY;
	PSystem->numBucketsXYZ = PSystem->numBucketsX*PSystem->numBucketsY*PSystem->numBucketsZ;					// Number of buckets in analysis area
	
	std::cout << std::endl << "DomainMIN: " << PSystem->domainMinX << " " << PSystem->domainMinY << " " << PSystem->domainMinZ;
	std::cout << std::endl << "DomainMAX: " << PSystem->domainMaxX << " " << PSystem->domainMaxY << " " << PSystem->domainMaxZ;
	if(PSystem->domainTypeBC == 1) {
		std::cout << std::endl << "PhysicMIN: " << PSystem->physDomMinX << " " << PSystem->physDomMinY << " " << PSystem->physDomMinZ;
		std::cout << std::endl << "PhysicMAX: " << PSystem->physDomMaxX << " " << PSystem->physDomMaxY << " " << PSystem->physDomMaxZ;
	}
	std::cout << std::endl << "Num Buckt: " << PSystem->numBucketsX << " " << PSystem->numBucketsY << " " << PSystem->numBucketsZ;
	std::cout << std::endl << "BucktSide: " << PSystem->bucketSide << " 1/BucketSide: " << PSystem->invBucketSide;
	std::cout << std::endl << "PeriodicL: " << PSystem->periodicLength[0] << " " << PSystem->periodicLength[1] << " " << PSystem->periodicLength[2];
	std::cout << std::endl;

	Particles->firstParticleInBucket = 	(int*)malloc(sizeof(int) * PSystem->numBucketsXYZ);	// First particle number stored in the bucket
	Particles->lastParticleInBucket = 	(int*)malloc(sizeof(int) * PSystem->numBucketsXYZ);	// Last particle number stored in the bucket
	Particles->nextParticleInSameBucket  = (int*)malloc(sizeof(int) * Particles->numParticles);	// Next particle number in the same bucket
	Particles->bucketPeriodicBC =	 	(int*)malloc(sizeof(int) * PSystem->numBucketsXYZ);	// Periodic Boundary Condition of the bucket
}


// Compute domain limits
void MpsBucket::calcDomainLimits(MpsParticleSystem *PSystem, MpsParticle *Particles)
{
	double **limDom;
	limDom = new double *[3];
	for(int i=0; i<3; i++) limDom[i] = new double[3];

	// limitTypeBC = 0: Border particle positions
	// limitTypeBC = 1: Domain limits min and max
	if(PSystem->limitTypeBC == 0) {
		// wall_type = 0: Use border particle positions to define all domain limits
		// wall_type = 1: Use border particle positions to define periodic domain limits
		limDom[0][0] = limDom[0][1] = Particles->pos[0*3  ];
		limDom[1][0] = limDom[1][1] = Particles->pos[0*3+1];
		limDom[2][0] = limDom[2][1] = Particles->pos[0*3+2];
		
		for(int i=0; i<Particles->numParticles; i++) {
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			limDom[0][0] = min(limDom[0][0], posXi);
			limDom[0][1] = max(limDom[0][1], posXi);
			limDom[1][0] = min(limDom[1][0], posYi);
			limDom[1][1] = max(limDom[1][1], posYi);
			limDom[2][0] = min(limDom[2][0], posZi);
			limDom[2][1] = max(limDom[2][1], posZi);
		}
		int testdim = 0;
		if(limDom[0][0] != limDom[0][1])
			testdim++;
		if(limDom[1][0] != limDom[1][1])
			testdim++;
		if(limDom[2][0] != limDom[2][1])
			testdim++;
		if(testdim != PSystem->dim) {
			fprintf(stderr, "\n Dimensions in json [%d] and grid file [%d] do not match!\n\n", int(PSystem->dim), testdim);
			exit(10);
		}
	}
	else {
		// Adopt min and max values from json to define domain limits
		limDom[0][0] = PSystem->domainMinX; limDom[0][1] = PSystem->domainMaxX;
		limDom[1][0] = PSystem->domainMinY; limDom[1][1] = PSystem->domainMaxY;
		limDom[2][0] = PSystem->domainMinZ; limDom[2][1] = PSystem->domainMaxZ;
	}

	// domainTypeBC == 0: None
	// domainTypeBC == 1: Periodic
	if(PSystem->domainTypeBC == 0 && PSystem->wallType == boundaryWallType::PARTICLE) { // Whithout any special domain boundary condition
		
		// Shift half particle distance
		limDom[0][0] -= 0.5*PSystem->partDist;	limDom[0][1] += 0.5*PSystem->partDist;
		limDom[1][0] -= 0.5*PSystem->partDist;	limDom[1][1] += 0.5*PSystem->partDist;
		if(PSystem->dim==3) {
			limDom[2][0] -= 0.5*PSystem->partDist;	limDom[2][1] += 0.5*PSystem->partDist;
		}
		
		for(int k=0; k<PSystem->dim; k++) {
			limDom[k][0] -= PSystem->bucketSide;
			limDom[k][1] += PSystem->bucketSide;
			if(k == 0) {
				PSystem->numBucketsX = (long)((limDom[k][1]-limDom[k][0])/PSystem->bucketSide+1);
				limDom[k][1] = limDom[k][0] + PSystem->numBucketsX*PSystem->bucketSide; // Adjust the maximum limit o X
			}
			else if(k == 1) {
				PSystem->numBucketsY = (long)((limDom[k][1]-limDom[k][0])/PSystem->bucketSide+1);
				limDom[k][1] = limDom[k][0] + PSystem->numBucketsY*PSystem->bucketSide; // Adjust the maximum limit o Y
			}
			else if(k == 2) {
				PSystem->numBucketsZ = (long)((limDom[k][1]-limDom[k][0])/PSystem->bucketSide+1);
				limDom[k][1] = limDom[k][0] + PSystem->numBucketsZ*PSystem->bucketSide; // Adjust the maximum limit o Z
			}
		}
	}
	else if(PSystem->domainTypeBC == 1) {
		for(int b=0; b<PSystem->numBC; b++) {
			if(Particles->bucketTypeBC[b] == domainBC::PERIODIC) {
				bool periodicX =  Particles->periodicDirection[b*3] && !Particles->periodicDirection[b*3+1] && !Particles->periodicDirection[b*3+2];
				bool periodicY = !Particles->periodicDirection[b*3] &&  Particles->periodicDirection[b*3+1] && !Particles->periodicDirection[b*3+2];
				bool periodicZ = !Particles->periodicDirection[b*3] && !Particles->periodicDirection[b*3+1] &&  Particles->periodicDirection[b*3+2];
				// Periodic in X
				if(periodicX) {
					PSystem->periodicLength[0] = limDom[0][1] - limDom[0][0] + PSystem->partDist;
					limDom[0][0] -= PSystem->partDist*0.5;
					limDom[0][1] += PSystem->partDist*0.5;
					// Physical domain
					PSystem->physDomMinX = limDom[0][0];
					PSystem->physDomMaxX = limDom[0][1];
					PSystem->numBucketsX = (long)((limDom[0][1]-limDom[0][0])/PSystem->bucketSide);
					// Adjust bucketSide to perfectly divide the domain without leftovers
					PSystem->bucketSide = PSystem->periodicLength[0]/(long)(PSystem->periodicLength[0]/PSystem->bucketSide);
					// Analysis domain is extended from the boundaries by one bucket width
					limDom[0][0] -= PSystem->bucketSide;
					limDom[0][1] += PSystem->bucketSide;
					PSystem->numBucketsX += 2;
					// Adjust limits of Y
					if(PSystem->wallType == boundaryWallType::PARTICLE) {
						limDom[1][0] -= PSystem->bucketSide;
						limDom[1][1] += PSystem->bucketSide;
					}
					else {
						limDom[1][0] = PSystem->domainMinY - PSystem->bucketSide;
						limDom[1][1] = PSystem->domainMaxY + PSystem->bucketSide;
					}
					// Physical domain
					PSystem->physDomMinY = limDom[1][0];
					PSystem->physDomMaxY = limDom[1][1];
					PSystem->numBucketsY = (long)((limDom[1][1]-limDom[1][0])/PSystem->bucketSide+1);
					limDom[1][1] = limDom[1][0] + PSystem->numBucketsY*PSystem->bucketSide;
					if(PSystem->dim == 3) {
						// Adjust limits of Z
						if(PSystem->wallType == boundaryWallType::PARTICLE) {
							limDom[2][0] -= PSystem->bucketSide;
							limDom[2][1] += PSystem->bucketSide;
						}
						else {
							limDom[2][0] = PSystem->domainMinZ - PSystem->bucketSide;
							limDom[2][1] = PSystem->domainMaxZ + PSystem->bucketSide;
						}
						// Physical domain
						PSystem->physDomMinZ = limDom[2][0];
						PSystem->physDomMaxZ = limDom[2][1];
						PSystem->numBucketsZ = (long)((limDom[2][1]-limDom[2][0])/PSystem->bucketSide+1);
						limDom[2][1] = limDom[2][0] + PSystem->numBucketsZ*PSystem->bucketSide;
					}
				}
				// Periodic in Y
				if(periodicY) {
					PSystem->periodicLength[1] = limDom[1][1] - limDom[1][0] + PSystem->partDist;
					limDom[1][0] -= PSystem->partDist*0.5;
					limDom[1][1] += PSystem->partDist*0.5;
					// Physical domain
					PSystem->physDomMinY = limDom[1][0];
					PSystem->physDomMaxY = limDom[1][1];
					PSystem->numBucketsY = (long)((limDom[1][1]-limDom[1][0])/PSystem->bucketSide);
					// Adjust bucketSide to perfectly divide the domain without leftovers
					PSystem->bucketSide = PSystem->periodicLength[1]/(long)(PSystem->periodicLength[1]/PSystem->bucketSide);
					// Analysis domain is extended from the boundaries by one bucket width
					limDom[1][0] -= PSystem->bucketSide;
					limDom[1][1] += PSystem->bucketSide;
					PSystem->numBucketsY += 2;
					// Adjust limits of X
					if(PSystem->wallType == boundaryWallType::PARTICLE) {
						limDom[0][0] -= PSystem->bucketSide;
						limDom[0][1] += PSystem->bucketSide;
					}
					else {
						limDom[0][0] = PSystem->domainMinX - PSystem->bucketSide;
						limDom[0][1] = PSystem->domainMaxX + PSystem->bucketSide;
					}
					// Physical domain
					PSystem->physDomMinX = limDom[0][0];
					PSystem->physDomMaxX = limDom[0][1];
					PSystem->numBucketsX = (long)((limDom[0][1]-limDom[0][0])/PSystem->bucketSide+1);
					limDom[0][1] = limDom[0][0] + PSystem->numBucketsX*PSystem->bucketSide;
					if(PSystem->dim == 3) {
						// Adjust limits of Z
						if(PSystem->wallType == boundaryWallType::PARTICLE) {
							limDom[2][0] -= PSystem->bucketSide;
							limDom[2][1] += PSystem->bucketSide;
						}
						else {
							limDom[2][0] = PSystem->domainMinZ - PSystem->bucketSide;
							limDom[2][1] = PSystem->domainMaxZ + PSystem->bucketSide;
						}
						// Physical domain
						PSystem->physDomMinZ = limDom[2][0];
						PSystem->physDomMaxZ = limDom[2][1];
						PSystem->numBucketsZ = (long)((limDom[2][1]-limDom[2][0])/PSystem->bucketSide+1);
						limDom[2][1] = limDom[2][0] + PSystem->numBucketsZ*PSystem->bucketSide;
					}
				}
				// Periodic in Z
				if(periodicZ) {
					PSystem->periodicLength[2] = limDom[2][1] - limDom[2][0] + PSystem->partDist;
					limDom[2][0] -= PSystem->partDist*0.5;
					limDom[2][1] += PSystem->partDist*0.5;
					// Physical domain
					PSystem->physDomMinZ = limDom[2][0];
					PSystem->physDomMaxZ = limDom[2][1];
					PSystem->numBucketsZ = (long)((limDom[2][1]-limDom[2][0])/PSystem->bucketSide);
					// Adjust bucketSide to perfectly divide the domain without leftovers
					PSystem->bucketSide = PSystem->periodicLength[2]/(long)(PSystem->periodicLength[2]/PSystem->bucketSide);
					// Analysis domain is extended from the boundaries by one bucket width
					limDom[2][0] -= PSystem->bucketSide;
					limDom[2][1] += PSystem->bucketSide;
					PSystem->numBucketsZ += 2;
					// Adjust limits of X
					if(PSystem->wallType == boundaryWallType::PARTICLE) {
						limDom[0][0] -= PSystem->bucketSide;
						limDom[0][1] += PSystem->bucketSide;
					}
					else {
						limDom[0][0] = PSystem->domainMinX - PSystem->bucketSide;
						limDom[0][1] = PSystem->domainMaxX + PSystem->bucketSide;
					}
					// Physical domain
					PSystem->physDomMinX = limDom[0][0];
					PSystem->physDomMaxX = limDom[0][1];
					PSystem->numBucketsX = (long)((limDom[0][1]-limDom[0][0])/PSystem->bucketSide+1);
					limDom[0][1] = limDom[0][0] + PSystem->numBucketsX*PSystem->bucketSide;
					// Adjust limits of Y
					if(PSystem->wallType == boundaryWallType::PARTICLE) {
						limDom[1][0] -= PSystem->bucketSide;
						limDom[1][1] += PSystem->bucketSide;
					}
					else {
						limDom[1][0] = PSystem->domainMinY - PSystem->bucketSide;
						limDom[1][1] = PSystem->domainMaxY + PSystem->bucketSide;
					}
					// Physical domain
					PSystem->physDomMinY = limDom[1][0];
					PSystem->physDomMaxY = limDom[1][1];
					PSystem->numBucketsY = (long)((limDom[1][1]-limDom[1][0])/PSystem->bucketSide+1);
					limDom[1][1] = limDom[1][0] + PSystem->numBucketsY*PSystem->bucketSide;
				}
				
				std::cout << std::endl << "Periodic Domain Limits";
				break;
			}
		}
	}

	PSystem->domainMinX = limDom[0][0];	PSystem->domainMaxX = limDom[0][1];
	PSystem->domainMinY = limDom[1][0];	PSystem->domainMaxY = limDom[1][1];
	PSystem->domainMinZ = limDom[2][0];	PSystem->domainMaxZ = limDom[2][1];

}

// Update particle ID's in buckets
// Murotani et al., 2015. Performance improvements of differential operators code for MPS method on GPU.
void MpsBucket::updateParticlesID(MpsParticleSystem *PSystem, MpsParticle *Particles) {
	if((int)PSystem->dim == 2) {
#pragma omp parallel for
		for(int i=0; i<PSystem->numBucketsXY; i++) {	
			Particles->firstParticleInBucket[i] = -1;
			Particles->lastParticleInBucket[i] = -1;
		}
#pragma omp parallel for
		for(int i=0; i<Particles->numParticlesZero; i++) {	
			Particles->nextParticleInSameBucket[i] = -1;
		}
		for(int i=0; i<Particles->numParticles; i++) {
			if(Particles->particleType[i] == PSystem->ghost) continue;
			int ix = (int)((Particles->pos[i*3  ] - PSystem->domainMinX)*PSystem->invBucketSide + PSystem->epsilonZero);
			int iy = (int)((Particles->pos[i*3+1] - PSystem->domainMinY)*PSystem->invBucketSide + PSystem->epsilonZero);
			int ib = iy*PSystem->numBucketsX + ix;
			int j = Particles->lastParticleInBucket[ib];
			Particles->lastParticleInBucket[ib] = i;
			if(j == -1) {	Particles->firstParticleInBucket[ib] = i;	}
			else 		{	Particles->nextParticleInSameBucket[j] = i;}
		}
	}
	else {
#pragma omp parallel for
		for(int i=0; i<PSystem->numBucketsXYZ; i++) {
			Particles->firstParticleInBucket[i] = -1;
			Particles->lastParticleInBucket[i] = -1;
		}
#pragma omp parallel for
		for(int i=0; i<Particles->numParticlesZero; i++) {
			Particles->nextParticleInSameBucket[i] = -1;
		}
		for(int i=0; i<Particles->numParticles; i++) {
			if(Particles->particleType[i] == PSystem->ghost) continue;
			int ix = (int)((Particles->pos[i*3  ] - PSystem->domainMinX)*PSystem->invBucketSide);
			int iy = (int)((Particles->pos[i*3+1] - PSystem->domainMinY)*PSystem->invBucketSide);
			int iz = (int)((Particles->pos[i*3+2] - PSystem->domainMinZ)*PSystem->invBucketSide);
			int ib = iz*PSystem->numBucketsXY + iy*PSystem->numBucketsX + ix;
			int j = Particles->lastParticleInBucket[ib];
			Particles->lastParticleInBucket[ib] = i;
			if(j == -1) {	Particles->firstParticleInBucket[ib] = i;	}
			else 		{	Particles->nextParticleInSameBucket[j] = i;}
		}
	}
	// Copy data from periodic buckets to border buckets
	for(int b=0; b<PSystem->numBC; b++) {
		if(Particles->bucketTypeBC[b] == domainBC::PERIODIC) {
			copyDataBetweenBuckets(b, PSystem, Particles);
		}
	}
#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Copy data from periodic buckets to border buckets
void MpsBucket::copyDataBetweenBuckets(const int b, MpsParticleSystem *PSystem, MpsParticle *Particles) {
	bool periodicX =  Particles->periodicDirection[b*3] && !Particles->periodicDirection[b*3+1] && !Particles->periodicDirection[b*3+2];
	bool periodicY = !Particles->periodicDirection[b*3] &&  Particles->periodicDirection[b*3+1] && !Particles->periodicDirection[b*3+2];
	bool periodicZ = !Particles->periodicDirection[b*3] && !Particles->periodicDirection[b*3+1] &&  Particles->periodicDirection[b*3+2];
	// Periodic in X
	if(periodicX) {
		//for(int iz=0; iz<PSystem->numBucketsZ; iz++) {
		//	for(int iy=0; iy<PSystem->numBucketsY; iy++) {
		//		for(int ix=0; ix<PSystem->numBucketsX; ix+=PSystem->numBucketsX-1) {
		//			int ib = iz*PSystem->numBucketsXY + iy*PSystem->numBucketsX + ix;
#pragma omp parallel for
		for (int ib=0; ib<PSystem->numBucketsXYZ; ib++) {
			if(Particles->bucketPeriodicBC[ib] == -1) {
				Particles->firstParticleInBucket[ib] = Particles->firstParticleInBucket[ib+PSystem->numBucketsX-2];
				Particles->lastParticleInBucket[ib] = Particles->lastParticleInBucket[ib+PSystem->numBucketsX-2];
			}
			if(Particles->bucketPeriodicBC[ib] == 1) {
				Particles->firstParticleInBucket[ib] = Particles->firstParticleInBucket[ib-PSystem->numBucketsX+2];
				Particles->lastParticleInBucket[ib] = Particles->lastParticleInBucket[ib-PSystem->numBucketsX+2];
			}
		}
		//}}}
	}
	// Periodic in Y
	if(periodicY) {
		//for(int ix=0; ix<PSystem->numBucketsX; ix++) {
		//	for(int iz=0; iz<PSystem->numBucketsZ; iz++) {
		//		for(int iy=0; iy<PSystem->numBucketsY; iy+=PSystem->numBucketsY-1) {
		//			int ib = iz*PSystem->numBucketsXY + iy*PSystem->numBucketsX + ix;
#pragma omp parallel for
		for (int ib=0; ib<PSystem->numBucketsXYZ; ib++) {
			if(Particles->bucketPeriodicBC[ib] == -1) {
				Particles->firstParticleInBucket[ib] = Particles->firstParticleInBucket[ib+PSystem->numBucketsX*(PSystem->numBucketsY-2)];
				Particles->lastParticleInBucket[ib] = Particles->lastParticleInBucket[ib+PSystem->numBucketsX*(PSystem->numBucketsY-2)];
			}
			if(Particles->bucketPeriodicBC[ib] == 1) {
				Particles->firstParticleInBucket[ib] = Particles->firstParticleInBucket[ib-PSystem->numBucketsX*(PSystem->numBucketsY-2)];
				Particles->lastParticleInBucket[ib] = Particles->lastParticleInBucket[ib-PSystem->numBucketsX*(PSystem->numBucketsY-2)];
			}
		}
		//}}}
	}
	// Periodic in Z
	if(periodicZ) {
		//for(int iy=0; iy<PSystem->numBucketsY; iy++) {
		//	for(int ix=0; ix<PSystem->numBucketsX; ix++) {
		//		for(int iz=0; iz<PSystem->numBucketsZ; iz+=PSystem->numBucketsZ-1) {
		//			int ib = iz*PSystem->numBucketsXY + iy*PSystem->numBucketsX + ix;
#pragma omp parallel for
		for (int ib=0; ib<PSystem->numBucketsXYZ; ib++) {
			if(Particles->bucketPeriodicBC[ib] == -1) {
				Particles->firstParticleInBucket[ib] = Particles->firstParticleInBucket[ib+PSystem->numBucketsXY*(PSystem->numBucketsZ-2)];
				Particles->lastParticleInBucket[ib] = Particles->lastParticleInBucket[ib+PSystem->numBucketsXY*(PSystem->numBucketsZ-2)];
			}
			if(Particles->bucketPeriodicBC[ib] == 1) {
				Particles->firstParticleInBucket[ib] = Particles->firstParticleInBucket[ib-PSystem->numBucketsXY*(PSystem->numBucketsZ-2)];
				Particles->lastParticleInBucket[ib] = Particles->lastParticleInBucket[ib-PSystem->numBucketsXY*(PSystem->numBucketsZ-2)];
			}
		}
		//}}}
	}
}

// Return the bucket coordinates for particle "i"
void MpsBucket::bucketCoordinates(int &bx, int &by, int &bz,
	const double rxi, const double ryi, const double rzi, MpsParticleSystem *PSystem) {
	bx = (int)((rxi - PSystem->domainMinX)*PSystem->invBucketSide + PSystem->epsilonZero);
	by = (int)((ryi - PSystem->domainMinY)*PSystem->invBucketSide + PSystem->epsilonZero);
	bz = (int)((rzi - PSystem->domainMinZ)*PSystem->invBucketSide + PSystem->epsilonZero);
}