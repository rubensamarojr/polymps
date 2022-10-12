// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <experimental/filesystem> 	///< numeric_limits
#include <iostream>					///< cout
#include "MpsBoundaryCondition.h"

using namespace std;

// Constructor declaration
MpsBoundaryCondition::MpsBoundaryCondition()
{
}
// Destructor declaration
MpsBoundaryCondition::~MpsBoundaryCondition()
{
}

// Set Periodic Boundary Condition of the bucket
void MpsBoundaryCondition::setBucketBC(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	for(int iz=0; iz<PSystem->numBucketsZ; iz++) {
		for(int iy=0; iy<PSystem->numBucketsY; iy++) {
			for(int ix=0; ix<PSystem->numBucketsX; ix++) {
				int ib = iz*PSystem->numBucketsXY + iy*PSystem->numBucketsX + ix;
				Particles->bucketPeriodicBC[ib] = 0;
				for(int b=0; b<PSystem->numBC; b++) {
					if(Particles->bucketTypeBC[b] == domainBC::PERIODIC) {
						// Periodic in X
						if(Particles->periodicDirection[b*3]) {
							if(ix == 0) { // Left border of domain
								Particles->bucketPeriodicBC[ib] = -1;
							}
							if(ix == PSystem->numBucketsX-1) { // Right border of domain
								Particles->bucketPeriodicBC[ib] = 1;
							}
						}
						// Periodic in Y
						if(Particles->periodicDirection[b*3+1]) {
							if(iy == 0) { // Bottom border of domain
								Particles->bucketPeriodicBC[ib] = -1;
							}
							if(iy == PSystem->numBucketsY-1) { // Top border of domain
								Particles->bucketPeriodicBC[ib] = 1;
							}
						}
						// Periodic in Z
						if(Particles->periodicDirection[b*3+2]) {
							if(iz == 0) { // Back border of domain
								Particles->bucketPeriodicBC[ib] = -1;
							}
							if(iz == PSystem->numBucketsZ-1) { // Front border of domain
								Particles->bucketPeriodicBC[ib] = 1;
							}
						}
					}
					//std::cout << std::endl << "ib: " << ib << " ix: " << ix
					//<< " iy: " << iy << " iz: " << iz;
					//std::cout << " bc: " << Particles->bucketPeriodicBC[ib];
				}
			}
		}
	}
}

// Update type of particle
void MpsBoundaryCondition::updateParticleBC(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
/*
	double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
	double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
	double ni = 0.0; double wSum = 0.0;
	Particles->numNeigh[i] = 0;
	// Add Number of neighbors due Wall polygon
	Particles->numNeigh[i] += numNeighWallContribution[i];
	
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
			if(dstij2 < PSystem->reL2 && (dstij2 < dstimj2 || PSystem->wallType == boundaryWallType::PARTICLE)) {
			if(j != i) {
				Particles->numNeigh[i] += 1;
				//double dst = sqrt(dstij2);
				//double wL = weight(dst, PSystem->reL*PSystem->invPartDist, PSystem->weightType);
				//npcdDeviation[i*3  ] += v0ij*wL*PSystem->invPartDist;
				//npcdDeviation[i*3+1] += v1ij*wL*PSystem->invPartDist;
				//npcdDeviation[i*3+2] += v2ij*wL*PSystem->invPartDist;
				//wSum += wL;
				if(dstij2 < PSystem->reS2) {
					double dst = sqrt(dstij2);
					double wS = weight(dst, PSystem->reS, PSystem->weightType);
					ni += wS;
					dst = dst*PSystem->invPartDist;
					wS = weight(dst, PSystem->reS*PSystem->invPartDist, PSystem->weightType);;
					npcdDeviation[i*3  ] += v0ij*wS*PSystem->invPartDist;
					npcdDeviation[i*3+1] += v1ij*wS*PSystem->invPartDist;
					npcdDeviation[i*3+2] += v2ij*wS*PSystem->invPartDist;
					wSum += wS;
			}}}
			j = Particles->nextParticleInSameBucket[j];
			if(j == -1) break;
		}
	}}}

	double mi;
	if(PTYPE[i] == 1) mi = DNS_FL1;
	else mi = DNS_FL2;
	//if(Particles->particleType[i]==PSystem->fluid)
	//	mi = Dns[partType::FLUID];
	//else
	//	mi = Dns[partType::WALL];

	if(PSystem->pndType == calcPNDType::SUM_WIJ || PSystem->pndType == calcPNDType::MEAN_SUM_WIJ)
//		if(PSystem->pndType == calcPNDType::SUM_WIJ)
	{
		// PND due particles and Wall polygon
		Particles->pndi[i] = ni + pndWallContribution[i];
	}

//		if(Particles->particleType[i] == PSystem->wall) {
		// PND due particles and Wall polygon
//			Particles->pndi[i] = ni + pndWallContribution[i];
//			if(Particles->pndi[i] < PSystem->pndSmallZero)
//				Particles->pndi[i] = PSystem->pndSmallZero;

//				Particles->pndi[i] = PSystem->pndSmallZero*pow((press[i]*gamma/(mi*PSystem->coeffPressWCMPS)+1),gamma);
//		}

	// Add PND due Wall polygon
	Particles->pndSmall[i] = ni + pndWallContribution[i];
	// Prevent Particles->pndSmall[i] = 0
//		if(Particles->numNeigh[i] > 1) {
	if(wSum > PSystem->epsilonZero) {
		//npcdDeviation[i*3  ] /= Particles->pndSmall[i];
		//npcdDeviation[i*3+1] /= Particles->pndSmall[i];
		//npcdDeviation[i*3+2] /= Particles->pndSmall[i];
		npcdDeviation[i*3  ] /= wSum;
		npcdDeviation[i*3+1] /= wSum;
		npcdDeviation[i*3+2] /= wSum;
	}

	Particles->npcdDeviation2[i] = npcdDeviation[i*3]*npcdDeviation[i*3]+npcdDeviation[i*3+1]*npcdDeviation[i*3+1]+npcdDeviation[i*3+2]*npcdDeviation[i*3+2];

	//deviationDotPolygonNormal[i] = npcdDeviation[i*3]*polygonNormal[i*3]+npcdDeviation[i*3+1]*polygonNormal[i*3+1]+npcdDeviation[i*3+2]*polygonNormal[i*3+2];
	if(npcdDeviation[i*3]*polygonNormal[i*3]+npcdDeviation[i*3+1]*polygonNormal[i*3+1]+npcdDeviation[i*3+2]*polygonNormal[i*3+2] < 0.0)
		deviationDotPolygonNormal[i] = 1;
	else
		deviationDotPolygonNormal[i] = -1;

	// PSystem->coeffPressWCMPS = soundSpeed*soundSpeed
	double pressure = 0.0;
*/		
		if(Particles->particleType[i] == PSystem->dummyWall) {
			Particles->particleBC[i] = PSystem->other;
			continue;
		}
		// First check based on particle number density
		if(Particles->pndSmall[i] < PSystem->betaPnd) {
		//if(Particles->pndi[i] < PSystem->betaPnd) {
			Particles->particleBC[i] = PSystem->surface;
		}
		else {
			Particles->particleBC[i] = PSystem->inner;
		}

		if(PSystem->freeSurfType == calcBCType::PND_NEIGH) {
			if(Particles->pndSmall[i] < PSystem->betaPnd && Particles->numNeigh[i] < PSystem->betaNeigh) {
			//if(Particles->pndi[i] < PSystem->betaPnd && Particles->numNeigh[i] < PSystem->betaNeigh) {
				Particles->particleBC[i] = PSystem->surface;
			}
			else {
				Particles->particleBC[i] = PSystem->inner;
			}
		}
		else if(PSystem->freeSurfType == calcBCType::PND_NPCD) {
			// Boundary particle verification based on relative distance and weight (NPCD)
			// 2016 - Fluid interface detection technique based on neighborhood particles 
			// centroid deviation (NPCD) for particle methods
			if(Particles->particleBC[i] == PSystem->surface) {
				if(Particles->numNeigh[i] > 4 && Particles->npcdDeviation2[i] < PSystem->delta2) {
					Particles->particleBC[i] = PSystem->inner;
				}
			}
		}
		else if(PSystem->freeSurfType == calcBCType::PND_ARC) {
			double normalXi = Particles->normal[i*3  ];	double normalYi = Particles->normal[i*3+1];	double normalZi = Particles->normal[i*3+2];
			double norm2 = normalXi*normalXi + normalYi*normalYi + normalZi*normalZi;
			// 2017 - A multiphase MPS solver for modeling multi-fluid interaction with 
			// free surface and its application in oil spill
			//if((Particles->pndSmall[i] >= PSystem->betaPnd || Particles->numNeigh[i] >= PSystem->betaNeigh) && norm2 <= PSystem->normThreshold2) {
			/*if(Particles->pndSmall[i] >= PSystem->betaPnd && Particles->numNeigh[i] >= PSystem->betaNeigh && norm2 <= PSystem->normThreshold2) {
			//if(Particles->pndSmall[i] >= PSystem->betaPnd || Particles->numNeigh[i] >= PSystem->betaNeigh) {
			*/
			if(Particles->pndi[i] >= PSystem->betaPnd && Particles->numNeigh[i] >= PSystem->betaNeigh && norm2 <= PSystem->normThreshold2) {
				Particles->particleBC[i] = PSystem->inner;
			}
			else {
				double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
				double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
				//double normalXi = Particles->normal[i*3  ];	double normalYi = Particles->normal[i*3+1];	double normalZi = Particles->normal[i*3+2];
				if(norm2 > PSystem->epsilonZero) {
					double norm = sqrt(norm2);
					normalXi /= norm;	normalYi /= norm;	normalZi /= norm;
				}
				else {
					normalXi = 0.0;	normalYi = 0.0;	normalZi = 0.0;
				}
				// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
				double Rref_i[9], normaliw[3], normalMod2;
				// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
				normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
				normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];
				if(normalMod2 > PSystem->epsilonZero) {
					double normalMod = sqrt(normalMod2);
					normaliw[0] /= normalMod; normaliw[1] /= normalMod; normaliw[2] /= normalMod;
				}
				else {
					normaliw[0] = 0.0; normaliw[1] = 0.0; normaliw[2] = 0;
				}

				// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
				Rref_i[0] = 1.0 - 2.0*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2.0*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2.0*normaliw[0]*normaliw[2];
				Rref_i[3] = 0.0 - 2.0*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2.0*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2.0*normaliw[1]*normaliw[2];
				Rref_i[6] = 0.0 - 2.0*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2.0*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2.0*normaliw[2]*normaliw[2];

				// Normal mirror i
				double normalMirrorXi = Rref_i[0]*normalXi + Rref_i[1]*normalYi + Rref_i[2]*normalZi;
				double normalMirrorYi = Rref_i[3]*normalXi + Rref_i[4]*normalYi + Rref_i[5]*normalZi;
				double normalMirrorZi = Rref_i[6]*normalXi + Rref_i[7]*normalYi + Rref_i[8]*normalZi;
				double normMirror2 = normalMirrorXi*normalMirrorXi + normalMirrorYi*normalMirrorYi + normalMirrorZi*normalMirrorZi;
				if(normMirror2 > PSystem->epsilonZero) {
					double normMirror = sqrt(normMirror2);
					normalMirrorXi /= normMirror;	normalMirrorYi /= normMirror;	normalMirrorZi /= normMirror;
				}
				else {
					normalMirrorXi = 0.0;	normalMirrorYi = 0.0;	normalMirrorZi = 0.0;
				}

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

						// Real particle i and neighbor j
						// If j is inside the neighborhood of i and 
						// is not at the same side of im (avoid real j in the virtual neihborhood)
						if(dstij2 < PSystem->reS2 && (dstij2 < dstimj2 || PSystem->wallType == boundaryWallType::PARTICLE)) {
						if(j != i) {
							double v0inj, v1inj, v2inj, dstinj2;
							Particles->sqrDistBetweenParticles(j, posXi + PSystem->partDist*normalXi, posYi + PSystem->partDist*normalYi, 
								posZi + PSystem->partDist*normalZi, v0inj, v1inj, v2inj, dstinj2, plx, ply, plz);
							
							double rijn = v0ij*normalXi + v1ij*normalYi + v2ij*normalZi;

							double ang = acos(rijn/sqrt(dstij2));

							/*
							if (dstij2 >= PSystem->dstThreshold2 && dstinj2 < PSystem->hThreshold2) {
								Particles->particleBC[i] = PSystem->inner;
								//goto endloop;
							}
							else if (dstij2 < PSystem->dstThreshold2 && ang < 3.14159265/4.0) {
							//else if (dstij2 < PSystem->dstThreshold2 && rijn*rijn < dstij2*0.5) {
								Particles->particleBC[i] = PSystem->inner;
								//goto endloop;
							}
							else {
								Particles->particleBC[i] = PSystem->surface;
								goto endloop;
							}
							*/
							if (dstij2 >= PSystem->dstThreshold2 && dstinj2 < PSystem->hThreshold2) {
								Particles->particleBC[i] = PSystem->inner;
								goto endloop;
							}
							else if (dstij2 < PSystem->dstThreshold2 && ang < PSystem->thetaArc) {
								Particles->particleBC[i] = PSystem->inner;
								goto endloop;
							}
							else {
								Particles->particleBC[i] = PSystem->surface;
							}
/*
							//if ((dstij2 < PSystem->dstThreshold2 && ang < PSystem->thetaArc) || (Particles->numNeigh[i] >= PSystem->betaNeigh)) {
							if (ang < PSystem->thetaArc && Particles->numNeigh[i] >= 4) {
							//else if (dstij2 < PSystem->dstThreshold2 && rijn*rijn < dstij2*0.5) {
								Particles->particleBC[i] = PSystem->inner;
								goto endloop;
							}
							else {
								Particles->particleBC[i] = PSystem->surface;
								//goto endloop;
							}*/
						}}

						if(PSystem->wallType == boundaryWallType::POLYGON){
							// Virtual particle i and real neighbor j
							// If j is inside the neighborhood of i and im (intersection) and 
							// is not at the same side of im (avoid real j in the virtual neihborhood)
							if(dstij2 < PSystem->reS2 && dstimj2 < PSystem->reS2 && dstij2 < dstimj2) {
								if(j != i) {
									double v0imnj, v1imnj, v2imnj, dstimnj2;
									Particles->sqrDistBetweenParticles(j, posMirrorXi + PSystem->partDist*normalMirrorXi, posMirrorYi + PSystem->partDist*normalMirrorYi, 
										posMirrorZi + PSystem->partDist*normalMirrorZi, v0imnj, v1imnj, v2imnj, dstimnj2, plx, ply, plz);
									
									double rimjn = v0imj*normalMirrorXi + v1imj*normalMirrorYi + v2imj*normalMirrorZi;

									double angm = acos(rimjn/sqrt(dstimj2));

									/*
									if (dstimj2 >= PSystem->dstThreshold2 && dstimnj2 < PSystem->hThreshold2) {
										Particles->particleBC[i] = PSystem->inner;
										//goto endloop;
									}
									//else if (dstimj2 < PSystem->dstThreshold2 && rimjn*rimjn < dstimj2*0.6675*0.6675) {
									else if (dstimj2 < PSystem->dstThreshold2 && angm < PSystem->thetaArc) {
									//else if (dstimj2 < PSystem->dstThreshold2 && rimjn*rimjn < dstimj2*0.5) {
										Particles->particleBC[i] = PSystem->inner;
										//goto endloop;
									}
									else {
										Particles->particleBC[i] = PSystem->surface;
										goto endloop;
									}

									*/

									if (dstimj2 >= PSystem->dstThreshold2 && dstimnj2 < PSystem->hThreshold2) {
										Particles->particleBC[i] = PSystem->inner;
										goto endloop;
									}
									else if (dstimj2 < PSystem->dstThreshold2 && angm < PSystem->thetaArc) {
										Particles->particleBC[i] = PSystem->inner;
										goto endloop;
									}
									else {
										Particles->particleBC[i] = PSystem->surface;
									}
									/*
									//if ((dstimj2 < PSystem->dstThreshold2 && angm < 3.14159265/4.0) || (Particles->numNeigh[i] >= PSystem->betaNeigh)) {
									if (angm < PSystem->thetaArc && Particles->numNeigh[i] >= 4) {
										Particles->particleBC[i] = PSystem->inner;
										goto endloop;
									}
									else {
										Particles->particleBC[i] = PSystem->surface;
										//goto endloop;
									}*/
								}
							}
						}
						j = Particles->nextParticleInSameBucket[j];
						if(j == -1) break;
					}
				}}}

				endloop: ;
			}
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}