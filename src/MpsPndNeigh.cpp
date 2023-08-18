// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <experimental/filesystem> 	///< numeric_limits
#include <iostream>					///< cout
#include "MpsPndNeigh.h"

using namespace std;

// Constructor declaration
MpsPndNeigh::MpsPndNeigh()
{
}
// Destructor declaration
MpsPndNeigh::~MpsPndNeigh()
{
}


// Set initial PND, number of neighbors and weighted deviation
void MpsPndNeigh::setInitialPndNumberOfNeighNPCD(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double wSum = 0.0;
		Particles->numNeigh[i] = 0;
		Particles->npcdDeviation[i*3] = Particles->npcdDeviation[i*3+1] = Particles->npcdDeviation[i*3+2] = 0.0;
		
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
						if(dstij2 < PSystem->reS2) {
							double dst = sqrt(dstij2);
							double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
							Particles->pndi[i] += wS;
							//dst = dst*PSystem->invPartDist;
							//wS = Particles->weight(dst, PSystem->reS*PSystem->invPartDist, PSystem->weightType);
							//Particles->npcdDeviation[i*3  ] += v0ij*wS*PSystem->invPartDist;
							//Particles->npcdDeviation[i*3+1] += v1ij*wS*PSystem->invPartDist;
							//Particles->npcdDeviation[i*3+2] += v2ij*wS*PSystem->invPartDist;
							Particles->npcdDeviation[i*3  ] += v0ij*wS;
							Particles->npcdDeviation[i*3+1] += v1ij*wS;
							Particles->npcdDeviation[i*3+2] += v2ij*wS;
							wSum += wS;
						}
					}
				}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// Add PND due wall polygon
		Particles->pndi[i] += Particles->pndWallContribution[i];
		if(Particles->particleType[i] == PSystem->wall)
			Particles->pndi[i] = PSystem->pndSmallZero;
		Particles->pndSmall[i] = Particles->pndi[i];
		Particles->pndki[i] = Particles->pndi[i];
		// Add Number of neighbors due wall polygon
		Particles->numNeigh[i] += Particles->numNeighWallContribution[i];

		if(wSum > PSystem->epsilonZero) {
			Particles->npcdDeviation[i*3  ] /= Particles->pndSmall[i];
			Particles->npcdDeviation[i*3+1] /= Particles->pndSmall[i];
			Particles->npcdDeviation[i*3+2] /= Particles->pndSmall[i];
			//Particles->npcdDeviation[i*3  ] /= wSum;
			//Particles->npcdDeviation[i*3+1] /= wSum;
			//Particles->npcdDeviation[i*3+2] /= wSum;
		}

		Particles->npcdDeviation2[i] = Particles->npcdDeviation[i*3]*Particles->npcdDeviation[i*3] + Particles->npcdDeviation[i*3+1]*Particles->npcdDeviation[i*3+1] +
			Particles->npcdDeviation[i*3+2]*Particles->npcdDeviation[i*3+2];

		//Particles->deviationDotPolygonNormal[i] = Particles->npcdDeviation[i*3]*Particles->polygonNormal[i*3]+Particles->npcdDeviation[i*3+1]*Particles->polygonNormal[i*3+1]+Particles->npcdDeviation[i*3+2]*Particles->polygonNormal[i*3+2];
		if(Particles->npcdDeviation[i*3]*Particles->polygonNormal[i*3]+Particles->npcdDeviation[i*3+1]*Particles->polygonNormal[i*3+1]+Particles->npcdDeviation[i*3+2]*Particles->polygonNormal[i*3+2] < 0.0)
			Particles->deviationDotPolygonNormal[i] = 1;
		else
			Particles->deviationDotPolygonNormal[i] = -1;
	}
}


// Calculates weighted deviation due to the closest polygon wall
void MpsPndNeigh::calcWallNPCD(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
			double npcdDeviationXi = 0.0;	double npcdDeviationYi = 0.0;	double npcdDeviationZi = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			// Wall gradient Mitsume`s model
			double Rref_i[9], normaliw[3], normaliwSqrt;
			// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
			normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
			normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);

			if(normaliwSqrt > PSystem->epsilonZero) {
				normaliw[0] = normaliw[0]/normaliwSqrt;
				normaliw[1] = normaliw[1]/normaliwSqrt;
				normaliw[2] = normaliw[2]/normaliwSqrt;
			}
			else {
				normaliw[0] = 0;
				normaliw[1] = 0;
				normaliw[2] = 0;
			}

		    // Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
			Rref_i[0] = 1.0 - 2.0*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2.0*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2.0*normaliw[0]*normaliw[2];
			Rref_i[3] = 0.0 - 2.0*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2.0*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2.0*normaliw[1]*normaliw[2];
			Rref_i[6] = 0.0 - 2.0*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2.0*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2.0*normaliw[2]*normaliw[2];

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

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < PSystem->reS2 && dstimj2 < PSystem->reS2 && dstij2 < dstimj2) {
						if(j != i) {
							double dst = sqrt(dstimj2);
							double wS = Particles->weightGradient(dst, PSystem->reS, PSystem->weightType);
							npcdDeviationXi += v0imj*wS;
							npcdDeviationYi += v1imj*wS;
							npcdDeviationZi += v2imj*wS;
						}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}

			}}}
			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			Particles->sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
			
			if(dstimi2 < PSystem->reS2) {
				double dst = sqrt(dstimi2);
				double wS = Particles->weightGradient(dst, PSystem->reS, PSystem->weightType);
				npcdDeviationXi += v0imi*wS;
				npcdDeviationYi += v1imi*wS;
				npcdDeviationZi += v2imi*wS;
			}
			Particles->npcdDeviation[i*3  ] += Rref_i[0]*npcdDeviationXi + Rref_i[1]*npcdDeviationYi + Rref_i[2]*npcdDeviationZi;
			Particles->npcdDeviation[i*3+1] += Rref_i[3]*npcdDeviationXi + Rref_i[4]*npcdDeviationYi + Rref_i[5]*npcdDeviationZi;
			Particles->npcdDeviation[i*3+2] += Rref_i[6]*npcdDeviationXi + Rref_i[7]*npcdDeviationYi + Rref_i[8]*npcdDeviationZi;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}


// Calculates PND, number of neighbors and weighted deviation
void MpsPndNeigh::calcPndnNeighNPCD(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double ni = 0.0; double wSum = 0.0;
		Particles->numNeigh[i] = 0;
		// Add Number of neighbors due Wall polygon
		Particles->numNeigh[i] += Particles->numNeighWallContribution[i];
		
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
						//double dst = sqrt(dst2);
						//double wL = Particles->weight(dst, PSystem->reL*PSystem->invPartDist, PSystem->weightType);
						//Particles->npcdDeviation[i*3  ] += v0*wL*PSystem->invPartDist;
						//Particles->npcdDeviation[i*3+1] += v1*wL*PSystem->invPartDist;
						//Particles->npcdDeviation[i*3+2] += v2*wL*PSystem->invPartDist;
						//wSum += wL;
						if(dstij2 < PSystem->reS2) {
							double dst = sqrt(dstij2);
							double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
							ni += wS;
							//dst = dst*PSystem->invPartDist;
							//wS = Particles->weight(dst, PSystem->reS*PSystem->invPartDist, PSystem->weightType);
							//Particles->npcdDeviation[i*3  ] += v0ij*wS*PSystem->invPartDist;
							//Particles->npcdDeviation[i*3+1] += v1ij*wS*PSystem->invPartDist;
							//Particles->npcdDeviation[i*3+2] += v2ij*wS*PSystem->invPartDist;
							Particles->npcdDeviation[i*3  ] += v0ij*wS;
							Particles->npcdDeviation[i*3+1] += v1ij*wS;
							Particles->npcdDeviation[i*3+2] += v2ij*wS;
							wSum += wS;
						}
					}
				}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

//		double mi;
//		if(PTYPE[i] == 1) mi = PSystem->DNS_FL1;
//		else mi = PSystem->DNS_FL2;
		//if(Particles->particleType[i]==fluid)
		//	mi = Dns[partType::FLUID];
		//else
		//	mi = Dns[partType::WALL];

		if(PSystem->pndType == calcPNDType::SUM_WIJ || PSystem->pndType == calcPNDType::MEAN_SUM_WIJ) {
//		if(PSystem->pndType == calcPNDType::SUM_WIJ) {
			// PND at initial of step k
			Particles->pndki[i] = Particles->pndi[i];
			// New PND due particles and Wall polygon
			Particles->pndi[i] = ni + Particles->pndWallContribution[i];
		}

//		if(Particles->particleType[i] == PSystem->wall) {
			// PND due particles and Wall polygon
//			Particles->pndi[i] = ni + Particles->pndWallContribution[i];
//			if(Particles->pndi[i] < PSystem->pndSmallZero)
//				Particles->pndi[i] = PSystem->pndSmallZero;

//				Particles->pndi[i] = PSystem->pndSmallZero*pow((Particles->press[i]*PSystem->gamma/(mi*PSystem->coeffPressWCMPS)+1),PSystem->gamma);
//		}

		// Add PND due Wall polygon
		Particles->pndSmall[i] = ni + Particles->pndWallContribution[i];
		// Prevent Particles->pndSmall[i] = 0
//		if(Particles->numNeigh[i]>1) {
		if(wSum > PSystem->epsilonZero) {
			Particles->npcdDeviation[i*3  ] /= Particles->pndSmall[i];
			Particles->npcdDeviation[i*3+1] /= Particles->pndSmall[i];
			Particles->npcdDeviation[i*3+2] /= Particles->pndSmall[i];
			//Particles->npcdDeviation[i*3  ] /= wSum;
			//Particles->npcdDeviation[i*3+1] /= wSum;
			//Particles->npcdDeviation[i*3+2] /= wSum;
		}

		Particles->npcdDeviation2[i] = Particles->npcdDeviation[i*3]*Particles->npcdDeviation[i*3]+Particles->npcdDeviation[i*3+1]*Particles->npcdDeviation[i*3+1]+Particles->npcdDeviation[i*3+2]*Particles->npcdDeviation[i*3+2];

		//Particles->deviationDotPolygonNormal[i] = Particles->npcdDeviation[i*3]*Particles->polygonNormal[i*3]+Particles->npcdDeviation[i*3+1]*Particles->polygonNormal[i*3+1]+Particles->npcdDeviation[i*3+2]*Particles->polygonNormal[i*3+2];
		if(Particles->npcdDeviation[i*3]*Particles->polygonNormal[i*3]+Particles->npcdDeviation[i*3+1]*Particles->polygonNormal[i*3+1]+Particles->npcdDeviation[i*3+2]*Particles->polygonNormal[i*3+2]< 0.0)
			Particles->deviationDotPolygonNormal[i] = 1;
		else
			Particles->deviationDotPolygonNormal[i] = -1;
		
		// First check based on particle number density
//		if(Particles->pndSmall[i] < PSystem->pndThreshold*PSystem->pndSmallZero)
//			Particles->particleBC[i] = PSystem->surface;
//		else
//			Particles->particleBC[i] = PSystem->inner;

		// Boundary particle verification based on relative distance and weight (NPCD)
//		if(Particles->particleBC[i] == PSystem->surface) {
//			if(Particles->numNeigh[i] > 4 && Particles->npcdDeviation2[i] < PSystem->delta2)
//			{
//				Particles->particleBC[i] = PSystem->inner;
				//printf(" inner %d \n", i);
//			}
//		}

//		if(Particles->pndSmall[i] < PSystem->pndThreshold*PSystem->pndSmallZero && Particles->numNeigh[i] < PSystem->neighThreshold*PSystem->numNeighZero)
//			Particles->particleBC[i] = PSystem->surface;
//		else
//			Particles->particleBC[i] = PSystem->inner;
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Calculates PND based on continuity equation including a diffusive term. Contribution of the neighboring particles.
void MpsPndNeigh::calcPndDiffusiveTerm(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	// PSystem->coeffViscMultiphase = 2.0*PSystem->dim/(PSystem->pndLargeZero*PSystem->lambdaZero);
	// double C1 = PSystem->diffusiveCoef*PSystem->timeStep*PSystem->soundSpeed*PSystem->soundSpeed*PSystem->coeffViscMultiphase/(PSystem->pndLargeZero);
	// double C2 = PSystem->diffusiveCoef*PSystem->partDist*PSystem->soundSpeed*PSystem->coeffViscMultiphase/(PSystem->pndLargeZero);
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
//	if(Particles->particleType[i] == PSystem->fluid) {
		double Di = 0.0; double DivV = 0.0; double flagDi = 1.0; 
		// double pndAux = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double MC[9];
		MC[0] = Particles->correcMatrixRow1[i*3];	MC[1] = Particles->correcMatrixRow1[i*3+1];	MC[2] = Particles->correcMatrixRow1[i*3+2];
		MC[3] = Particles->correcMatrixRow2[i*3];	MC[4] = Particles->correcMatrixRow2[i*3+1];	MC[5] = Particles->correcMatrixRow2[i*3+2];
		MC[6] = Particles->correcMatrixRow3[i*3];	MC[7] = Particles->correcMatrixRow3[i*3+1];	MC[8] = Particles->correcMatrixRow3[i*3+2];

		double ni = Particles->pndi[i];
		if(ni < PSystem->epsilonZero) continue;
		double Pi = Particles->press[i];
		
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
	//				if(j != i && Particles->particleType[j] == PSystem->fluid) {
						double dst = sqrt(dstij2);
						double wL = Particles->weight(dst, PSystem->reL, PSystem->weightType);
						// double nj = Particles->pndi[j];
						if(Particles->particleType[i] == PSystem->fluid && Particles->particleType[j] == PSystem->fluid) {
						//if(Particles->particleType[i] == PSystem->inner) {
	//					if(Particles->particleType[i] == PSystem->fluid && Particles->particleType[j] == PSystem->fluid && Particles->particleBC[i] == PSystem->inner) {
							// PSystem->coeffViscMultiphase = 2.0*PSystem->dim/(PSystem->pndLargeZero*PSystem->lambdaZero);
	//						double pgh = Particles->RHO[i]*(PSystem->gravityX*v0+PSystem->gravityY*v1+PSystem->gravityZ*v2);
	//						double CB = Particles->RHO[i]*PSystem->soundSpeed*PSystem->soundSpeed/PSystem->gamma;
	//						double nijH = PSystem->pndSmallZero*(pow((pgh+1.0)/CB,1.0/PSystem->gamma)-1.0);
	//						double nijH = PSystem->pndSmallZero*(pow((pgh)/CB+1.0,1.0/PSystem->gamma)-1.0);
	//
						//	pow( ( pgh + 1.0 ) / CB - 1.0 , 1.0  )

							//double CB = PSystem->soundSpeed*PSystem->soundSpeed*Particles->RHO[i];
							//double nijH = PSystem->pndSmallZero*((PijH+1.0)/CB-1);
	//						if(isnan(nijH) == 0)
	//							Di += C1*(nj - ni - nijH)*wL;
	//						else
								////////Di += C1*(nj - ni)*wL;
							//if(isnan(nijH) == 1)
							//if(i == 200)
							//	printf(" pgh %e CB %e ni %e nj %e nijH %e PSystem->reS %e \n", pgh, CB, ni, nj, nijH, PSystem->pndSmallZero*(pow((pgh+1)/CB,1/PSystem->gamma)-1));
							// PND
	//						Di += C1*(nj-ni)*wL;
							//Di += C2*(nj-ni)*wL;
							// Pressure
							// Delta Voronoi smoothed particle hydrodynamics, δ-VSPH
							// https://doi.org/10.1016/j.jcp.2019.109000
							double pgh = -2.0*(Particles->RHO[i]*Particles->RHO[j]/(Particles->RHO[i]+Particles->RHO[j]))*(PSystem->gravityX*v0ij+PSystem->gravityY*v1ij+PSystem->gravityZ*v2ij);
							double Pj = Particles->press[j];
							Di += PSystem->timeStep/Particles->RHO[i]*PSystem->coeffViscMultiphase*(Pj-Pi-pgh)*wL;
							//Di += (PSystem->partDist/PSystem->soundSpeed)/Particles->RHO[i]*PSystem->coeffViscMultiphase*(Pj-Pi+pgh)*wL;
						}
						//else
						//	flagDi = 0.0;
						if(dstij2 < PSystem->reS2) {
							double vijx = Particles->vel[j*3  ]-velXi;
							double vijy = Particles->vel[j*3+1]-velYi;
							double vijz = Particles->vel[j*3+2]-velZi;
							double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
							//if(ni > PSystem->epsilonZero)
							//{
							//	DivV += (PSystem->dim/PSystem->pndSmallZero)*(nj/ni)*(vijx*v0ij+vijy*v1ij+vijz*v2ij)*wS/dstij2;
							//}
							if(PSystem->gradientCorrection == false) {
								if(ni > PSystem->epsilonZero) {
									DivV += (PSystem->dim/PSystem->pndSmallZero)*(Particles->pndi[j]/ni)*(vijx*v0ij+vijy*v1ij+vijz*v2ij)*wS/dstij2;
								}
							}
							else {
								double v0ijC = (v0ij*MC[0] + v1ij*MC[1] + v2ij*MC[2]);
								double v1ijC = (v0ij*MC[3] + v1ij*MC[4] + v2ij*MC[5]);
								double v2ijC = (v0ij*MC[6] + v1ij*MC[7] + v2ij*MC[8]);
								DivV += (PSystem->dim/PSystem->pndSmallZero)*(vijx*v0ijC+vijy*v1ijC+vijz*v2ijC)*wS/dstij2;
							}

	//						M1[0][0] += (PSystem->dim/PSystem->pndSmallZero)*(nj/ni)*(v0*vijx)*wS/dst2; M1[0][1] += (PSystem->dim/PSystem->pndSmallZero)*(nj/ni)*(v0*vijy)*wS/dst2; M1[0][2] += (PSystem->dim/PSystem->pndSmallZero)*(nj/ni)*(v0*vijz)*wS/dst2;
	//						M1[1][0] += (PSystem->dim/PSystem->pndSmallZero)*(nj/ni)*(v1*vijx)*wS/dst2; M1[1][1] += (PSystem->dim/PSystem->pndSmallZero)*(nj/ni)*(v1*vijy)*wS/dst2; M1[1][2] += (PSystem->dim/PSystem->pndSmallZero)*(nj/ni)*(v1*vijz)*wS/dst2;
	//						M1[2][0] += (PSystem->dim/PSystem->pndSmallZero)*(nj/ni)*(v2*vijx)*wS/dst2; M1[2][1] += (PSystem->dim/PSystem->pndSmallZero)*(nj/ni)*(v2*vijy)*wS/dst2; M1[2][2] += (PSystem->dim/PSystem->pndSmallZero)*(nj/ni)*(v2*vijz)*wS/dst2;

	//						if(Particles->particleType[i] == PSystem->wall)
	//						if(Particles->particleType[i] == PSystem->wall || Particles->particleBC[i] == PSystem->surface)
	//							pndAux += wS;
						}
					}
				}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}


//		DivV = 	Particles->correcMatrixRow1[i*3  ]*M1[0][0] + Particles->correcMatrixRow1[i*3+1]*M1[1][0] + Particles->correcMatrixRow1[i*3+2]*M1[2][0] +
//				Particles->correcMatrixRow2[i*3  ]*M1[0][1] + Particles->correcMatrixRow2[i*3+1]*M1[1][1] + Particles->correcMatrixRow2[i*3+2]*M1[2][1] +
//				Particles->correcMatrixRow3[i*3  ]*M1[0][2] + Particles->correcMatrixRow3[i*3+1]*M1[1][2] + Particles->correcMatrixRow3[i*3+2]*M1[2][2];

		Particles->pndAuxVar1[i] = Particles->pndi[i]*(1.0+PSystem->timeStep*(-DivV+Di*flagDi));


//		if(isnan(DivV) || isnan(Di))
//			printf(" i %d \n", i);
		Particles->velDivergence[i] = DivV;
		Particles->diffusiveTerm[i] = Di;

	//Particles->pndAuxVar1[i] = Particles->pndi[i]*(1.0+PSystem->timeStep*(-(1.0-PSystem->diffusiveCoef)*DivV+Di*flagDi));
	//Particles->pndAuxVar1[i] = Particles->pndSmall[i]*(1.0+PSystem->timeStep*(-DivV+Di*flagDi));
	//Particles->pndAuxVar1[i] =     PSystem->pndSmallZero*(1.0+PSystem->timeStep*(-DivV+Di*flagDi)); // Ruim
//		if(Particles->particleType[i] == PSystem->wall)
//		{
//			if(pndAux < PSystem->pndSmallZero)
//				pndAux = PSystem->pndSmallZero;
//			Particles->pndAuxVar1[i] = pndAux;
//		}
//		if(Particles->particleBC[i] == PSystem->surface)
//		{
//			Particles->pndAuxVar1[i] = pndAux;
//		}
	}
//#pragma omp parallel for
//	for(int i=0; i<Particles->numParticles; i++) {
/////	if(Particles->particleType[i] == PSystem->fluid) {
//		Particles->pndi[i] = Particles->pndAuxVar1[i];
//		Particles->pndAuxVar1[i]=0.0;
//	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Calculates PND based on continuity equation including a diffusive term. Contribution of the closest polygon wall. Free-slip boundary condition.
void MpsPndNeigh::calcWallSlipPndDiffusiveTerm(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
		double ni = Particles->pndi[i];
		//if(Particles->particleType[i] == PSystem->fluid && ni > PSystem->epsilonZero) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true && ni > PSystem->epsilonZero) {
	//	if(Particles->particleType[i] == PSystem->fluid) {
			double DivV = 0.0;
			//double Pi = Particles->press[i];
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
	//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2

			// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
			double Rref_i[9], normaliw[3], normalMod2;
			// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
			normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
			normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];
			if(normalMod2 > PSystem->epsilonZero) {
				double normalMod = sqrt(normalMod2);
				normaliw[0] = normaliw[0]/normalMod;
				normaliw[1] = normaliw[1]/normalMod;
				normaliw[2] = normaliw[2]/normalMod;
			}
			else {
				normaliw[0] = 0;
				normaliw[1] = 0;
				normaliw[2] = 0;
			}

		    //  Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
		    Rref_i[0] = 1.0 - 2.0*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2.0*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2.0*normaliw[0]*normaliw[2];
			Rref_i[3] = 0.0 - 2.0*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2.0*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2.0*normaliw[1]*normaliw[2];
			Rref_i[6] = 0.0 - 2.0*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2.0*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2.0*normaliw[2]*normaliw[2];

			// Mirror particle velocity vi' = Rref_i * vi
	      	double velMirrorXi = (Rref_i[0]*velXi + Rref_i[1]*velYi + Rref_i[2]*velZi);
			double velMirrorYi = (Rref_i[3]*velXi + Rref_i[4]*velYi + Rref_i[5]*velZi);
			double velMirrorZi = (Rref_i[6]*velXi + Rref_i[7]*velYi + Rref_i[8]*velZi);

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

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < PSystem->reS2 && dstimj2 < PSystem->reS2 && dstij2 < dstimj2) {
						if(j != i) {
							double dst = sqrt(dstimj2);
							double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
							double vijx = Particles->vel[j*3  ]-velMirrorXi;
							double vijy = Particles->vel[j*3+1]-velMirrorYi;
							double vijz = Particles->vel[j*3+2]-velMirrorZi;
							DivV += (PSystem->dim/PSystem->pndSmallZero)*(Particles->pndi[j]/ni)*(vijx*v0imj+vijy*v1imj+vijz*v2imj)*wS/dstimj2;
						}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}

			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			Particles->sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
			
			if(dstimi2 < PSystem->reS2) {
				double dst = sqrt(dstimi2);
				double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
				double vijx = velXi-velMirrorXi;
				double vijy = velYi-velMirrorYi;
				double vijz = velZi-velMirrorZi;
				DivV += (PSystem->dim/PSystem->pndSmallZero)*(Particles->pndi[i]/ni)*(vijx*v0imi+vijy*v1imi+vijz*v2imi)*wS/dstimi2;
		  	}

			Particles->pndAuxVar1[i] += -Particles->pndi[i]*PSystem->timeStep*DivV;
			//Particles->pndAuxVar1[i] = PSystem->pndSmallZero*(1.0+PSystem->timeStep*(-DivV+Di*flagDi));
		}
	}
//#pragma omp parallel for
//	for(int i=0; i<Particles->numParticles; i++) {
//	if(Particles->particleType[i] == PSystem->fluid) {
//		Particles->pndi[i] += Particles->pndAuxVar1[i];
//		Particles->pndAuxVar1[i]=0.0;
//	}}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Calculates PND based on continuity equation including a diffusive term. Contribution of the closest polygon wall. No-slip boundary condition.
void MpsPndNeigh::calcWallNoSlipPndDiffusiveTerm(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
		double ni = Particles->pndi[i];
		//if(Particles->particleType[i] == PSystem->fluid && ni > PSystem->epsilonZero) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true && ni > PSystem->epsilonZero) {
	//	if(Particles->particleType[i] == PSystem->fluid) {
			double DivV = 0.0;
			//double Pi = Particles->press[i];
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
	//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2

			// Inverse matrix Rinv_i = - I
			double Rinv_i[9], Rref_i[9], normaliw[3], normalMod2;
		    // normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		    normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
		    normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];
		    if(normalMod2 > PSystem->epsilonZero) {
		    	double normalMod = sqrt(normalMod2);
		    	normaliw[0] = normaliw[0]/normalMod;
		    	normaliw[1] = normaliw[1]/normalMod;
		    	normaliw[2] = normaliw[2]/normalMod;
		    }
		    else {
		    	normaliw[0] = 0;
		    	normaliw[1] = 0;
		    	normaliw[2] = 0;
		    }

		    //  Inverse transformation matrix Rinv_i = - I
		    Rinv_i[0] = -1.0; Rinv_i[1] =  0.0; Rinv_i[2] =  0.0;
			Rinv_i[3] =  0.0; Rinv_i[4] = -1.0; Rinv_i[5] =  0.0;
			Rinv_i[6] =  0.0; Rinv_i[7] =  0.0; Rinv_i[8] = -1.0;

		    //  Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
		    Rref_i[0] = 1.0 - 2.0*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2.0*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2.0*normaliw[0]*normaliw[2];
			Rref_i[3] = 0.0 - 2.0*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2.0*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2.0*normaliw[1]*normaliw[2];
			Rref_i[6] = 0.0 - 2.0*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2.0*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2.0*normaliw[2]*normaliw[2];

			double viwall[3], vtil[3];
			// Wall velocity (0 if fixed)
			viwall[0]=viwall[1]=viwall[2]=0.0;

			if(Particles->nearMeshType[i] == meshType::FORCED) {
				viwall[0] = PSystem->uniformVelWall[0];
				viwall[1] = PSystem->uniformVelWall[1];
				viwall[2] = PSystem->uniformVelWall[2];
			}

			// normal_iwall*v_iwall
			double dotnv = normaliw[0]*viwall[0] + normaliw[1]*viwall[1] + normaliw[2]*viwall[2];
			// vtil = vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}
			vtil[0] = velXi - 2.0*(viwall[0] - dotnv*normaliw[0]);
			vtil[1] = velYi - 2.0*(viwall[1] - dotnv*normaliw[1]);
			vtil[2] = velZi - 2.0*(viwall[2] - dotnv*normaliw[2]);
			// Mirror particle velocity vi' = Rinv_i * [vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}] 
	      	double velMirrorXi = (Rinv_i[0]*vtil[0] + Rinv_i[1]*vtil[1] + Rinv_i[2]*vtil[2]);
			double velMirrorYi = (Rinv_i[3]*vtil[0] + Rinv_i[4]*vtil[1] + Rinv_i[5]*vtil[2]);
			double velMirrorZi = (Rinv_i[6]*vtil[0] + Rinv_i[7]*vtil[1] + Rinv_i[8]*vtil[2]);

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

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < PSystem->reS2 && dstimj2 < PSystem->reS2 && dstij2 < dstimj2) {
						if(j != i) {
							double dst = sqrt(dstimj2);
							double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
							double vijx = -(Particles->vel[j*3  ]-velMirrorXi);
							double vijy = -(Particles->vel[j*3+1]-velMirrorYi);
							double vijz = -(Particles->vel[j*3+2]-velMirrorZi);
							// Refelected rij' = Rref_i * ri'j
	      					double v0m = (Rref_i[0]*v0imj + Rref_i[1]*v1imj + Rref_i[2]*v2imj);
							double v1m = (Rref_i[3]*v0imj + Rref_i[4]*v1imj + Rref_i[5]*v2imj);
							double v2m = (Rref_i[6]*v0imj + Rref_i[7]*v1imj + Rref_i[8]*v2imj);
							DivV += (PSystem->dim/PSystem->pndSmallZero)*(Particles->pndi[j]/ni)*(vijx*v0m+vijy*v1m+vijz*v2m)*wS/dstimj2;
						}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
					}
			}}}

			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			Particles->sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
			
			if(dstimi2 < PSystem->reS2) {
				double dst = sqrt(dstimi2);
				double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
				double vijx = -(velXi-velMirrorXi);
				double vijy = -(velYi-velMirrorYi);
				double vijz = -(velZi-velMirrorZi);
				// Refelected rij' = Rref_i * ri'j
				double v0m = (Rref_i[0]*v0imi + Rref_i[1]*v1imi + Rref_i[2]*v2imi);
				double v1m = (Rref_i[3]*v0imi + Rref_i[4]*v1imi + Rref_i[5]*v2imi);
				double v2m = (Rref_i[6]*v0imi + Rref_i[7]*v1imi + Rref_i[8]*v2imi);
				DivV += (PSystem->dim/PSystem->pndSmallZero)*(Particles->pndi[i]/ni)*(vijx*v0m+vijy*v1m+vijz*v2m)*wS/dstimi2;
		  	}

			Particles->pndAuxVar1[i] += -Particles->pndi[i]*PSystem->timeStep*DivV;
			//Particles->pndAuxVar1[i] = PSystem->pndSmallZero*(1.0+PSystem->timeStep*(-DivV+Di*flagDi));
		}
	}
//#pragma omp parallel for
//	for(int i=0; i<Particles->numParticles; i++) {
//	if(Particles->particleType[i] == PSystem->fluid) {
//		Particles->pndi[i] += Particles->pndAuxVar1[i];
//		Particles->pndAuxVar1[i]=0.0;
//	}}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Updates diffusive PND
void MpsPndNeigh::updatePnd(MpsParticle *Particles) {
#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
	//		if(Particles->particleType[i] == PSystem->fluid)
			Particles->pndi[i] = Particles->pndAuxVar1[i];
//		else
//		{
//			Particles->pndi[i] = Particles->pndAuxVar1[i];
			/*
			double mi;
			if(PTYPE[i] == 1) 
				mi = PSystem->DNS_FL1;
			else 
				mi = PSystem->DNS_FL2;
			if(mpsType == calcPressType::EXPLICIT)
				Particles->pndi[i] = PSystem->pndSmallZero*(Particles->press[i]/(mi*PSystem->coeffPressWCMPS)+1);
			else if(mpsType == calcPressType::WEAKLY)
				Particles->pndi[i] = PSystem->pndSmallZero*pow(Particles->press[i]*PSystem->gamma/(mi*PSystem->coeffPressWCMPS)+1,PSystem->gamma);
				*/
//		}
		Particles->pndAuxVar1[i]=0.0;
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Mean diffusive PND at wall, dummy or free-surface particles
void MpsPndNeigh::meanPndParticlesWallDummySurface(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		// Set variables to zero
		Particles->pndAuxVar1[i] = 0.0;
		Particles->pndAuxVar2[i] = 0.0;
//	if(Particles->particleType[i] == PSystem->wall) {
		if(Particles->particleType[i] == PSystem->wall || Particles->particleBC[i] == PSystem->surface) {
//	if(Particles->particleBC[i] == PSystem->surface) {
			double PNDup = 0.0;
			double PNDdo = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
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
					if(dstij2 < PSystem->reS2 && (dstij2 < dstimj2 || PSystem->wallType == boundaryWallType::PARTICLE)) {
						if(j != i) {
							double dst = sqrt(dstij2);
							double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
							PNDup += Particles->pndi[j]*wS;
							PNDdo += wS;
						}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			Particles->pndAuxVar1[i] = PNDup;
			Particles->pndAuxVar2[i] = PNDdo;
			//Particles->pndAuxVar1[i] = PSystem->pndSmallZero*(1.0+PSystem->timeStep*(-DivV+Di*flagDi));
	//	}}}
		}
	}
#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
//	if(Particles->particleType[i] == PSystem->wall) {
		if(Particles->particleType[i] == PSystem->wall || Particles->particleBC[i] == PSystem->surface) {
//	if(Particles->particleBC[i] == PSystem->surface) {
		// Prevent PNDdo = 0
//		if(Particles->numNeigh[i] < 1)
			if(Particles->pndAuxVar2[i] < PSystem->epsilonZero)
				Particles->pndi[i] = Particles->pndAuxVar1[i];
			else
				Particles->pndi[i] = Particles->pndAuxVar1[i]/(Particles->pndAuxVar2[i]);
//			Particles->pndi[i] = Particles->pndAuxVar1[i]/(Particles->pndAuxVar2[i] + 0.01*PSystem->reS2/4.0);
		}
//	}
		Particles->pndAuxVar1[i]=0.0;
		Particles->pndAuxVar2[i]=0.0;
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Mean PND computed as the sum of the weight function. Contribution of the neighboring particles.
void MpsPndNeigh::meanPnd(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		double PNDup = 0.0;
		double PNDdo = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
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
				if(dstij2 < PSystem->reS2 && (dstij2 < dstimj2 || PSystem->wallType == boundaryWallType::PARTICLE)) {
					if(j != i) {
						double dst = sqrt(dstij2);
						double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
						PNDup += Particles->pndi[j]*wS;
						PNDdo += wS;
					}
				}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// Sum to polygon wall contribution
		Particles->pndAuxVar1[i] += PNDup;
		Particles->pndAuxVar2[i] += PNDdo;
		//Particles->pndAuxVar1[i] = PSystem->pndSmallZero*(1.0+PSystem->timeStep*(-DivV+Di*flagDi));
	}
#pragma omp parallel for
for(int i=0; i<Particles->numParticles; i++) {
//	if(Particles->particleType[i] == PSystem->fluid) {
	// Prevent PNDdo = 0
		if(Particles->numNeigh[i] < 1) {
			Particles->pndi[i] = Particles->pndAuxVar1[i];
		}
		else {
			Particles->pndi[i] = Particles->pndAuxVar1[i]/Particles->pndAuxVar2[i];
		}
		// Set variables to zero
		Particles->pndAuxVar1[i] = 0.0;
		Particles->pndAuxVar2[i] = 0.0;
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Mean PND computed as the sum of the weight function. Contribution of the closest polygon wall.
void MpsPndNeigh::meanWallPnd(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
		// Set variables to zero
		Particles->pndAuxVar1[i] = 0.0;
		Particles->pndAuxVar2[i] = 0.0;
//	if(Particles->particleType[i] == PSystem->fluid) {
		if(Particles->particleNearWall[i] == true) {
			double PNDup = 0.0;
			double PNDdo = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
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

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < PSystem->reS2 && dstimj2 < PSystem->reS2 && dstij2 < dstimj2) {
						if(j != i) {
							double dst = sqrt(dstimj2);
							double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
							PNDup += Particles->pndi[j]*wS;
							PNDdo += wS;
						}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}

			}}}
			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			Particles->sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
			
			if(dstimi2 < PSystem->reS2) {
				double dst = sqrt(dstimi2);
				double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
				PNDup += Particles->pndi[i]*wS;
				PNDdo += wS;
			}
			Particles->pndAuxVar1[i] = PNDup;
			Particles->pndAuxVar2[i] = PNDdo;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Mean PND computed as the sum of the weight function. Contribution of the fluid neighboring particles. // CHANGED
void MpsPndNeigh::meanNeighFluidPnd(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		double PNDup = 0.0;
		double PNDdo = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
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
				if(dstij2 < PSystem->reS2 && (dstij2 < dstimj2 || PSystem->wallType == boundaryWallType::PARTICLE)) {
					if(j != i && Particles->particleType[j] == PSystem->fluid) {
						double dst = sqrt(dstij2);
						double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
						PNDup += Particles->pndki[j]*wS;
						PNDdo += wS;
					}
				}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		if(PNDdo > PSystem->epsilonZero) {
			Particles->pndski[i] = PNDup/PNDdo;
		}
		else {
			Particles->pndski[i] = PNDup;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}