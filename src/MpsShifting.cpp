// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <experimental/filesystem> 	///< numeric_limits
#include <iostream>					///< cout
#include "MpsShifting.h"

using namespace std;

// Constructor declaration
MpsShifting::MpsShifting()
{
}
// Destructor declaration
MpsShifting::~MpsShifting()
{
}


// Adjustment of particle velocity
// Improvements for accuracy and stability in a weakly-compressible particle method
// https://www.sciencedirect.com/science/article/pii/S0045793016302250
void MpsShifting::calcShifting(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
	if(Particles->particleType[i] == PSystem->fluid) {
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double duXi = 0.0;	double duYi = 0.0;	double duZi = 0.0;
		//double ni = 0.0;
		
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
					double dw = Particles->delWeight(dst, PSystem->reS, PSystem->weightType);
					duXi += dw*(Particles->vel[j*3  ]-velXi);
					duYi += dw*(Particles->vel[j*3+1]-velYi);
					duZi += dw*(Particles->vel[j*3+2]-velZi);
					//double w = Particles->weight(dst, r, PSystem->weightType);
					//ni += w;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

//		if(PSystem->wallType == boundaryWallType::POLYGON) {
//			if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero)
//		}
//		else if(PSystem->wallType == boundaryWallType::PARTICLE) 
//			if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
//		}
//		if(pndSmall[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
		if(Particles->particleBC[i] == PSystem->inner) {
			Particles->vel[i*3  ] -= PSystem->coeffShifting1*duXi;
			Particles->vel[i*3+1] -= PSystem->coeffShifting1*duYi;
			Particles->vel[i*3+2] -= PSystem->coeffShifting1*duZi;
		}
		// else {
		// 	double Inn[9], duAux[3];
		// 	// I - nxn
		// 	Inn[0] = 1.0 - Particles->normal[i*3  ]*Particles->normal[i*3  ]; Inn[1] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+1]; Inn[2] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+2];
		// 	Inn[3] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3  ]; Inn[4] = 1.0 - Particles->normal[i*3+1]*Particles->normal[i*3+1]; Inn[5] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3+2];
		// 	Inn[6] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3  ]; Inn[7] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3+1]; Inn[8] = 1.0 - Particles->normal[i*3+2]*Particles->normal[i*3+2];
		// 	// (I - nxn)dr
		// 	duAux[0] = Inn[0]*duXi + Inn[1]*duYi + Inn[2]*duZi;
		// 	duAux[1] = Inn[3]*duXi + Inn[4]*duYi + Inn[5]*duZi;
		// 	duAux[2] = Inn[6]*duXi + Inn[7]*duYi + Inn[8]*duZi;
		// 	Particles->vel[i*3  ] -= PSystem->coeffShifting1*duAux[0];
		// 	Particles->vel[i*3+1] -= PSystem->coeffShifting1*duAux[1];
		// 	Particles->vel[i*3+2] -= PSystem->coeffShifting1*duAux[2];
		// }
	}}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Normal vector on the fluid
// An accurate and stable multiphase moving particle semi-implicit method based on a corrective matrix for all particle interaction models
// https://onlinelibrary.wiley.com/doi/full/10.1002/nme.5844
void MpsShifting::calcNormalParticles(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double dr_ix = 0.0;	double dr_iy = 0.0;	double dr_iz = 0.0;
		
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
					dr_ix += v0ij*wS/dst;
					dr_iy += v1ij*wS/dst;
					dr_iz += v2ij*wS/dst;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		Particles->normal[i*3  ] = - dr_ix/PSystem->pndSmallZero;
		Particles->normal[i*3+1] = - dr_iy/PSystem->pndSmallZero;
		Particles->normal[i*3+2] = - dr_iz/PSystem->pndSmallZero;
		
		/*
		if(PSystem->wallType == boundaryWallType::PARTICLE)
		{
			// Normalize
			double norm2 = dr_ix*dr_ix + dr_iy*dr_iy + dr_iz*dr_iz;
			if(norm2 > 0.0) {
				double norm = sqrt(norm2);
				Particles->normal[i*3  ] /= norm;
				Particles->normal[i*3+1] /= norm;
				Particles->normal[i*3+2] /= norm;
			}
			else {
				Particles->normal[i*3  ] = 0.0;
				Particles->normal[i*3+1] = 0.0;
				Particles->normal[i*3+2] = 0.0;
			}
		}
		*/
	}
#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Adjustment of particle velocity (Polygon wall)
void MpsShifting::calcWallShifting(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
	//if(Particles->particleType[i] == PSystem->fluid) {
	if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
		double duXi = 0.0;	double duYi = 0.0;	double duZi = 0.0;

		// No-slip

		// Inverse matrix Rinv_i = - I
		double Rinv_i[9], normaliw[3], normalMod2;
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

		// Inverse transformation matrix Rinv_i = - I
		Rinv_i[0] = -1.0; Rinv_i[1] =  0.0; Rinv_i[2] =  0.0;
		Rinv_i[3] =  0.0; Rinv_i[4] = -1.0; Rinv_i[5] =  0.0;
		Rinv_i[6] =  0.0; Rinv_i[7] =  0.0; Rinv_i[8] = -1.0;

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
					double dw = Particles->delWeight(dst, PSystem->reS, PSystem->weightType);
					duXi += dw*(Particles->vel[j*3  ]-velMirrorXi);
					duYi += dw*(Particles->vel[j*3+1]-velMirrorYi);
					duZi += dw*(Particles->vel[j*3+2]-velMirrorZi);
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// Add "i" contribution ("i" is a neighbor of "mirror i")
		double v0imi, v1imi, v2imi, dstimi2;
		Particles->sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
	  	
		if(dstimi2 < PSystem->reS2) {
			double dst = sqrt(dstimi2);
			double dw = Particles->delWeight(dst, PSystem->reS, PSystem->weightType);
			duXi += dw*(velXi-velMirrorXi);
			duYi += dw*(velYi-velMirrorYi);
			duZi += dw*(velZi-velMirrorZi);
		}
	
//		if(PSystem->wallType == boundaryWallType::POLYGON) {
//			if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero)
//		}
//		else if(PSystem->wallType == boundaryWallType::PARTICLE) 
//			if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
//		}
//		if(pndSmall[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
		if(Particles->particleBC[i] == PSystem->inner) {
			double dux = Rinv_i[0]*duXi + Rinv_i[1]*duYi + Rinv_i[2]*duZi;
			double duy = Rinv_i[3]*duXi + Rinv_i[4]*duYi + Rinv_i[5]*duZi;
			double duz = Rinv_i[6]*duXi + Rinv_i[7]*duYi + Rinv_i[8]*duZi;
			Particles->vel[i*3  ] -= PSystem->coeffShifting1*dux;
			Particles->vel[i*3+1] -= PSystem->coeffShifting1*duy;
			Particles->vel[i*3+2] -= PSystem->coeffShifting1*duz;
		}
		// else {
		// 	double Inn[9], duAux[3];
		// 	// I - nxn
		// 	Inn[0] = 1.0 - Particles->normal[i*3  ]*Particles->normal[i*3  ]; Inn[1] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+1]; Inn[2] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+2];
		// 	Inn[3] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3  ]; Inn[4] = 1.0 - Particles->normal[i*3+1]*Particles->normal[i*3+1]; Inn[5] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3+2];
		// 	Inn[6] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3  ]; Inn[7] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3+1]; Inn[8] = 1.0 - Particles->normal[i*3+2]*Particles->normal[i*3+2];
		// 	double dux = Rinv_i[0]*duXi + Rinv_i[1]*duYi + Rinv_i[2]*duZi;
		//	double duy = Rinv_i[3]*duXi + Rinv_i[4]*duYi + Rinv_i[5]*duZi;
		//	double duz = Rinv_i[6]*duXi + Rinv_i[7]*duYi + Rinv_i[8]*duZi;
		// 	// (I - nxn)dr
		// 	duAux[0] = Inn[0]*dux + Inn[1]*duy + Inn[2]*duz;
		// 	duAux[1] = Inn[3]*dux + Inn[4]*duy + Inn[5]*duz;
		// 	duAux[2] = Inn[6]*dux + Inn[7]*duy + Inn[8]*duz;
		// 	Particles->vel[i*3  ] -= PSystem->coeffShifting1*duAux[0];
		// 	Particles->vel[i*3+1] -= PSystem->coeffShifting1*duAux[1];
		// 	Particles->vel[i*3+2] -= PSystem->coeffShifting1*duAux[2];
		// }
	}}
}

// Normal vector on the fluid (Polygon wall)
void MpsShifting::calcWallNormalParticles(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
	//if(Particles->particleType[i] == PSystem->fluid) {
	if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double dr_ix = 0.0;	double dr_iy = 0.0;	double dr_iz = 0.0;
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

		//  Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
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
					double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
					dr_ix += v0imj*wS/dst;
					dr_iy += v1imj*wS/dst;
					dr_iz += v2imj*wS/dst;
				}}
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
			dr_ix += v0imi*wS/dst;
			dr_iy += v1imi*wS/dst;
			dr_iz += v2imi*wS/dst;
		}

		double drx = Rref_i[0]*dr_ix + Rref_i[1]*dr_iy + Rref_i[2]*dr_iz;
		double dry = Rref_i[3]*dr_ix + Rref_i[4]*dr_iy + Rref_i[5]*dr_iz;
		double drz = Rref_i[6]*dr_ix + Rref_i[7]*dr_iy + Rref_i[8]*dr_iz;
		
		Particles->normal[i*3  ] += - drx/PSystem->pndSmallZero;
		Particles->normal[i*3+1] += - dry/PSystem->pndSmallZero;
		Particles->normal[i*3+2] += - drz/PSystem->pndSmallZero;

/*
		// Normalize
		double norm2 = Particles->normal[i*3  ]*Particles->normal[i*3  ] + Particles->normal[i*3+1]*Particles->normal[i*3+1] + Particles->normal[i*3+2]*Particles->normal[i*3+2];
		if(norm2 > 0.0) {
			double norm = sqrt(norm2);
			Particles->normal[i*3  ] /= norm;
			Particles->normal[i*3+1] /= norm;
			Particles->normal[i*3+2] /= norm;
		}
		else {
			Particles->normal[i*3  ] = 0.0;
			Particles->normal[i*3+1] = 0.0;
			Particles->normal[i*3+2] = 0.0;
		}
*/	
	}}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Concentration and Gradient of concentration. Adjustment of particle position.
void MpsShifting::calcConcAndConcGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	// Concentration
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		Particles->concentration[i] = 0.0;
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
					Particles->concentration[i] += wS;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Add PND due Wall polygon
		Particles->concentration[i] += Particles->pndWallContribution[i];

		Particles->concentration[i] /= PSystem->pndSmallZero;
	}
	// Gradient of concentration
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
	if(Particles->particleType[i] == PSystem->fluid) {
		Particles->gradConcentration[i*3  ] = 0.0;	Particles->gradConcentration[i*3+1] = 0.0;	Particles->gradConcentration[i*3+2] = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double conc_i = Particles->concentration[i];
		
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
					Particles->gradConcentration[i*3  ] += (conc_i + Particles->concentration[j])*v0ij*wS/dstij2;
					Particles->gradConcentration[i*3+1] += (conc_i + Particles->concentration[j])*v1ij*wS/dstij2;
					Particles->gradConcentration[i*3+2] += (conc_i + Particles->concentration[j])*v2ij*wS/dstij2;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		//PSystem->coeffPressGrad = -PSystem->dim/PSystem->pndGradientZero
		Particles->gradConcentration[3*i  ] *= -PSystem->coeffPressGrad;
		Particles->gradConcentration[3*i+1] *= -PSystem->coeffPressGrad;
		Particles->gradConcentration[3*i+2] *= -PSystem->coeffPressGrad;

//		if(PSystem->wallType == boundaryWallType::POLYGON) {
//			if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero)
//		}
//		else if(PSystem->wallType == boundaryWallType::PARTICLE) 
//			if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
//		}
//		if(pndSmall[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
		if(Particles->particleBC[i] == PSystem->inner) {
			// PSystem->coeffShifting2 = PSystem->coefA*PSystem->partDist*PSystem->partDist*PSystem->cflNumber*PSystem->machNumber;	// Coefficient used to adjust velocity
			Particles->pos[i*3  ] -= PSystem->coeffShifting2*Particles->gradConcentration[3*i  ];
			Particles->pos[i*3+1] -= PSystem->coeffShifting2*Particles->gradConcentration[3*i+1];
			Particles->pos[i*3+2] -= PSystem->coeffShifting2*Particles->gradConcentration[3*i+2];
		}
		/*
		else
		{
			double Inn[9], drAux[3];
			// I - nxn
			Inn[0] = 1.0 - Particles->normal[i*3  ]*Particles->normal[i*3  ]; Inn[1] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+1]; Inn[2] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+2];
			Inn[3] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3  ]; Inn[4] = 1.0 - Particles->normal[i*3+1]*Particles->normal[i*3+1]; Inn[5] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3+2];
			Inn[6] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3  ]; Inn[7] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3+1]; Inn[8] = 1.0 - Particles->normal[i*3+2]*Particles->normal[i*3+2];
			// (I - nxn)dr
			drAux[0] = Inn[0]*Particles->gradConcentration[3*i  ] + Inn[1]*Particles->gradConcentration[3*i+1] + Inn[2]*Particles->gradConcentration[3*i+2];
			drAux[1] = Inn[3]*Particles->gradConcentration[3*i  ] + Inn[4]*Particles->gradConcentration[3*i+1] + Inn[5]*Particles->gradConcentration[3*i+2];
			drAux[2] = Inn[6]*Particles->gradConcentration[3*i  ] + Inn[7]*Particles->gradConcentration[3*i+1] + Inn[8]*Particles->gradConcentration[3*i+2];
			// PSystem->coeffShifting2 = PSystem->coefA*PSystem->partDist*PSystem->partDist*PSystem->cflNumber*PSystem->machNumber;	// Coefficient used to adjust Velocity
			Particles->pos[i*3  ] -= PSystem->coeffShifting2*drAux[0];
			Particles->pos[i*3+1] -= PSystem->coeffShifting2*drAux[1];
			Particles->pos[i*3+2] -= PSystem->coeffShifting2*drAux[2];
		}
		*/
		// else {
		// 	double Inn[9], duAux[3];
		// 	// I - nxn
		// 	Inn[0] = 1.0 - Particles->normal[i*3  ]*Particles->normal[i*3  ]; Inn[1] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+1]; Inn[2] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+2];
		// 	Inn[3] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3  ]; Inn[4] = 1.0 - Particles->normal[i*3+1]*Particles->normal[i*3+1]; Inn[5] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3+2];
		// 	Inn[6] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3  ]; Inn[7] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3+1]; Inn[8] = 1.0 - Particles->normal[i*3+2]*Particles->normal[i*3+2];
		// 	// (I - nxn)dr
		// 	duAux[0] = Inn[0]*duXi + Inn[1]*duYi + Inn[2]*duZi;
		// 	duAux[1] = Inn[3]*duXi + Inn[4]*duYi + Inn[5]*duZi;
		// 	duAux[2] = Inn[6]*duXi + Inn[7]*duYi + Inn[8]*duZi;
		// 	Particles->vel[i*3  ] -= PSystem->coeffShifting1*duAux[0];
		// 	Particles->vel[i*3+1] -= PSystem->coeffShifting1*duAux[1];
		// 	Particles->vel[i*3+2] -= PSystem->coeffShifting1*duAux[2];
		// }
	}}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Concentration and Gradient of concentration. Adjustment of particle position. (Polygon wall)
void MpsShifting::calcWallConcAndConcGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	// Gradient of concentration due Polygon wall
	//int nPartNearMesh = partNearMesh.size();
	// double VolumeForce = pow(PSystem->partDist,PSystem->dim);
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
	//if(Particles->particleType[i] == PSystem->fluid) {
	if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double drX = 0.0;			double drY = 0.0;			double drZ = 0.0;
		double gradCiWallX = 0.0;			double gradCiWallY = 0.0;			double gradCiWallZ = 0.0;
		double conc_i = Particles->concentration[i];

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

		//  Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
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

					drX += (conc_i + Particles->concentration[j])*v0imj*wS/dstimj2;
					drY += (conc_i + Particles->concentration[j])*v1imj*wS/dstimj2;
					drZ += (conc_i + Particles->concentration[j])*v2imj*wS/dstimj2;
				}}
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

			drX += (conc_i + conc_i)*v0imi*wS/dstimi2;
			drY += (conc_i + conc_i)*v1imi*wS/dstimi2;
			drZ += (conc_i + conc_i)*v2imi*wS/dstimi2;
		}

		// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/noGrad)
		// Original
//		acc[i*3  ] += (PSystem->relaxPress*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*PSystem->coeffPressGrad - rpsForce[0])*invDns[partType::FLUID];
//		acc[i*3+1] += (PSystem->relaxPress*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*PSystem->coeffPressGrad - rpsForce[1])*invDns[partType::FLUID];
//		acc[i*3+2] += (PSystem->relaxPress*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*PSystem->coeffPressGrad - rpsForce[2])*invDns[partType::FLUID];
		// Modified
		gradCiWallX = -(Rref_i[0]*drX + Rref_i[1]*drY + Rref_i[2]*drZ)*PSystem->coeffPressGrad;
		gradCiWallY = -(Rref_i[3]*drX + Rref_i[4]*drY + Rref_i[5]*drZ)*PSystem->coeffPressGrad;
		gradCiWallZ = -(Rref_i[6]*drX + Rref_i[7]*drY + Rref_i[8]*drZ)*PSystem->coeffPressGrad;

		Particles->gradConcentration[i*3  ] += gradCiWallX;
		Particles->gradConcentration[i*3+1] += gradCiWallY;
		Particles->gradConcentration[i*3+2] += gradCiWallZ;

		if(Particles->particleBC[i] == PSystem->inner) {
			// PSystem->coeffShifting2 = PSystem->coefA*PSystem->partDist*PSystem->partDist*PSystem->cflNumber*PSystem->machNumber;	// Coefficient used to adjust velocity
			Particles->pos[i*3  ] -= PSystem->coeffShifting2*gradCiWallX;
			Particles->pos[i*3+1] -= PSystem->coeffShifting2*gradCiWallY;
			Particles->pos[i*3+2] -= PSystem->coeffShifting2*gradCiWallZ;
		}
	}}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Normal vector on the fluid
void MpsShifting::calcNormalConcentration(MpsParticle *Particles) {
#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		double norm2GradCi = Particles->gradConcentration[3*i]*Particles->gradConcentration[3*i] + Particles->gradConcentration[3*i+1]*Particles->gradConcentration[3*i+1] + Particles->gradConcentration[3*i+2]*Particles->gradConcentration[3*i+2];
		
		if(norm2GradCi > 0.0) {
			double norm = sqrt(norm2GradCi);
			Particles->normal[i*3  ] = -Particles->gradConcentration[i*3  ]/norm;
			Particles->normal[i*3+1] = -Particles->gradConcentration[i*3+1]/norm;
			Particles->normal[i*3+2] = -Particles->gradConcentration[i*3+2]/norm;
		}
		else {
			Particles->normal[i*3  ] = 0.0;
			Particles->normal[i*3+1] = 0.0;
			Particles->normal[i*3+2] = 0.0;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}