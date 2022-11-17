// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <iostream>		///< cout
#include <experimental/filesystem> 	///< numeric_limits
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

#define ADJ_VEL		// Adjustment of velocity
#define CONC_GRAD 	// Gradient of concentration

// Adjustment of particle velocity
// Improvements for accuracy and stability in a weakly-compressible particle method
// https://www.sciencedirect.com/science/article/pii/S0045793016302250
void MpsShifting::calcShifting(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		Particles->dVelShift[i*3  ] = 0.0;	Particles->dVelShift[i*3+1] = 0.0;	Particles->dVelShift[i*3+2] = 0.0;
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
#ifdef ADJ_VEL
						duXi += dw*(Particles->vel[j*3  ]-velXi);
						duYi += dw*(Particles->vel[j*3+1]-velYi);
						duZi += dw*(Particles->vel[j*3+2]-velZi);
#else
						duXi += dw*v0ij;
						duYi += dw*v1ij;
						duZi += dw*v2ij;
#endif
						//double w = Particles->weight(dst, r, PSystem->weightType);
						//ni += w;
					}}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}

			Particles->dVelShift[i*3  ] = PSystem->coeffShifting1 * duXi;
			Particles->dVelShift[i*3+1] = PSystem->coeffShifting1 * duYi;
			Particles->dVelShift[i*3+2] = PSystem->coeffShifting1 * duZi;

			// if(PSystem->wallType == boundaryWallType::POLYGON) {
			// 	if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero)
			// }
			// else if(PSystem->wallType == boundaryWallType::PARTICLE) 
			// 	if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
			// }
			// if(pndSmall[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
			// if(Particles->particleBC[i] == PSystem->inner) {
			// 	double du_x = PSystem->coeffShifting1 * duXi;
			// 	double du_y = PSystem->coeffShifting1 * duYi;
			// 	double du_z = PSystem->coeffShifting1 * duZi;
			// 	double duMod2 = du_x * du_x + du_y * du_y + du_z * du_z;

			// 	// double dr_x = du_x * PSystem->timeStep;
			// 	// double dr_y = du_y * PSystem->timeStep;
			// 	// double dr_z = du_z * PSystem->timeStep;
			// 	// double drMod2 = dr_x * dr_x + dr_y * dr_y + dr_z * dr_z;

			// 	if (duMod2 > PSystem->epsilonZero) {
			// 		double duMod = sqrt(duMod2);
			// 		double duMin = min(duMod, PSystem->maxDu);
			// 		Particles->dVelShift[i*3  ] = - duMin*du_x/duMod;
			// 		Particles->dVelShift[i*3+1] = - duMin*du_y/duMod;
			// 		Particles->dVelShift[i*3+2] = - duMin*du_z/duMod;
			// 		// Particles->vel[i*3  ] -= duMin*du_x/duMod;
			// 		// Particles->vel[i*3+1] -= duMin*du_y/duMod;
			// 		// Particles->vel[i*3+2] -= duMin*du_z/duMod;
			// 	}

			// 	// if (drMod2 > PSystem->epsilonZero) {
			// 	// 	double drMod = sqrt(drMod2);
			// 	// 	double drMin = min(drMod, PSystem->maxDr);
			// 	// 	Particles->dPosShift[i*3  ] = - drMin*dr_x/drMod;
			// 	// 	Particles->dPosShift[i*3+1] = - drMin*dr_y/drMod;
			// 	// 	Particles->dPosShift[i*3+2] = - drMin*dr_z/drMod;
			// 	// 	Particles->pos[i*3  ] += Particles->dPosShift[i*3  ];
			// 	// 	Particles->pos[i*3+1] += Particles->dPosShift[i*3+1];
			// 	// 	Particles->pos[i*3+2] += Particles->dPosShift[i*3+2];
			// 	// }
			// }
			// else {
			// 	double Inn[9];
			// 	// I - nxn
			// 	Inn[0] = 1.0 - Particles->normal[i*3  ]*Particles->normal[i*3  ]; Inn[1] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+1]; Inn[2] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+2];
			// 	Inn[3] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3  ]; Inn[4] = 1.0 - Particles->normal[i*3+1]*Particles->normal[i*3+1]; Inn[5] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3+2];
			// 	Inn[6] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3  ]; Inn[7] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3+1]; Inn[8] = 1.0 - Particles->normal[i*3+2]*Particles->normal[i*3+2];
			// 	// (I - nxn)dr
			// 	double du_x = PSystem->coeffShifting1 * (Inn[0]*duXi + Inn[1]*duYi + Inn[2]*duZi);
			// 	double du_y = PSystem->coeffShifting1 * (Inn[3]*duXi + Inn[4]*duYi + Inn[5]*duZi);
			// 	double du_z = PSystem->coeffShifting1 * (Inn[6]*duXi + Inn[7]*duYi + Inn[8]*duZi);
			// 	double duMod2 = du_x * du_x + du_y * du_y + du_z * du_z;

			// 	if (duMod2 > PSystem->epsilonZero) {
			// 		double duMod = sqrt(duMod2);
			// 		double duMin = min(duMod, PSystem->maxDu);
			// 		Particles->vel[i*3  ] -= duMin*du_x/duMod;
			// 		Particles->vel[i*3+1] -= duMin*du_y/duMod;
			// 		Particles->vel[i*3+2] -= duMin*du_z/duMod;
			// 	}
			// }
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Adjustment of particle velocity (Polygon wall)
// Improvements for accuracy and stability in a weakly-compressible particle method
// https://www.sciencedirect.com/science/article/pii/S0045793016302250
void MpsShifting::calcWallShifting(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//  Inverse transformation matrix Rinv_i = - I
	double Rinv_i[9];
	Rinv_i[0] = -1.0; Rinv_i[1] =  0.0; Rinv_i[2] =  0.0;
	Rinv_i[3] =  0.0; Rinv_i[4] = -1.0; Rinv_i[5] =  0.0;
	Rinv_i[6] =  0.0; Rinv_i[7] =  0.0; Rinv_i[8] = -1.0;
	//int nPartNearMesh = partNearMesh.size();
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		//if(Particles->particleType[i] == PSystem->fluid) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double duXi = 0.0;	double duYi = 0.0;	double duZi = 0.0;

			// No-slip

			double normaliw[3], normalMod2;
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
#ifdef ADJ_VEL
						duXi += dw*(Particles->vel[j*3  ]-velMirrorXi);
						duYi += dw*(Particles->vel[j*3+1]-velMirrorYi);
						duZi += dw*(Particles->vel[j*3+2]-velMirrorZi);
#else
						duXi += dw*v0imj;
						duYi += dw*v1imj;
						duZi += dw*v2imj;
#endif
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
#ifdef ADJ_VEL
				duXi += dw*(velXi-velMirrorXi);
				duYi += dw*(velYi-velMirrorYi);
				duZi += dw*(velZi-velMirrorZi);
#else
				duXi += dw*v0imi;
				duYi += dw*v1imi;
				duZi += dw*v2imi;
#endif
			}

			Particles->dVelShift[i*3  ] += PSystem->coeffShifting1 * (Rinv_i[0]*duXi + Rinv_i[1]*duYi + Rinv_i[2]*duZi);
			Particles->dVelShift[i*3+1] += PSystem->coeffShifting1 * (Rinv_i[3]*duXi + Rinv_i[4]*duYi + Rinv_i[5]*duZi);
			Particles->dVelShift[i*3+2] += PSystem->coeffShifting1 * (Rinv_i[6]*duXi + Rinv_i[7]*duYi + Rinv_i[8]*duZi);

			// // if(PSystem->wallType == boundaryWallType::POLYGON) {
			// // 	if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero)
			// // }
			// // else if(PSystem->wallType == boundaryWallType::PARTICLE) 
			// // 	if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
			// // }
			// // if(pndSmall[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
			// if(Particles->particleBC[i] == PSystem->inner) {
			// 	double du_x = PSystem->coeffShifting1 * (Rinv_i[0]*duXi + Rinv_i[1]*duYi + Rinv_i[2]*duZi);
			// 	double du_y = PSystem->coeffShifting1 * (Rinv_i[3]*duXi + Rinv_i[4]*duYi + Rinv_i[5]*duZi);
			// 	double du_z = PSystem->coeffShifting1 * (Rinv_i[6]*duXi + Rinv_i[7]*duYi + Rinv_i[8]*duZi);
			// 	double duMod2 = du_x * du_x + du_y * du_y + du_z * du_z;

			// 	if (duMod2 > PSystem->epsilonZero) {
			// 		double duMod = sqrt(duMod2);
			// 		double duMin = min(duMod, PSystem->maxDu);
			// 		Particles->dVelShift[i*3  ] = - duMin*du_x/duMod;
			// 		Particles->dVelShift[i*3+1] = - duMin*du_y/duMod;
			// 		Particles->dVelShift[i*3+2] = - duMin*du_z/duMod;
			// 		// Particles->vel[i*3  ] -= duMin*du_x/duMod;
			// 		// Particles->vel[i*3+1] -= duMin*du_y/duMod;
			// 		// Particles->vel[i*3+2] -= duMin*du_z/duMod;
			// 	}
			// }
			// else {
			// 	double Inn[9];
			// 	// I - nxn
			// 	Inn[0] = 1.0 - Particles->normal[i*3  ]*Particles->normal[i*3  ]; Inn[1] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+1]; Inn[2] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+2];
			// 	Inn[3] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3  ]; Inn[4] = 1.0 - Particles->normal[i*3+1]*Particles->normal[i*3+1]; Inn[5] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3+2];
			// 	Inn[6] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3  ]; Inn[7] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3+1]; Inn[8] = 1.0 - Particles->normal[i*3+2]*Particles->normal[i*3+2];
			// 	double dux = Rinv_i[0]*duXi + Rinv_i[1]*duYi + Rinv_i[2]*duZi;
			// 	double duy = Rinv_i[3]*duXi + Rinv_i[4]*duYi + Rinv_i[5]*duZi;
			// 	double duz = Rinv_i[6]*duXi + Rinv_i[7]*duYi + Rinv_i[8]*duZi;
			// 	// (I - nxn)dr
			// 	double du_x = PSystem->coeffShifting1 * (Inn[0]*dux + Inn[1]*duy + Inn[2]*duz);
			// 	double du_y = PSystem->coeffShifting1 * (Inn[3]*dux + Inn[4]*duy + Inn[5]*duz);
			// 	double du_z = PSystem->coeffShifting1 * (Inn[6]*dux + Inn[7]*duy + Inn[8]*duz);
			// 	double duMod2 = du_x * du_x + du_y * du_y + du_z * du_z;

			// 	if (duMod2 > PSystem->epsilonZero) {
			// 		double duMod = sqrt(duMod2);
			// 		double duMin = min(duMod, PSystem->maxDu);
			// 		Particles->vel[i*3  ] -= duMin*du_x/duMod;
			// 		Particles->vel[i*3+1] -= duMin*du_y/duMod;
			// 		Particles->vel[i*3+2] -= duMin*du_z/duMod;
			// 	}
			// }
		}
	}
}

// Adjustment of particle velocity based on the shifting type 1
// Improvements for accuracy and stability in a weakly-compressible particle method
// https://www.sciencedirect.com/science/article/pii/S0045793016302250
void MpsShifting::updateVelocity(MpsParticleSystem *PSystem, MpsParticle *Particles) {
#pragma omp parallel for
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		// if(Particles->particleBC[i] == PSystem->inner) {
		if(Particles->particleType[i] == PSystem->fluid) {
		// if(Particles->particleBC[i] == PSystem->inner && Particles->numNeigh[i] > 29) {
		
			double du_x = Particles->dVelShift[i*3  ];
			double du_y = Particles->dVelShift[i*3+1];
			double du_z = Particles->dVelShift[i*3+2];
			double duMod2 = du_x * du_x + du_y * du_y + du_z * du_z;

			if (duMod2 > PSystem->epsilonZero) {
				double duMod = sqrt(duMod2);
				double duMin = min(duMod, PSystem->maxDu);
				Particles->dVelShift[i*3  ] = - duMin*du_x/duMod;
				Particles->dVelShift[i*3+1] = - duMin*du_y/duMod;
				Particles->dVelShift[i*3+2] = - duMin*du_z/duMod;
#ifdef ADJ_VEL
				Particles->vel[i*3  ] -= duMin*du_x/duMod;
				Particles->vel[i*3+1] -= duMin*du_y/duMod;
				Particles->vel[i*3+2] -= duMin*du_z/duMod;
#else
				Particles->pos[i*3  ] -= duMin*du_x/duMod;
				Particles->pos[i*3+1] -= duMin*du_y/duMod;
				Particles->pos[i*3+2] -= duMin*du_z/duMod;
#endif
			}

			// double dr_x = PSystem->timeStep * Particles->dVelShift[i*3  ];
			// double dr_y = PSystem->timeStep * Particles->dVelShift[i*3+1];
			// double dr_z = PSystem->timeStep * Particles->dVelShift[i*3+2];
			// double drMod2 = dr_x * dr_x + dr_y * dr_y + dr_z * dr_z;

			// if (drMod2 > PSystem->epsilonZero) {
			// 	double drMod = sqrt(drMod2);
			// 	double drMin = min(drMod, PSystem->maxDr);
			// 	Particles->dPosShift[i*3  ] = - drMin*dr_x/drMod;
			// 	Particles->dPosShift[i*3+1] = - drMin*dr_y/drMod;
			// 	Particles->dPosShift[i*3+2] = - drMin*dr_z/drMod;
			// 	Particles->pos[i*3  ] += Particles->dPosShift[i*3  ];
			// 	Particles->pos[i*3+1] += Particles->dPosShift[i*3+1];
			// 	Particles->pos[i*3+2] += Particles->dPosShift[i*3+2];

			// 	// Particles->vel[i*3  ] += Particles->dPosShift[i*3  ];
			// 	// Particles->vel[i*3+1] += Particles->dPosShift[i*3+1];
			// 	// Particles->vel[i*3+2] += Particles->dPosShift[i*3+2];
			// }
		}
	}

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
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
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
		
		
		// if(PSystem->wallType == boundaryWallType::PARTICLE)
		// {
		// 	// Normalize
		// 	double norm2 = dr_ix*dr_ix + dr_iy*dr_iy + dr_iz*dr_iz;
		// 	if(norm2 > 0.0) {
		// 		double norm = sqrt(norm2);
		// 		Particles->normal[i*3  ] /= norm;
		// 		Particles->normal[i*3+1] /= norm;
		// 		Particles->normal[i*3+2] /= norm;
		// 	}
		// 	else {
		// 		Particles->normal[i*3  ] = 0.0;
		// 		Particles->normal[i*3+1] = 0.0;
		// 		Particles->normal[i*3+2] = 0.0;
		// 	}
		// }
		
	}
#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Normal vector on the fluid (Polygon wall)
void MpsShifting::calcWallNormalParticles(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		// if(Particles->particleType[i] == PSystem->fluid) {
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

	
			// // Normalize
			// double norm2 = Particles->normal[i*3  ]*Particles->normal[i*3  ] + Particles->normal[i*3+1]*Particles->normal[i*3+1] + Particles->normal[i*3+2]*Particles->normal[i*3+2];
			// if(norm2 > 0.0) {
			// 	double norm = sqrt(norm2);
			// 	Particles->normal[i*3  ] /= norm;
			// 	Particles->normal[i*3+1] /= norm;
			// 	Particles->normal[i*3+2] /= norm;
			// }
			// else {
			// 	Particles->normal[i*3  ] = 0.0;
			// 	Particles->normal[i*3+1] = 0.0;
			// 	Particles->normal[i*3+2] = 0.0;
			// }
		
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Concentration. Adjustment of particle position.
void MpsShifting::calcConcentration(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	// Concentration
#pragma omp parallel for schedule(dynamic,64)
	for(int ip=0; ip<Particles->numParticles; ip++) {

		int i = Particles->particleID[ip];
		// if(Particles->particleType[i] != PSystem->inOutflowParticle) {
		if(Particles->particleType[i] == PSystem->fluid) {
			double conc_i = 0.0;
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
						conc_i += wS;
					}}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}

			// Add PND due Wall polygon
			// Particles->concentration[i] += Particles->pndWallContribution[i];

			Particles->concentration[i] = conc_i / PSystem->pndSmallZero;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Concentration (Polygon wall)
void MpsShifting::calcWallConcentration(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		//if(Particles->particleType[i] == PSystem->fluid) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {

			// Add PND due Wall polygon
			// Particles->concentration[i] += Particles->pndWallContribution[i] / PSystem->pndSmallZero;

			double conc_i = 0.0;
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
						conc_i += wS;
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
				conc_i += wS;
			}

			Particles->concentration[i] += conc_i / PSystem->pndSmallZero;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}


// Gradient of concentration. Adjustment of particle position.
void MpsShifting::calcConcentrationGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	// Concentration
// #pragma omp parallel for schedule(dynamic,64)
// 	for(int ip=0; ip<Particles->numParticles; ip++) {
// 		int i = Particles->particleID[ip];
// 		Particles->concentration[i] = 0.0;
// 		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
// 		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		
// 		int ix, iy, iz;
// 		Buckets->bucketCoordinates(ix, iy, iz, posXi, posYi, posZi, PSystem);
// 		int minZ = (iz-1)*((int)(PSystem->dim-2.0)); int maxZ = (iz+1)*((int)(PSystem->dim-2.0));
// 		for(int jz=minZ;jz<=maxZ;jz++) {
// 		for(int jy=iy-1;jy<=iy+1;jy++) {
// 		for(int jx=ix-1;jx<=ix+1;jx++) {
// 			int jb = jz*PSystem->numBucketsXY + jy*PSystem->numBucketsX + jx;
// 			int j = Particles->firstParticleInBucket[jb];
// 			if(j == -1) continue;
// 			double plx, ply, plz;
// 			Particles->getPeriodicLengths(jb, plx, ply, plz, PSystem);
// 			while(true) {
// 				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
// 				// Particle square distance r_ij^2 = (Xj - Xi_temporary_position)^2
// 				Particles->sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2, plx, ply, plz);
// 				// Mirror particle square distance r_imj^2 = (Xj - Xim_temporary_position)^2
// 				Particles->sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2, plx, ply, plz);

// 				// If j is inside the neighborhood of i and 
// 				// is not at the same side of im (avoid real j in the virtual neihborhood)
// 				if(dstij2 < PSystem->reS2 && (dstij2 < dstimj2 || PSystem->wallType == boundaryWallType::PARTICLE)) {
// 				if(j != i) {
// 					double dst = sqrt(dstij2);
// 					double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
// 					Particles->concentration[i] += wS;
// 				}}
// 				j = Particles->nextParticleInSameBucket[j];
// 				if(j == -1) break;
// 			}
// 		}}}

// 		// Add PND due Wall polygon
// 		Particles->concentration[i] += Particles->pndWallContribution[i];

// 		Particles->concentration[i] /= PSystem->pndSmallZero;
// 	}
	// Gradient of concentration
#pragma omp parallel for schedule(dynamic,64)
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		if(Particles->particleType[i] == PSystem->fluid) {
			Particles->gradConcentration[i*3  ] = 0.0;	Particles->gradConcentration[i*3+1] = 0.0;	Particles->gradConcentration[i*3+2] = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			double conc_i = Particles->concentration[i];

			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			
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
#ifdef CONC_GRAD
						double conc_ij = conc_i + Particles->concentration[j];
						// Particles->gradConcentration[i*3  ] += conc_ij*v0ij*wS/dstij2;
						// Particles->gradConcentration[i*3+1] += conc_ij*v1ij*wS/dstij2;
						// Particles->gradConcentration[i*3+2] += conc_ij*v2ij*wS/dstij2;
						
						// if(conc_ij > 2.0) conc_ij -= 2.0;
						Particles->gradConcentration[i*3  ] += conc_ij*v0ij*wS/dstij2;
						Particles->gradConcentration[i*3+1] += conc_ij*v1ij*wS/dstij2;
						Particles->gradConcentration[i*3+2] += conc_ij*v2ij*wS/dstij2;

						// Particles->gradConcentration[i*3  ] += conc_ij*(Particles->vel[j*3  ]-velXi)*wS/dstij2;
						// Particles->gradConcentration[i*3+1] += conc_ij*(Particles->vel[j*3+1]-velYi)*wS/dstij2;
						// Particles->gradConcentration[i*3+2] += conc_ij*(Particles->vel[j*3+2]-velZi)*wS/dstij2;
#else
						Particles->gradConcentration[i*3  ] += v0ij*wS/dstij2;
						Particles->gradConcentration[i*3+1] += v1ij*wS/dstij2;
						Particles->gradConcentration[i*3+2] += v2ij*wS/dstij2;
#endif
					}}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			//PSystem->coeffPressGrad = -PSystem->dim/PSystem->pndGradientZero
			Particles->gradConcentration[i*3  ] *= -PSystem->coeffPressGrad;
			Particles->gradConcentration[i*3+1] *= -PSystem->coeffPressGrad;
			Particles->gradConcentration[i*3+2] *= -PSystem->coeffPressGrad;

			// if(PSystem->wallType == boundaryWallType::POLYGON) {
			// 	if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero)
			// }
			// else if(PSystem->wallType == boundaryWallType::PARTICLE) 
			// 	if(Particles->pndi[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
			// }
			// if(pndSmall[i] > pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero)
			
			// if(Particles->particleBC[i] == PSystem->inner) {
			// 	// PSystem->coeffShifting2 = PSystem->coefA*PSystem->partDist*PSystem->partDist*PSystem->cflNumber*PSystem->machNumber;	// Coefficient used to adjust velocity
			// 	Particles->pos[i*3  ] -= PSystem->coeffShifting2*Particles->gradConcentration[i*3  ];
			// 	Particles->pos[i*3+1] -= PSystem->coeffShifting2*Particles->gradConcentration[i*3+1];
			// 	Particles->pos[i*3+2] -= PSystem->coeffShifting2*Particles->gradConcentration[i*3+2];
			// }
			
			// else
			// {
			// 	double Inn[9], drAux[3];
			// 	// I - nxn
			// 	Inn[0] = 1.0 - Particles->normal[i*3  ]*Particles->normal[i*3  ]; Inn[1] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+1]; Inn[2] = 0.0 - Particles->normal[i*3  ]*Particles->normal[i*3+2];
			// 	Inn[3] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3  ]; Inn[4] = 1.0 - Particles->normal[i*3+1]*Particles->normal[i*3+1]; Inn[5] = 0.0 - Particles->normal[i*3+1]*Particles->normal[i*3+2];
			// 	Inn[6] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3  ]; Inn[7] = 0.0 - Particles->normal[i*3+2]*Particles->normal[i*3+1]; Inn[8] = 1.0 - Particles->normal[i*3+2]*Particles->normal[i*3+2];
			// 	// (I - nxn)dr
			// 	drAux[0] = Inn[0]*Particles->gradConcentration[i*3  ] + Inn[1]*Particles->gradConcentration[i*3+1] + Inn[2]*Particles->gradConcentration[i*3+2];
			// 	drAux[1] = Inn[3]*Particles->gradConcentration[i*3  ] + Inn[4]*Particles->gradConcentration[i*3+1] + Inn[5]*Particles->gradConcentration[i*3+2];
			// 	drAux[2] = Inn[6]*Particles->gradConcentration[i*3  ] + Inn[7]*Particles->gradConcentration[i*3+1] + Inn[8]*Particles->gradConcentration[i*3+2];
			// 	// PSystem->coeffShifting2 = PSystem->coefA*PSystem->partDist*PSystem->partDist*PSystem->cflNumber*PSystem->machNumber;	// Coefficient used to adjust Velocity
			// 	Particles->pos[i*3  ] -= PSystem->coeffShifting2*drAux[0];
			// 	Particles->pos[i*3+1] -= PSystem->coeffShifting2*drAux[1];
			// 	Particles->pos[i*3+2] -= PSystem->coeffShifting2*drAux[2];
			// }
			
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
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Gradient of concentration. Adjustment of particle position. (Polygon wall)
void MpsShifting::calcWallConcentrationGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//  Inverse transformation matrix Rinv_i = - I
	double Rinv_i[9];
	Rinv_i[0] = -1.0; Rinv_i[1] =  0.0; Rinv_i[2] =  0.0;
	Rinv_i[3] =  0.0; Rinv_i[4] = -1.0; Rinv_i[5] =  0.0;
	Rinv_i[6] =  0.0; Rinv_i[7] =  0.0; Rinv_i[8] = -1.0;

	// Gradient of concentration due Polygon wall
	//int nPartNearMesh = partNearMesh.size();
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		//if(Particles->particleType[i] == PSystem->fluid) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			double drX = 0.0;			double drY = 0.0;			double drZ = 0.0;
			double gradCiWallX = 0.0;	double gradCiWallY = 0.0;	double gradCiWallZ = 0.0;
			double conc_i = Particles->concentration[i];

			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];

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


			// No-slip

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
#ifdef CONC_GRAD
						double conc_ij = conc_i + Particles->concentration[j];
						// drX += conc_ij*v0imj*wS/dstimj2;
						// drY += conc_ij*v1imj*wS/dstimj2;
						// drZ += conc_ij*v2imj*wS/dstimj2;

						// if(conc_ij > 2.0) conc_ij -= 2.0;
						drX += conc_ij*v0imj*wS/dstimj2;
						drY += conc_ij*v1imj*wS/dstimj2;
						drZ += conc_ij*v2imj*wS/dstimj2;

						// drX += conc_ij*(Particles->vel[j*3  ]-velMirrorXi)*wS/dstimj2;
						// drY += conc_ij*(Particles->vel[j*3+1]-velMirrorYi)*wS/dstimj2;
						// drZ += conc_ij*(Particles->vel[j*3+2]-velMirrorZi)*wS/dstimj2;
#else
						drX += v0imj*wS/dstimj2;
						drY += v1imj*wS/dstimj2;
						drZ += v2imj*wS/dstimj2;
#endif
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
#ifdef CONC_GRAD
				double conc_ii = conc_i + conc_i;
				// drX += conc_ii*v0imi*wS/dstimi2;
				// drY += conc_ii*v1imi*wS/dstimi2;
				// drZ += conc_ii*v2imi*wS/dstimi2;

				// if(conc_ii > 2.0) conc_ii -= 2.0;
				drX += conc_ii*v0imi*wS/dstimi2;
				drY += conc_ii*v1imi*wS/dstimi2;
				drZ += conc_ii*v2imi*wS/dstimi2;

				// drX += conc_ii*(Particles->vel[i*3  ]-velMirrorXi)*wS/dstimi2;
				// drY += conc_ii*(Particles->vel[i*3+1]-velMirrorYi)*wS/dstimi2;
				// drZ += conc_ii*(Particles->vel[i*3+2]-velMirrorZi)*wS/dstimi2;
#else
				drX += v0imi*wS/dstimi2;
				drY += v1imi*wS/dstimi2;
				drZ += v2imi*wS/dstimi2;
#endif
			}

			// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/noGrad)
			// Original
			// acc[i*3  ] += (PSystem->relaxPress*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*PSystem->coeffPressGrad - rpsForce[0])*invDns[partType::FLUID];
			// acc[i*3+1] += (PSystem->relaxPress*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*PSystem->coeffPressGrad - rpsForce[1])*invDns[partType::FLUID];
			// acc[i*3+2] += (PSystem->relaxPress*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*PSystem->coeffPressGrad - rpsForce[2])*invDns[partType::FLUID];
			// Modified
			gradCiWallX = -(Rref_i[0]*drX + Rref_i[1]*drY + Rref_i[2]*drZ)*PSystem->coeffPressGrad;
			gradCiWallY = -(Rref_i[3]*drX + Rref_i[4]*drY + Rref_i[5]*drZ)*PSystem->coeffPressGrad;
			gradCiWallZ = -(Rref_i[6]*drX + Rref_i[7]*drY + Rref_i[8]*drZ)*PSystem->coeffPressGrad;

			Particles->gradConcentration[i*3  ] += gradCiWallX;
			Particles->gradConcentration[i*3+1] += gradCiWallY;
			Particles->gradConcentration[i*3+2] += gradCiWallZ;

			// if(Particles->particleBC[i] == PSystem->inner) {
			// 	// PSystem->coeffShifting2 = PSystem->coefA*PSystem->partDist*PSystem->partDist*PSystem->cflNumber*PSystem->machNumber;	// Coefficient used to adjust velocity
			// 	Particles->pos[i*3  ] -= PSystem->coeffShifting2*gradCiWallX;
			// 	Particles->pos[i*3+1] -= PSystem->coeffShifting2*gradCiWallY;
			// 	Particles->pos[i*3+2] -= PSystem->coeffShifting2*gradCiWallZ;
			// }
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Adjustment of particle position based on the Gradient of concentration
void MpsShifting::updatePosition(MpsParticleSystem *PSystem, MpsParticle *Particles) {
#pragma omp parallel for
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		Particles->dPosShift[i*3  ] = 0.0;	Particles->dPosShift[i*3+1] = 0.0;	Particles->dPosShift[i*3+2] = 0.0;
		// if(Particles->particleBC[i] == PSystem->inner) {
		if(Particles->particleType[i] == PSystem->fluid) {
		// if(Particles->particleBC[i] == PSystem->inner && Particles->numNeigh[i] > 29) {
			// PSystem->coeffShifting2 = PSystem->coefA*PSystem->partDist*PSystem->partDist*PSystem->cflNumber*PSystem->machNumber;	// Coefficient used to adjust velocity
#ifdef CONC_GRAD
			double dr_x = PSystem->coeffShifting2 * Particles->gradConcentration[i*3  ];
			double dr_y = PSystem->coeffShifting2 * Particles->gradConcentration[i*3+1];
			double dr_z = PSystem->coeffShifting2 * Particles->gradConcentration[i*3+2];

			// double dr_x = PSystem->reS * PSystem->coeffShifting2 * Particles->gradConcentration[i*3  ] / (PSystem->partDist * PSystem->partDist);
			// double dr_y = PSystem->reS * PSystem->coeffShifting2 * Particles->gradConcentration[i*3+1] / (PSystem->partDist * PSystem->partDist);
			// double dr_z = PSystem->reS * PSystem->coeffShifting2 * Particles->gradConcentration[i*3+2] / (PSystem->partDist * PSystem->partDist);
#else
			double dr_x = (0.1 / PSystem->dim) * PSystem->partDist * PSystem->partDist * Particles->gradConcentration[i*3  ];
			double dr_y = (0.1 / PSystem->dim) * PSystem->partDist * PSystem->partDist * Particles->gradConcentration[i*3+1];
			double dr_z = (0.1 / PSystem->dim) * PSystem->partDist * PSystem->partDist * Particles->gradConcentration[i*3+2];
#endif
			double drMod2 = dr_x * dr_x + dr_y * dr_y + dr_z * dr_z;

			if (drMod2 > PSystem->epsilonZero) {
				double drMod = sqrt(drMod2);
				double drMin = min(drMod, PSystem->maxDr);
				Particles->dPosShift[i*3  ] = - drMin*dr_x/drMod;
				Particles->dPosShift[i*3+1] = - drMin*dr_y/drMod;
				Particles->dPosShift[i*3+2] = - drMin*dr_z/drMod;
				Particles->pos[i*3  ] += Particles->dPosShift[i*3  ];
				Particles->pos[i*3+1] += Particles->dPosShift[i*3+1];
				Particles->pos[i*3+2] += Particles->dPosShift[i*3+2];

				// Particles->vel[i*3  ] += Particles->dPosShift[i*3  ];
				// Particles->vel[i*3+1] += Particles->dPosShift[i*3+1];
				// Particles->vel[i*3+2] += Particles->dPosShift[i*3+2];
			}
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Gradient of velocity. Adjustment of particle position. (Polygon wall)
void MpsShifting::calcVelGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		if(Particles->particleType[i] == PSystem->fluid) {
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			double MC[9];
			MC[0] = Particles->correcMatrixRow1[i*3];	MC[1] = Particles->correcMatrixRow1[i*3+1];	MC[2] = Particles->correcMatrixRow1[i*3+2];
			MC[3] = Particles->correcMatrixRow2[i*3];	MC[4] = Particles->correcMatrixRow2[i*3+1];	MC[5] = Particles->correcMatrixRow2[i*3+2];
			MC[6] = Particles->correcMatrixRow3[i*3];	MC[7] = Particles->correcMatrixRow3[i*3+1];	MC[8] = Particles->correcMatrixRow3[i*3+2];
			double gradUxx = 0.0;	double gradUxy = 0.0;	double gradUxz = 0.0;
			double gradUyx = 0.0;	double gradUyy = 0.0;	double gradUyz = 0.0;
			double gradUzx = 0.0;	double gradUzy = 0.0;	double gradUzz = 0.0;

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
						double vijx = Particles->vel[j*3  ] - velXi;	
						double vijy = Particles->vel[j*3+1] - velYi;	
						double vijz = Particles->vel[j*3+2] - velZi;
						double invDstij2 = 1.0/dstij2;

						if(PSystem->gradientCorrection == false) {
							gradUxx += vijx*v0ij*wS*invDstij2;
							gradUxy += vijx*v1ij*wS*invDstij2;
							gradUxz += vijx*v2ij*wS*invDstij2;
							
							gradUyx += vijy*v0ij*wS*invDstij2;
							gradUyy += vijy*v1ij*wS*invDstij2;
							gradUyz += vijy*v2ij*wS*invDstij2;
							
							gradUzx += vijz*v0ij*wS*invDstij2;
							gradUzy += vijz*v1ij*wS*invDstij2;
							gradUzz += vijz*v2ij*wS*invDstij2;
						}
						else {
							double v0ijC = (v0ij*MC[0] + v1ij*MC[1] + v2ij*MC[2]);
							double v1ijC = (v0ij*MC[3] + v1ij*MC[4] + v2ij*MC[5]);
							double v2ijC = (v0ij*MC[6] + v1ij*MC[7] + v2ij*MC[8]);
							gradUxx += vijx*v0ijC*wS*invDstij2;
							gradUxy += vijx*v1ijC*wS*invDstij2;
							gradUxz += vijx*v2ijC*wS*invDstij2;
							
							gradUyx += vijy*v0ijC*wS*invDstij2;
							gradUyy += vijy*v1ijC*wS*invDstij2;
							gradUyz += vijy*v2ijC*wS*invDstij2;
							
							gradUzx += vijz*v0ijC*wS*invDstij2;
							gradUzy += vijz*v1ijC*wS*invDstij2;
							gradUzz += vijz*v2ijC*wS*invDstij2;
						}
					}}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}

			// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/noGrad)
			Particles->velGradientRow1[i*3  ] = -PSystem->coeffPressGrad*gradUxx;
			Particles->velGradientRow1[i*3+1] = -PSystem->coeffPressGrad*gradUxy;
			Particles->velGradientRow1[i*3+2] = -PSystem->coeffPressGrad*gradUxz;
			Particles->velGradientRow2[i*3  ] = -PSystem->coeffPressGrad*gradUyx;
			Particles->velGradientRow2[i*3+1] = -PSystem->coeffPressGrad*gradUyy;
			Particles->velGradientRow2[i*3+2] = -PSystem->coeffPressGrad*gradUyz;
			Particles->velGradientRow3[i*3  ] = -PSystem->coeffPressGrad*gradUzx;
			Particles->velGradientRow3[i*3+1] = -PSystem->coeffPressGrad*gradUzy;
			Particles->velGradientRow3[i*3+2] = -PSystem->coeffPressGrad*gradUzz;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Gradient of velocity (Polygon wall) - Free-slip
void MpsShifting::calcWallSlipVelGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {

	//int nPartNearMesh = partNearMesh.size();
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
	// if(Particles->particleType[i] == PSystem->fluid) {
		double ni = Particles->pndi[i];
		if(ni < PSystem->epsilonZero) continue;	// Avoid division by zero

		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			double MC[9];
			MC[0] = Particles->correcMatrixRow1[i*3];	MC[1] = Particles->correcMatrixRow1[i*3+1];	MC[2] = Particles->correcMatrixRow1[i*3+2];
			MC[3] = Particles->correcMatrixRow2[i*3];	MC[4] = Particles->correcMatrixRow2[i*3+1];	MC[5] = Particles->correcMatrixRow2[i*3+2];
			MC[6] = Particles->correcMatrixRow3[i*3];	MC[7] = Particles->correcMatrixRow3[i*3+1];	MC[8] = Particles->correcMatrixRow3[i*3+2];
			double gradUxx = 0.0;	double gradUxy = 0.0;	double gradUxz = 0.0;
			double gradUyx = 0.0;	double gradUyy = 0.0;	double gradUyz = 0.0;
			double gradUzx = 0.0;	double gradUzy = 0.0;	double gradUzz = 0.0;

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
							double vimjx = Particles->vel[j*3  ]-velMirrorXi;
							double vimjy = Particles->vel[j*3+1]-velMirrorYi;
							double vimjz = Particles->vel[j*3+2]-velMirrorZi;

							// Refelected rij' = Rref_i * ri'j
							double v0ijm = (Rref_i[0]*v0imj + Rref_i[1]*v1imj + Rref_i[2]*v2imj);
							double v1ijm = (Rref_i[3]*v0imj + Rref_i[4]*v1imj + Rref_i[5]*v2imj);
							double v2ijm = (Rref_i[6]*v0imj + Rref_i[7]*v1imj + Rref_i[8]*v2imj);
							double invDstimj2 = 1.0/dstimj2;

							// Refelected uij' = Rref_i * ui'j
							double vijmx = (Rref_i[0]*vimjx + Rref_i[1]*vimjy + Rref_i[2]*vimjz);
							double vijmy = (Rref_i[3]*vimjx + Rref_i[4]*vimjy + Rref_i[5]*vimjz);
							double vijmz = (Rref_i[6]*vimjx + Rref_i[7]*vimjy + Rref_i[8]*vimjz);

							if(PSystem->gradientCorrection == false) {
								gradUxx += vijmx*v0ijm*wS*invDstimj2;
								gradUxy += vijmx*v1ijm*wS*invDstimj2;
								gradUxz += vijmx*v2ijm*wS*invDstimj2;
								
								gradUyx += vijmy*v0ijm*wS*invDstimj2;
								gradUyy += vijmy*v1ijm*wS*invDstimj2;
								gradUyz += vijmy*v2ijm*wS*invDstimj2;
								
								gradUzx += vijmz*v0ijm*wS*invDstimj2;
								gradUzy += vijmz*v1ijm*wS*invDstimj2;
								gradUzz += vijmz*v2ijm*wS*invDstimj2;
							}
							else {
								double v0ijmC = (v0ijm*MC[0] + v1ijm*MC[1] + v2ijm*MC[2]);
								double v1ijmC = (v0ijm*MC[3] + v1ijm*MC[4] + v2ijm*MC[5]);
								double v2ijmC = (v0ijm*MC[6] + v1ijm*MC[7] + v2ijm*MC[8]);
								gradUxx += vijmx*v0ijmC*wS*invDstimj2;
								gradUxy += vijmx*v1ijmC*wS*invDstimj2;
								gradUxz += vijmx*v2ijmC*wS*invDstimj2;
								
								gradUyx += vijmy*v0ijmC*wS*invDstimj2;
								gradUyy += vijmy*v1ijmC*wS*invDstimj2;
								gradUyz += vijmy*v2ijmC*wS*invDstimj2;
								
								gradUzx += vijmz*v0ijmC*wS*invDstimj2;
								gradUzy += vijmz*v1ijmC*wS*invDstimj2;
								gradUzz += vijmz*v2ijmC*wS*invDstimj2;
							}
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
				double vimix = velXi-velMirrorXi;
				double vimiy = velYi-velMirrorYi;
				double vimiz = velZi-velMirrorZi;

				// Refelected rii' = Rref_i * ri'i
				double v0iim = (Rref_i[0]*v0imi + Rref_i[1]*v1imi + Rref_i[2]*v2imi);
				double v1iim = (Rref_i[3]*v0imi + Rref_i[4]*v1imi + Rref_i[5]*v2imi);
				double v2iim = (Rref_i[6]*v0imi + Rref_i[7]*v1imi + Rref_i[8]*v2imi);
				double invDstimi2 = 1.0/dstimi2;

				// Refelected uii' = Rref_i * ui'i
				double viimx = (Rref_i[0]*vimix + Rref_i[1]*vimiy + Rref_i[2]*vimiz);
				double viimy = (Rref_i[3]*vimix + Rref_i[4]*vimiy + Rref_i[5]*vimiz);
				double viimz = (Rref_i[6]*vimix + Rref_i[7]*vimiy + Rref_i[8]*vimiz);

				if(PSystem->gradientCorrection == false) {
					gradUxx += viimx*v0iim*wS*invDstimi2;
					gradUxy += viimx*v1iim*wS*invDstimi2;
					gradUxz += viimx*v2iim*wS*invDstimi2;
					
					gradUyx += viimy*v0iim*wS*invDstimi2;
					gradUyy += viimy*v1iim*wS*invDstimi2;
					gradUyz += viimy*v2iim*wS*invDstimi2;
					
					gradUzx += viimz*v0iim*wS*invDstimi2;
					gradUzy += viimz*v1iim*wS*invDstimi2;
					gradUzz += viimz*v2iim*wS*invDstimi2;
				}
				else {
					double v0iimC = (v0iim*MC[0] + v1iim*MC[1] + v2iim*MC[2]);
					double v1iimC = (v0iim*MC[3] + v1iim*MC[4] + v2iim*MC[5]);
					double v2iimC = (v0iim*MC[6] + v1iim*MC[7] + v2iim*MC[8]);
					gradUxx += viimx*v0iimC*wS*invDstimi2;
					gradUxy += viimx*v1iimC*wS*invDstimi2;
					gradUxz += viimx*v2iimC*wS*invDstimi2;
					
					gradUyx += viimy*v0iimC*wS*invDstimi2;
					gradUyy += viimy*v1iimC*wS*invDstimi2;
					gradUyz += viimy*v2iimC*wS*invDstimi2;
					
					gradUzx += viimz*v0iimC*wS*invDstimi2;
					gradUzy += viimz*v1iimC*wS*invDstimi2;
					gradUzz += viimz*v2iimC*wS*invDstimi2;
				}
			}

			// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/noGrad)
			Particles->velGradientRow1[i*3  ] += -PSystem->coeffPressGrad*gradUxx;
			Particles->velGradientRow1[i*3+1] += -PSystem->coeffPressGrad*gradUxy;
			Particles->velGradientRow1[i*3+2] += -PSystem->coeffPressGrad*gradUxz;
			Particles->velGradientRow2[i*3  ] += -PSystem->coeffPressGrad*gradUyx;
			Particles->velGradientRow2[i*3+1] += -PSystem->coeffPressGrad*gradUyy;
			Particles->velGradientRow2[i*3+2] += -PSystem->coeffPressGrad*gradUyz;
			Particles->velGradientRow3[i*3  ] += -PSystem->coeffPressGrad*gradUzx;
			Particles->velGradientRow3[i*3+1] += -PSystem->coeffPressGrad*gradUzy;
			Particles->velGradientRow3[i*3+2] += -PSystem->coeffPressGrad*gradUzz;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Gradient of velocity (Polygon wall) - No-slip
void MpsShifting::calcWallNoSlipVelGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {

	//  Inverse transformation matrix Rinv_i = - I
	double Rinv_i[9];
	Rinv_i[0] = -1.0; Rinv_i[1] =  0.0; Rinv_i[2] =  0.0;
	Rinv_i[3] =  0.0; Rinv_i[4] = -1.0; Rinv_i[5] =  0.0;
	Rinv_i[6] =  0.0; Rinv_i[7] =  0.0; Rinv_i[8] = -1.0;
	
	//int nPartNearMesh = partNearMesh.size();
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		double ni = Particles->pndi[i];
		if(ni < PSystem->epsilonZero) continue;	// Avoid division by zero

		// if(Particles->particleType[i] == PSystem->fluid && ni > PSystem->epsilonZero) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
		// if(Particles->particleType[i] == PSystem->fluid) {
			double DivV = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			double MC[9];
			MC[0] = Particles->correcMatrixRow1[i*3];	MC[1] = Particles->correcMatrixRow1[i*3+1];	MC[2] = Particles->correcMatrixRow1[i*3+2];
			MC[3] = Particles->correcMatrixRow2[i*3];	MC[4] = Particles->correcMatrixRow2[i*3+1];	MC[5] = Particles->correcMatrixRow2[i*3+2];
			MC[6] = Particles->correcMatrixRow3[i*3];	MC[7] = Particles->correcMatrixRow3[i*3+1];	MC[8] = Particles->correcMatrixRow3[i*3+2];
			double gradUxx = 0.0;	double gradUxy = 0.0;	double gradUxz = 0.0;
			double gradUyx = 0.0;	double gradUyy = 0.0;	double gradUyz = 0.0;
			double gradUzx = 0.0;	double gradUzy = 0.0;	double gradUzz = 0.0;

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
							double vimjx = -(Particles->vel[j*3  ]-velMirrorXi);
							double vimjy = -(Particles->vel[j*3+1]-velMirrorYi);
							double vimjz = -(Particles->vel[j*3+2]-velMirrorZi);

							// Refelected rij' = Rref_i * ri'j
							double v0ijm = (Rref_i[0]*v0imj + Rref_i[1]*v1imj + Rref_i[2]*v2imj);
							double v1ijm = (Rref_i[3]*v0imj + Rref_i[4]*v1imj + Rref_i[5]*v2imj);
							double v2ijm = (Rref_i[6]*v0imj + Rref_i[7]*v1imj + Rref_i[8]*v2imj);
							double invDstimj2 = 1.0/dstimj2;

							if(PSystem->gradientCorrection == false) {
								gradUxx += vimjx*v0ijm*wS*invDstimj2;
								gradUxy += vimjx*v1ijm*wS*invDstimj2;
								gradUxz += vimjx*v2ijm*wS*invDstimj2;
								
								gradUyx += vimjy*v0ijm*wS*invDstimj2;
								gradUyy += vimjy*v1ijm*wS*invDstimj2;
								gradUyz += vimjy*v2ijm*wS*invDstimj2;
								
								gradUzx += vimjz*v0ijm*wS*invDstimj2;
								gradUzy += vimjz*v1ijm*wS*invDstimj2;
								gradUzz += vimjz*v2ijm*wS*invDstimj2;
							}
							else {
								double v0ijmC = (v0ijm*MC[0] + v1ijm*MC[1] + v2ijm*MC[2]);
								double v1ijmC = (v0ijm*MC[3] + v1ijm*MC[4] + v2ijm*MC[5]);
								double v2ijmC = (v0ijm*MC[6] + v1ijm*MC[7] + v2ijm*MC[8]);
								gradUxx += vimjx*v0ijmC*wS*invDstimj2;
								gradUxy += vimjx*v1ijmC*wS*invDstimj2;
								gradUxz += vimjx*v2ijmC*wS*invDstimj2;
								
								gradUyx += vimjy*v0ijmC*wS*invDstimj2;
								gradUyy += vimjy*v1ijmC*wS*invDstimj2;
								gradUyz += vimjy*v2ijmC*wS*invDstimj2;
								
								gradUzx += vimjz*v0ijmC*wS*invDstimj2;
								gradUzy += vimjz*v1ijmC*wS*invDstimj2;
								gradUzz += vimjz*v2ijmC*wS*invDstimj2;
							}
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
				double vimix = -(velXi-velMirrorXi);
				double vimiy = -(velYi-velMirrorYi);
				double vimiz = -(velZi-velMirrorZi);

				// Refelected rii' = Rref_i * ri'i
				double v0iim = (Rref_i[0]*v0imi + Rref_i[1]*v1imi + Rref_i[2]*v2imi);
				double v1iim = (Rref_i[3]*v0imi + Rref_i[4]*v1imi + Rref_i[5]*v2imi);
				double v2iim = (Rref_i[6]*v0imi + Rref_i[7]*v1imi + Rref_i[8]*v2imi);
				double invDstimi2 = 1.0/dstimi2;

				if(PSystem->gradientCorrection == false) {
					gradUxx += vimix*v0iim*wS*invDstimi2;
					gradUxy += vimix*v1iim*wS*invDstimi2;
					gradUxz += vimix*v2iim*wS*invDstimi2;
					
					gradUyx += vimiy*v0iim*wS*invDstimi2;
					gradUyy += vimiy*v1iim*wS*invDstimi2;
					gradUyz += vimiy*v2iim*wS*invDstimi2;
					
					gradUzx += vimiz*v0iim*wS*invDstimi2;
					gradUzy += vimiz*v1iim*wS*invDstimi2;
					gradUzz += vimiz*v2iim*wS*invDstimi2;
				}
				else {
					double v0iimC = (v0iim*MC[0] + v1iim*MC[1] + v2iim*MC[2]);
					double v1iimC = (v0iim*MC[3] + v1iim*MC[4] + v2iim*MC[5]);
					double v2iimC = (v0iim*MC[6] + v1iim*MC[7] + v2iim*MC[8]);
					gradUxx += vimix*v0iimC*wS*invDstimi2;
					gradUxy += vimix*v1iimC*wS*invDstimi2;
					gradUxz += vimix*v2iimC*wS*invDstimi2;
					
					gradUyx += vimiy*v0iimC*wS*invDstimi2;
					gradUyy += vimiy*v1iimC*wS*invDstimi2;
					gradUyz += vimiy*v2iimC*wS*invDstimi2;
					
					gradUzx += vimiz*v0iimC*wS*invDstimi2;
					gradUzy += vimiz*v1iimC*wS*invDstimi2;
					gradUzz += vimiz*v2iimC*wS*invDstimi2;
				}
			}

			// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/noGrad)
			Particles->velGradientRow1[i*3  ] += -PSystem->coeffPressGrad*gradUxx;
			Particles->velGradientRow1[i*3+1] += -PSystem->coeffPressGrad*gradUxy;
			Particles->velGradientRow1[i*3+2] += -PSystem->coeffPressGrad*gradUxz;
			Particles->velGradientRow2[i*3  ] += -PSystem->coeffPressGrad*gradUyx;
			Particles->velGradientRow2[i*3+1] += -PSystem->coeffPressGrad*gradUyy;
			Particles->velGradientRow2[i*3+2] += -PSystem->coeffPressGrad*gradUyz;
			Particles->velGradientRow3[i*3  ] += -PSystem->coeffPressGrad*gradUzx;
			Particles->velGradientRow3[i*3+1] += -PSystem->coeffPressGrad*gradUzy;
			Particles->velGradientRow3[i*3+2] += -PSystem->coeffPressGrad*gradUzz;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Interpolate velocity according to a first-order Taylor expansion
void MpsShifting::interpolateVelocity(MpsParticleSystem *PSystem, MpsParticle *Particles) {
#pragma omp parallel for
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		// if(Particles->particleBC[i] == PSystem->inner) {
		// if(Particles->particleBC[i] == PSystem->inner && Particles->numNeigh[i] > 29) {
		if(Particles->particleType[i] == PSystem->fluid) {
		
		Particles->vel[i*3  ]-= Particles->velGradientRow1[i*3  ]*Particles->dPosShift[i*3  ] + 
								Particles->velGradientRow1[i*3+1]*Particles->dPosShift[i*3+1] +
								Particles->velGradientRow1[i*3+2]*Particles->dPosShift[i*3+2];
		
		Particles->vel[i*3+1]-= Particles->velGradientRow2[i*3  ]*Particles->dPosShift[i*3  ] + 
								Particles->velGradientRow2[i*3+1]*Particles->dPosShift[i*3+1] +
								Particles->velGradientRow2[i*3+2]*Particles->dPosShift[i*3+2];
		
		Particles->vel[i*3+2]-= Particles->velGradientRow3[i*3  ]*Particles->dPosShift[i*3  ] + 
								Particles->velGradientRow3[i*3+1]*Particles->dPosShift[i*3+1] +
								Particles->velGradientRow3[i*3+2]*Particles->dPosShift[i*3+2];

		// Particles->pos[i*3  ] += Particles->dPosShift[i*3  ];
		// Particles->pos[i*3+1] += Particles->dPosShift[i*3+1];
		// Particles->pos[i*3+2] += Particles->dPosShift[i*3+2];
		
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Normal vector on the fluid
void MpsShifting::calcNormalConcentration(MpsParticle *Particles) {
#pragma omp parallel for
	for(int ip=0; ip<Particles->numParticles; ip++) {
		int i = Particles->particleID[ip];
		double norm2GradCi = Particles->gradConcentration[i*3]*Particles->gradConcentration[i*3] + Particles->gradConcentration[i*3+1]*Particles->gradConcentration[i*3+1] + Particles->gradConcentration[i*3+2]*Particles->gradConcentration[i*3+2];
		
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