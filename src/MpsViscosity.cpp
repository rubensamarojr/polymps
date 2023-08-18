// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <experimental/filesystem> 	///< numeric_limits
#include <iostream>					///< cout
#include "MpsViscosity.h"

using namespace std;

// Constructor declaration
MpsViscosity::MpsViscosity()
{
}
// Destructor declaration
MpsViscosity::~MpsViscosity()
{
}

// Acceleration due Laplacian of velocity
void MpsViscosity::calcViscosity(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
//		if(Particles->particleType[i] == PSystem->fluid) {
		double meu_i = Particles->MEU[i];
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
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
				if(dstij2 < PSystem->reL2 && (dstij2 < dstimj2 || PSystem->wallType == boundaryWallType::PARTICLE)) {
					if(j != i) {
						double dst = sqrt(dstij2);
						double wL = Particles->weight(dst, PSystem->reL, PSystem->weightType);
						double neu_ij;
						if ((meu_i + Particles->MEU[j]) > PSystem->epsilonZero)
							neu_ij = 2.0 * meu_i * Particles->MEU[j] / (meu_i + Particles->MEU[j]);
						else
							neu_ij = 0.0;
						if(Particles->particleType[j] == PSystem->wall) neu_ij = 2.0 * meu_i; // Particles->MEU[j] -> oo
						//neu_ij = PSystem->KNM_VS2 * PSystem->DNS_FL2;
	//					if(Particles->PTYPE[i] == 1) neu_ij = neu_ij/PSystem->DNS_FL1;
	//					else neu_ij = neu_ij/PSystem->DNS_FL2;
						neu_ij = neu_ij/Particles->RHO[i];
						//if((NEUt[i] + NEUt[j]) > 0) neu_ij = neu_ij + (2.0 * NEUt[i] * Particles->RHO[j] * NEUt[j] * Particles->RHO[j] / (NEUt[i] * Particles->RHO[i] + NEUt[j] * Particles->RHO[j])) / Particles->RHO[i];
						// Original
	//					accX +=(Particles->vel[j*3  ]-velXi)*w;
	//					accY +=(Particles->vel[j*3+1]-velYi)*w;
	//					accZ +=(Particles->vel[j*3+2]-velZi)*w;
						// Modified
						accX +=(Particles->vel[j*3  ]-velXi)*wL*neu_ij;
						accY +=(Particles->vel[j*3+1]-velYi)*wL*neu_ij;
						accZ +=(Particles->vel[j*3+2]-velZi)*wL*neu_ij;
					}
				}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// Original
//		acc[i*3  ]=accX*PSystem->coeffViscosity + PSystem->gravityX;
//		acc[i*3+1]=accY*PSystem->coeffViscosity + PSystem->gravityY;
//		acc[i*3+2]=accZ*PSystem->coeffViscosity + PSystem->gravityZ;
		// Modified
		//if(PSystem->timeCurrent > 0.3) {
		// coeffViscMultiphase = 2.0*PSystem->dim/(PSystem->pndLargeZero*PSystem->lambdaZero);
		Particles->acc[i*3  ] = PSystem->coeffViscMultiphase*accX;
		Particles->acc[i*3+1] = PSystem->coeffViscMultiphase*accY;
		Particles->acc[i*3+2] = PSystem->coeffViscMultiphase*accZ;
		//}
				
		// Particles->accStar[i*3  ] = Particles->acc[i*3  ];
		// Particles->accStar[i*3+1] = Particles->acc[i*3+1];
		// Particles->accStar[i*3+2] = Particles->acc[i*3+2];
	}
	
#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}


// Calculation of the volume of fraction if phase II in the mixture
void MpsViscosity::calcVolumeFraction(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets)
{
	if(PSystem->Fraction_method == 1) {   //Linear distribution
#pragma omp parallel for schedule(dynamic,64)
		for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
			double sum1 = 0.0, sum2 = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			
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
					
					if(dstij2 < PSystem->reS2) {
					if(j != i && Particles->particleType[j] == PSystem->fluid) {
						sum1 += 1.0;
						if(Particles->PTYPE[j] >= 2) sum2 += 1.0;
						}}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			if(sum1 < PSystem->epsilonZero)
				Particles->Cv[i] = 0.0;
			else 
				Particles->Cv[i] = sum2/sum1;
		}}
	}
	else if(PSystem->Fraction_method == 2) {   //Non linear :  Smoothed using the weight funtion
#pragma omp parallel for schedule(dynamic,64)
		for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
			double sum1 = 0.0, sum2 = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			
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

					if(dstij2 < PSystem->reS2) {
					if(j != i && Particles->particleType[j] == PSystem->fluid) {
						double dst = sqrt(dstij2);
						double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
						sum1 += wS;
						if(Particles->PTYPE[j] >= 2) sum2 += wS;
						}}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			if(sum1 < PSystem->epsilonZero)
				Particles->Cv[i] = 0.0;
			else 
				Particles->Cv[i] = sum2/sum1;
		}}
	}
#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// 2018 - Meshfree particle numerical modelling of sub-aerial and submerged landslides
// Viscosity interaction values for "real" PSystem->fluid particles
void MpsViscosity::calcViscosityInteractionVal(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {

	double gravityMod = sqrt(PSystem->gravityX*PSystem->gravityX + PSystem->gravityY*PSystem->gravityY + PSystem->gravityZ*PSystem->gravityZ);

	//double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;//,grain_VF, *p_smooth;
	double phi = 0.0, phi2 = 0.0, meu_0, normal_stress;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
	// Xmin = 0.0 - PSystem->partDist*3.0; Xmax = 1.65 + PSystem->partDist*3.0;
	// Ymin = 0.0 - PSystem->partDist*3.0; Ymax = 0.15 + PSystem->partDist*3.0;
	// Zmin = 0.0 - PSystem->partDist*3.0; //Zmax = 0.7 + PSystem->partDist*30.0;
	// damErosion3D
	// Xmin = 0.0 - PSystem->partDist*3.0; Xmax = 2.00 + PSystem->partDist*3.0;
	// Ymin = 0.0 - PSystem->partDist*3.0; Ymax = 0.10 + PSystem->partDist*3.0;
	// damErosion2D
	// Xmin = 0.0 - PSystem->partDist*3.0; Xmax = 2.00 + PSystem->partDist*3.0;
	// Ymin = 0.0; Ymax = 0.00;
	// S1 2D
	// Xmin = 0.0 - PSystem->partDist*3.0; Xmax = 0.30 + PSystem->partDist*3.0;
	// Ymin = 0.0; Ymax = 0.00;
	// Subaquatic 0.016
	// Xmin = 0.0016 - PSystem->partDist*3.0; Xmax = 0.9792 + PSystem->partDist*3.0;
	// Ymin = 0.00; Ymax = 0.00;
	// Changed !!!

	Xmin = PSystem->domainMinX;
	Xmax = PSystem->domainMaxX;
	Ymin = PSystem->domainMinY;
	Ymax = Ymin;	Zmin = 0.0;
	if ((int)PSystem->dim == 3) {
		Ymax = PSystem->domainMaxY;
		Zmin = PSystem->domainMinZ;
	}

	// Search free-surface particles for each interval of aa = 2 particles in wall
	double aaL0 = 2.0*PSystem->partDist;
	int kx_max = int( (Xmax - Xmin)/aaL0 ) + 1;
	int ky_max = int( (Ymax - Ymin)/aaL0 ) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
/*
	S11 = new double[Particles->numParticles + 1];
	S22 = new double[Particles->numParticles + 1];
	S33 = new double[Particles->numParticles + 1];
	S12 = new double[Particles->numParticles + 1];
	S13 = new double[Particles->numParticles + 1];
	S23 = new double[Particles->numParticles + 1];
*/
	//p_smooth = new double[Particles->numParticles + 1];
	BL = new double*[kx_max + 1];  // bed level
	WL = new double*[kx_max + 1];  // water level
	PS = new double*[kx_max + 1];  // pressure sediment

#pragma omp parallel for
	for(int m = 1; m <= kx_max; m++) {
		BL[m] = new double[ky_max + 1];
		WL[m] = new double[ky_max + 1];
		PS[m] = new double[ky_max + 1];
	}

	// Determining the bed level
	if((int)PSystem->dim == 2) {
#pragma omp parallel for schedule(dynamic,64)
		for(int kx = 1; kx <= kx_max; kx++) {
			for(int ky = 1; ky <= ky_max; ky++) {
				BL[kx][ky] = Ymin;
				WL[kx][ky] = Ymin;
				PS[kx][ky] = 0.0;
			}
		}
	}
	else {
#pragma omp parallel for schedule(dynamic,64)
		for(int kx = 1; kx <= kx_max; kx++) {
			for(int ky = 1; ky <= ky_max; ky++) {
				BL[kx][ky] = Zmin;
				WL[kx][ky] = Zmin;
				PS[kx][ky] = 0.0;
			}
		}
	}

	if((int)PSystem->dim == 2) {
#pragma omp parallel for
		for(int i=0; i<Particles->numParticles; i++) {
			if(Particles->particleType[i] == PSystem->fluid) {
				double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];
			
				int kx = int( (posXi - Xmin)/aaL0 ) + 1;
				int ky = 1;
				//if(posYi > BL[kx][ky] && Particles->Cv[i] > 0.5) { BL[kx][ky] = posYi; PS[kx][ky] = pnew[i]; }
				if(posYi > BL[kx][ky] && Particles->Cv[i] > 0.5) { BL[kx][ky] = posYi; PS[kx][ky] = Particles->press[i]; }
				if(posYi > WL[kx][ky] && Particles->PTYPE[i] == 1) { WL[kx][ky] = posYi; }
			}
		}
	}
	else {
#pragma omp parallel for
		for(int i=0; i<Particles->numParticles; i++) {
			if(Particles->particleType[i] == PSystem->fluid) {
				double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			
				int kx = int( (posXi - Xmin)/aaL0 ) + 1;
				int ky = int( (posYi - Ymin)/aaL0 ) + 1;
				//if(posZi > BL[kx][ky] && Particles->Cv[i] > 0.5) { BL[kx][ky] = posZi; PS[kx][ky] = pnew[i]; }
				if(posZi > BL[kx][ky] && Particles->Cv[i] > 0.5) { BL[kx][ky] = posZi; PS[kx][ky] = Particles->press[i]; }
				if(posZi > WL[kx][ky] && Particles->PTYPE[i] == 1) { WL[kx][ky] = posZi; }
			}
		}
	}

	// Strain rate calculation
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
	if(Particles->particleType[i] == PSystem->fluid) {
		double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0, sum8 = 0.0, sum9 = 0.0, sum10 = 0.0;
		double sumWs = 0.0;
//		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double MC[9];
		MC[0] = Particles->correcMatrixRow1[i*3];	MC[1] = Particles->correcMatrixRow1[i*3+1];	MC[2] = Particles->correcMatrixRow1[i*3+2];
		MC[3] = Particles->correcMatrixRow2[i*3];	MC[4] = Particles->correcMatrixRow2[i*3+1];	MC[5] = Particles->correcMatrixRow2[i*3+2];
		MC[6] = Particles->correcMatrixRow3[i*3];	MC[7] = Particles->correcMatrixRow3[i*3+1];	MC[8] = Particles->correcMatrixRow3[i*3+2];
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
					double vec_ijx = Particles->vel[j*3  ] - velXi;	
					double vec_ijy = Particles->vel[j*3+1] - velYi;	
					double vec_ijz = Particles->vel[j*3+2] - velZi;
					double invDstij2 = 1.0/dstij2;

					if(PSystem->gradientCorrection == false) {
						sum1 += vec_ijx*v0ij*wS*invDstij2;
						sum2 += vec_ijx*v1ij*wS*invDstij2;
						sum3 += vec_ijx*v2ij*wS*invDstij2;
						
						sum4 += vec_ijy*v0ij*wS*invDstij2;
						sum5 += vec_ijy*v1ij*wS*invDstij2;
						sum6 += vec_ijy*v2ij*wS*invDstij2;
						
						sum7 += vec_ijz*v0ij*wS*invDstij2;
						sum8 += vec_ijz*v1ij*wS*invDstij2;
						sum9 += vec_ijz*v2ij*wS*invDstij2;
					}
					else {
						double v0ijC = (v0ij*MC[0] + v1ij*MC[1] + v2ij*MC[2]);
						double v1ijC = (v0ij*MC[3] + v1ij*MC[4] + v2ij*MC[5]);
						double v2ijC = (v0ij*MC[6] + v1ij*MC[7] + v2ij*MC[8]);
						sum1 += vec_ijx*v0ijC*wS*invDstij2;
						sum2 += vec_ijx*v1ijC*wS*invDstij2;
						sum3 += vec_ijx*v2ijC*wS*invDstij2;
						
						sum4 += vec_ijy*v0ijC*wS*invDstij2;
						sum5 += vec_ijy*v1ijC*wS*invDstij2;
						sum6 += vec_ijy*v2ijC*wS*invDstij2;
						
						sum7 += vec_ijz*v0ijC*wS*invDstij2;
						sum8 += vec_ijz*v1ijC*wS*invDstij2;
						sum9 += vec_ijz*v2ijC*wS*invDstij2;
					}
					
					sum10 += Particles->press[j]*wS;
					sumWs += wS;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/PSystem->pndGradientZero)
		Uxx = -PSystem->coeffPressGrad*sum1; Uxy = -PSystem->coeffPressGrad*sum2; Uxz = -PSystem->coeffPressGrad*sum3;
		Uyx = -PSystem->coeffPressGrad*sum4; Uyy = -PSystem->coeffPressGrad*sum5; Uyz = -PSystem->coeffPressGrad*sum6;
		Uzx = -PSystem->coeffPressGrad*sum7; Uzy = -PSystem->coeffPressGrad*sum8; Uzz = -PSystem->coeffPressGrad*sum9;

		/*if(PSystem->gradientCorrection == true) {
			double Uaux[9];
			Uaux[0] = Uxx*Particles->correcMatrixRow1[i*3] + Uyx*Particles->correcMatrixRow1[i*3+1] + Uzx*Particles->correcMatrixRow1[i*3+2];
			Uaux[1] = Uxy*Particles->correcMatrixRow1[i*3] + Uyy*Particles->correcMatrixRow1[i*3+1] + Uzy*Particles->correcMatrixRow1[i*3+2];
			Uaux[2] = Uxz*Particles->correcMatrixRow1[i*3] + Uyz*Particles->correcMatrixRow1[i*3+1] + Uzz*Particles->correcMatrixRow1[i*3+2];
			Uaux[3] = Uxx*Particles->correcMatrixRow2[i*3] + Uyx*Particles->correcMatrixRow2[i*3+1] + Uzx*Particles->correcMatrixRow2[i*3+2];
			Uaux[4] = Uxy*Particles->correcMatrixRow2[i*3] + Uyy*Particles->correcMatrixRow2[i*3+1] + Uzy*Particles->correcMatrixRow2[i*3+2];
			Uaux[5] = Uxz*Particles->correcMatrixRow2[i*3] + Uyz*Particles->correcMatrixRow2[i*3+1] + Uzz*Particles->correcMatrixRow2[i*3+2];
			Uaux[6] = Uxx*Particles->correcMatrixRow3[i*3] + Uyx*Particles->correcMatrixRow3[i*3+1] + Uzx*Particles->correcMatrixRow3[i*3+2];
			Uaux[7] = Uxy*Particles->correcMatrixRow3[i*3] + Uyy*Particles->correcMatrixRow3[i*3+1] + Uzy*Particles->correcMatrixRow3[i*3+2];
			Uaux[8] = Uxz*Particles->correcMatrixRow3[i*3] + Uyz*Particles->correcMatrixRow3[i*3+1] + Uzz*Particles->correcMatrixRow3[i*3+2];

			Uxx = Uaux[0];	Uxy = Uaux[1];	Uxz = Uaux[2];
			Uyx = Uaux[3];	Uyy = Uaux[4];	Uyz = Uaux[5];
			Uzx = Uaux[6];	Uzy = Uaux[7];	Uzz = Uaux[8];
		}
		*/
		if(sumWs > PSystem->epsilonZero) {
			Particles->p_smooth[i] = sum10/sumWs;
		}
		else {
			Particles->p_smooth[i] = Particles->press[i];
		}
		if(Particles->p_smooth[i] < PSystem->epsilonZero) Particles->p_smooth[i] = 0.0;

		Particles->S11[i] = 0.5*(Uxx + Uxx);
		Particles->S12[i] = 0.5*(Uxy + Uyx);
		Particles->S13[i] = 0.5*(Uxz + Uzx);
		Particles->S22[i] = 0.5*(Uyy + Uyy);
		Particles->S23[i] = 0.5*(Uyz + Uzy);
		Particles->S33[i] = 0.5*(Uzz + Uzz);

		// Square of the second invariant of the strain-rate tensor = 0.5*tr(SS^2)
		//Particles->II[i] = 0.5*Uxx*Uxx + 0.5*Uyy*Uyy + 0.25*(Uxy + Uyx)*(Uxy + Uyx);
		//Particles->II[i] = 0.5*(Particles->S11[i]*Particles->S11[i] + Particles->S12[i]*Particles->S12[i] + Particles->S13[i]*Particles->S13[i] + Particles->S12[i]*Particles->S12[i] + Particles->S22[i]*Particles->S22[i] + Particles->S23[i]*Particles->S23[i] + Particles->S13[i]*Particles->S13[i] + Particles->S23[i]*Particles->S23[i] + Particles->S33[i]*Particles->S33[i]);
		Particles->II[i] = 0.5*(Particles->S11[i]*Particles->S11[i] + 2.0*Particles->S12[i]*Particles->S12[i] + 2.0*Particles->S13[i]*Particles->S13[i] + Particles->S22[i]*Particles->S22[i] + 2.0*Particles->S23[i]*Particles->S23[i] + Particles->S33[i]*Particles->S33[i]);
//		Particles->II[i] = - (Particles->S11[i]*Particles->S22[i] + Particles->S22[i]*Particles->S33[i] + Particles->S11[i]*Particles->S33[i] - Particles->S12[i]*Particles->S12[i] - Particles->S13[i]*Particles->S13[i] - Particles->S23[i]*Particles->S23[i]);
//		Particles->II[i] = sqrt(Particles->II[i]*Particles->II[i]);
		if(Particles->II[i] < PSystem->epsilonZero || Particles->II[i]*0 != 0) Particles->II[i] = 0.0;
		//II=fabs(Particles->S11[i]*Particles->S22[i]-Particles->S12[i]*Particles->S12[i]);
		//std::cout << " II: " << Particles->II[i] << std::endl;
	}}

	// Newtonian viscosity
	if(PSystem->fluidType == viscType::NEWTONIAN)
	{
#pragma omp parallel for
		for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
			if(Particles->PTYPE[i] <= 1)Particles->MEU[i] = PSystem->KNM_VS1 * PSystem->DNS_FL1;
			if(Particles->PTYPE[i] != 1)Particles->MEU[i] = PSystem->KNM_VS2 * PSystem->DNS_FL2;
		}}

//		if(TURB > 0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2.0*sqrt(Particles->II[i]);

//			if(NEUt[i] * 0 != 0)  NEUt[i] = 0.0;
//			if(NEUt[i] > 1.0)     NEUt[i] = 1.0;
//		}
	}

	double mi_max = 0.0;
	// Granular PSystem->fluid
	if(PSystem->fluidType == viscType::NON_NEWTONIAN)
	{
	// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
//#pragma omp parallel for
		for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {

//			double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double vel2i = velXi*velXi + velYi*velYi + velZi*velZi;

			if(Particles->PTYPE[i] == 1) { // Newtonian PSystem->fluid
				// Mojtaba
				//Particles->MEU[i] = PSystem->KNM_VS1*PSystem->DNS_FL1*(1 + 2.5*Particles->Cv[i]);
				//////////////
				// Original //
				Particles->MEU[i] = PSystem->KNM_VS1*PSystem->DNS_FL1;
			}
			else if(Particles->PTYPE[i] == 2) { // Non-Newtonian mixture
				//////////////
				// Original //
				// phi: internal friction angle
				// phi2: maximum friction angle
				phi = (Particles->Cv[i] - 0.25)*PSystem->PHI_1/(1.0 - 0.25);
				phi2 = (Particles->Cv[i] - 0.25)*PSystem->PHI_2/(1.0 - 0.25);
				if(Particles->Cv[i] <= 0.25) { phi = PSystem->epsilonZero; phi2 = PSystem->epsilonZero; } // phi close to zero
				// if(Particles->PTYPE[i] <= 0) phi = PSystem->PHI_BED; // ghost

				// normal stress calculation (mechanical pressure)
				Particles->p_rheo_new[i] = Particles->p_smooth[i];

				int kx = int( (posXi - Xmin)/aaL0 ) + 1;
				//int ky = int( (posYi - Ymin)/aaL0 ) + 1;
				int ky;
				if((int)PSystem->dim == 2) {
					ky = 1;
					// Effective pressure = total pressure (from EOS) - hydrostatic pressure
					//normal_stress = (BL[kx][ky] - posYi + PSystem->partDist*0.5)*(PSystem->DNS_FL2)*gravityMod;	// normal_stress= Gama.H
					// SPH Simulation of Sediment Flushing Induced by a Rapid Water Flow
					normal_stress = (BL[kx][ky] - posYi + PSystem->partDist*0.5)*(PSystem->DNS_FL2 - PSystem->DNS_FL1)*gravityMod - vel2i*(PSystem->DNS_FL2 - PSystem->DNS_FL1)*0.5;	// normal_stress= Gama.H, Eq. (8)

//					if(Particles->p_smooth[i] < (WL[kx][ky] - posYi)*PSystem->DNS_FL1*gravityMod) Particles->p_smooth[i] = (WL[kx][ky] - posYi)*PSystem->DNS_FL1*gravityMod;
					//if(PSystem->timeCurrent <= 1.0) normal_stress = (1.0 - PSystem->timeCurrent)*(Particles->p_smooth[i] - (WL[kx][ky] - posYi)*PSystem->DNS_FL1*gravityMod) + PSystem->timeCurrent*normal_stress;

					if(WL[kx][ky] < BL[kx][ky]) {
						normal_stress = Particles->p_smooth[i];		// Free-fall (dry granular material)
					}
					else {
						//normal_stress = Particles->p_smooth[i] - (WL[kx][ky] - posYi)*PSystem->DNS_FL1*gravityMod;
						// A multiphase meshfree particle method for continuum-based modeling of dry and submerged granular flows
						normal_stress = (Particles->p_smooth[i] - PS[kx][ky])*(PSystem->DNS_FL2 - PSystem->DNS_FL1)*Particles->VF[i]/Particles->RHO[i];	// Grain inertia (submerged) Eq. (20)

						normal_stress = Particles->p_smooth[i] - PS[kx][ky];
//						normal_stress = Particles->p_smooth[i] - (WL[kx][ky] - posYi)*(PSystem->DNS_FL1)*gravityMod;	// Grain inertia (submerged)
					}
				}
				else {
					ky = int( (posYi - Ymin)/aaL0 ) + 1;
					// Effective pressure = total pressure (from EOS) - hydrostatic pressure
					//normal_stress = (BL[kx][ky] - posZi + PSystem->partDist*0.5)*(PSystem->DNS_FL2)*gravityMod;	// normal_stress= Gama.H
					normal_stress = (BL[kx][ky] - posZi + PSystem->partDist*0.5)*(PSystem->DNS_FL2 - PSystem->DNS_FL1)*gravityMod - vel2i*(PSystem->DNS_FL2 - PSystem->DNS_FL1)*0.5;	// normal_stress= Gama.H

					if(Particles->p_smooth[i] - (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod < PSystem->epsilonZero) Particles->p_smooth[i] = (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod;
					if(PSystem->timeCurrent <= 1.0) normal_stress = (1.0 - PSystem->timeCurrent)*(Particles->p_smooth[i] - (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod) + PSystem->timeCurrent*normal_stress;

//					normal_stress = Particles->p_smooth[i] - (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod;
				}

				
//				normal_stress = Particles->p_smooth[i];
				//normal_stress=normal_stress*0.61*1500/PSystem->DNS_FL2;
				if(normal_stress < 1.0 || Particles->Cv[i] < 0.5) normal_stress = 1.0;

				Particles->p_rheo_new[i] = normal_stress;

				double modE = sqrt(Particles->II[i]);

				// Yield stress calculation (Text below Eq. (8))
				if(WL[kx][ky] < BL[kx][ky]) {
					Particles->Inertia[i] = modE*PSystem->DG/sqrt(normal_stress/PSystem->DNS_SDT);		// Free-fall (dry granular material)
				}
				else {
					Particles->Inertia[i] = modE*PSystem->DG/sqrt(normal_stress/(PSystem->DNS_FL1*PSystem->Cd));	// Grain inertia (submerged)
				}
				//Particles->Inertia[i] = modE*(PSystem->KNM_VS1*PSystem->DNS_FL1)/normal_stress ;	// Viscous regime

//				Particles->Inertia[i] = 1.0;
				// PSystem->VF_max PSystem->VF_min
				Particles->VF[i] = PSystem->VF_max - (PSystem->VF_max - PSystem->VF_min)*Particles->Inertia[i];
				if(Particles->VF[i] < PSystem->VF_min) Particles->VF[i] = PSystem->VF_min;
				Particles->RHO[i] = PSystem->DNS_SDT*Particles->VF[i] + (1.0 - Particles->VF[i])*PSystem->DNS_FL1;
				phi *= (Particles->VF[i]/PSystem->VF_max);

				double yield_stress;

				// Drucker-Prager model
				double alpha_1 = 2.0*sqrt(3.0)*sin(phi)/(3.0 - sin(phi));
				double beta_1 = 2.0*sqrt(3.0)*cos(phi)/(3.0 - sin(phi));
				// Mohr-Coulomb model
				//double alpha_1 = sin(phi);
				//double beta_1 = cos(phi);
				//double alpha_1 = tan(phi);
				//double beta_1 = 1.0;

				yield_stress = PSystem->cohes*beta_1 + normal_stress*alpha_1;

				//yield_stress = normal_stress*tan(phi);

				if(yield_stress < PSystem->epsilonZero) yield_stress = 0.0;

				double visc_max = (PSystem->MEU0 + yield_stress*PSystem->mm*0.5); // Below Eq. (5)

				// Pre-failure portion
				//if(modE > PSystem->epsilonZero) {
				//	Particles->MEU_Y[i] = yield_stress*(1.0 - exp(-PSystem->mm*modE))/(2.0*modE); // Eq. (5)
				//}
				//else {
				//	Particles->MEU_Y[i] = visc_max;
				//}


				Particles->MEU_Y[i] = yield_stress/(2.0*sqrt(modE*modE + PSystem->epsilonZero));

				// H-B rheology

				//meu_0 = PSystem->MEU0;

//				phi = PSystem->PHI_1; phi2 = PSystem->PHI_2;
//				if(Particles->Cv[i] <= 0.25) { phi = PSystem->epsilonZero; phi2 = PSystem->epsilonZero; } // phi close to zero

				// Non-linear Meu(I) rheology
				if(WL[kx][ky] < BL[kx][ky]) {
//					meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*PSystem->DG/(PSystem->I0*sqrt(normal_stress/PSystem->DNS_FL2) + modE*PSystem->DG);			//free fall
					meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*PSystem->DG/(PSystem->I0*sqrt(normal_stress/PSystem->DNS_FL2) + modE*PSystem->DG)*modE/(sqrt(modE*modE + PSystem->epsilonZero));
					//meu_0 = (tan(phi) + (tan(phi2) - tan(phi))*modE*PSystem->DG/(PSystem->I0*sqrt(normal_stress/PSystem->DNS_FL2) + modE*PSystem->DG))*normal_stress/(sqrt(Particles->II[i] + PSystem->epsilonZero));
				}
				else {
//					meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*PSystem->DG/(PSystem->I0*sqrt(normal_stress/(PSystem->DNS_FL1*PSystem->Cd)) + modE*PSystem->DG);	//grain inertia
					meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*PSystem->DG/(PSystem->I0*sqrt(normal_stress/(PSystem->DNS_FL1*PSystem->Cd)) + modE*PSystem->DG)*modE/(sqrt(modE*modE + PSystem->epsilonZero));

					//meu_0 = (tan(phi) + (tan(phi2) - tan(phi))*modE*PSystem->DG/(PSystem->I0*sqrt(normal_stress/PSystem->DNS_FL1*PSystem->Cd) + modE*PSystem->DG))*normal_stress/(sqrt(Particles->II[i] + PSystem->epsilonZero));
				}

			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(PSystem->KNM_VS1*PSystem->DNS_FL1)/(PSystem->I0*normal_stress + modE*(PSystem->KNM_VS1*PSystem->DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*PSystem->DG*sqrt(normal_stress*PSystem->DNS_FL2)/PSystem->I0;		//free fall
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*PSystem->DG*sqrt(normal_stress*PSystem->DNS_FL1*PSystem->Cd)/PSystem->I0;	//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*(PSystem->KNM_VS1*PSystem->DNS_FL1)/PSystem->I0;					//viscous

				if(modE <= PSystem->epsilonZero || (meu_0*0) != 0) meu_0 = PSystem->MEU0;

				visc_max = (meu_0 + yield_stress*PSystem->mm*0.5); // Below Eq. (5)

				//if(isnan(Particles->II[i]) || isinf(Particles->II[i])) {
					//std::cout << " viscmax: " << Particles->II[i] << std::endl;
				//	assert(visc_max >= 0.0 || visc_max <= 0.0);
				//}
				
				// Herschel bulkley papanastasiou
//				Particles->MEU[i] = Particles->MEU_Y[i] + PSystem->MEU0*pow(4.0*Particles->II[i], (PSystem->N - 1.0)*0.5);

				// MEU_Y rheological model
				Particles->MEU[i] = Particles->MEU_Y[i] + meu_0;
				
				//if(Particles->II[i] <= PSystem->epsilonZero || Particles->MEU[i] > visc_max) {
				if(Particles->II[i] <= PSystem->epsilonZero) {
					//std::cout << " MEU>viscmax: " << yield_stress*PSystem->mm*0.5 << " PSystem->MEU0: " << meu_0 << " II: " << Particles->II[i] << std::endl;
					Particles->MEU[i] = visc_max;
				}
				// if(Particles->PTYPE[i] <= 0) Particles->MEU[i] = Particles->MEU[i]*Particles->Cv[i] + PSystem->DNS_FL1*PSystem->KNM_VS1*(1.0 - Particles->Cv[i]); // ghost

				//if(Particles->MEU[i]/Particles->RHO[i] > maxVIS) maxVIS = Particles->MEU[i]/Particles->RHO[i];
				if(Particles->MEU[i] > mi_max) mi_max = Particles->MEU[i];
			}
			
			if(Particles->PTYPE[i] >= 2) {
				if(Particles->Cv[i] > 0.5) Particles->RHO[i] = PSystem->DNS_FL2;
				else Particles->RHO[i] = Particles->Cv[i]*PSystem->DNS_FL2 + (1.0 - Particles->Cv[i])*PSystem->DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if(stress_cal_method == 2)
//		{
//			for(i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0, sum8 = 0.0, sum9 = 0.0, sum10 = 0.0;
//				for(l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if(i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2.0 * Particles->MEU[i] * Particles->MEU[j] / (Particles->MEU[i] + Particles->MEU[j]);
//						if((NEUt[i] + NEUt[j])>0) meuij = meuij + 2.0 * NEUt[i] * Particles->RHO[i] * NEUt[j] * Particles->RHO[j] / (NEUt[i] * Particles->RHO[i] + NEUt[j] * Particles->RHO[j]);

//						sum1 = sum1 + meuij * (x_vel[j] - x_vel[i])*DX(i, j)*w / d / d;
//						sum2 = sum2 + meuij * (x_vel[j] - x_vel[i])*DY(i, j)*w / d / d;
//						sum3 = sum3 + meuij * (x_vel[j] - x_vel[i])*DZ(i, j)*w / d / d;

//						sum4 = sum4 + meuij * (y_vel[j] - y_vel[i])*DX(i, j)*w / d / d;
//						sum5 = sum5 + meuij * (y_vel[j] - y_vel[i])*DY(i, j)*w / d / d;
//						sum6 = sum6 + meuij * (y_vel[j] - y_vel[i])*DZ(i, j)*w / d / d;

//						sum7 = sum7 + meuij * (z_vel[j] - z_vel[i])*DX(i, j)*w / d / d;
//						sum8 = sum8 + meuij * (z_vel[j] - z_vel[i])*DY(i, j)*w / d / d;
//						sum9 = sum9 + meuij * (z_vel[j] - z_vel[i])*DZ(i, j)*w / d / d;
//					}
//				}

//				Tau_xx[i] = (PSystem->dim / n0) * 2.0 * sum1;
//				Tau_yy[i] = (PSystem->dim / n0) * 2.0 * sum5;
//				Tau_zz[i] = (PSystem->dim / n0) * 2.0 * sum9;

//				Tau_xy[i] = (PSystem->dim / n0)*(sum2 + sum4);
//				Tau_xz[i] = (PSystem->dim / n0)*(sum3 + sum7);
//				Tau_yz[i] = (PSystem->dim / n0)*(sum6 + sum8);
//			}
//		}

	} // if(PSystem->fluidType == viscType::NON_NEWTONIAN)

	//---------------------------------------------------------------

//	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS; //delete[]p_smooth;
//	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL; //p_smooth = NULL;

	PSystem->CFLvisc = PSystem->timeStep*mi_max/((PSystem->DNS_SDT*PSystem->VF_max + (1.0 - PSystem->VF_max)*PSystem->DNS_FL1)*PSystem->partDist*PSystem->partDist);

	delete[]BL; delete[]WL; delete[]PS;
	BL = NULL; WL = NULL; PS = NULL;

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Free-slip condition. Viscosity interaction values
void MpsViscosity::calcWallSlipViscosityInteractionVal(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {

	double gravityMod = sqrt(PSystem->gravityX*PSystem->gravityX + PSystem->gravityY*PSystem->gravityY + PSystem->gravityZ*PSystem->gravityZ);

	//double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;//, grain_VF, *p_smooth;
	double phi = 0.0, phi2 = 0.0, meu_0, normal_stress;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
	// Xmin = 0.0 - PSystem->partDist*3.0; Xmax = 1.65 + PSystem->partDist*3.0;
	// Ymin = 0.0 - PSystem->partDist*3.0; Ymax = 0.15 + PSystem->partDist*3.0;
	// Zmin = 0.0 - PSystem->partDist*3.0; //Zmax = 0.7 + PSystem->partDist*30.0;
	// damErosion3D
	// Xmin = 0.0 - PSystem->partDist*3.0; Xmax = 2.00 + PSystem->partDist*3.0;
	// Ymin = 0.0 - PSystem->partDist*3.0; Ymax = 0.10 + PSystem->partDist*3.0;
	// Changed !!!

	Xmin = PSystem->domainMinX;
	Xmax = PSystem->domainMaxX;
	Ymin = PSystem->domainMinY;
	Ymax = Ymin;	Zmin = 0.0;
	if ((int)PSystem->dim == 3) {
		Ymax = PSystem->domainMaxY;
		Zmin = PSystem->domainMinZ;
	}

	// Search free-surface particles for each interval of aa = 2 particles in wall
	double aaL0 = 2.0*PSystem->partDist;
	int kx_max = int( (Xmax - Xmin)/aaL0 ) + 1;
	int ky_max = int( (Ymax - Ymin)/aaL0 ) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
	double aUxx, aUxy, aUxz, aUyx, aUyy, aUyz, aUzx, aUzy, aUzz;
/*
	S11 = new double[Particles->numParticles + 1];
	S22 = new double[Particles->numParticles + 1];
	S33 = new double[Particles->numParticles + 1];
	S12 = new double[Particles->numParticles + 1];
	S13 = new double[Particles->numParticles + 1];
	S23 = new double[Particles->numParticles + 1];
*/
	//p_smooth = new double[Particles->numParticles + 1];
	BL = new double*[kx_max + 1];  // bed level
	WL = new double*[kx_max + 1];  // water level
	PS = new double*[kx_max + 1];  // pressure sediment

#pragma omp parallel for
	for(int m = 1; m <= kx_max; m++)
	{
		BL[m] = new double[ky_max + 1];
		WL[m] = new double[ky_max + 1];
		PS[m] = new double[ky_max + 1];
	}

	// Determining the bed level
#pragma omp parallel for schedule(dynamic,64)
	for(int kx = 1; kx <= kx_max; kx++)
	{
		for(int ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			
			int kx = int( (posXi - Xmin)/aaL0 ) + 1;
			int ky = int( (posYi - Ymin)/aaL0 ) + 1;

			//if(posZi>BL[kx][ky] && Particles->Cv[i]>0.5) { BL[kx][ky] = posZi; PS[kx][ky] = pnew[i]; }
			if(posZi>BL[kx][ky] && Particles->Cv[i]>0.5) { BL[kx][ky] = posZi; PS[kx][ky] = Particles->press[i]; }
			if(posZi>WL[kx][ky] && Particles->PTYPE[i] == 1) { WL[kx][ky] = posZi; }
		}
	}

	// Strain rate calculation
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
	//if(Particles->particleType[i] == PSystem->fluid) {
	if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
		double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0, sum8 = 0.0, sum9 = 0.0, sum10 = 0.0;

//		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2

		// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
		double Rref_i[9], normaliw[3], normalMod2;
		// normal PSystem->fluid-wall particle = 0.5*(normal PSystem->fluid-mirror particle)
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

		//  Transformation matrix R_i = I - 2.0*normal_iwall*normal_iwall
		Rref_i[0] = 1.0 - 2.0*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2.0*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2.0*normaliw[0]*normaliw[2];
		Rref_i[3] = 0.0 - 2.0*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2.0*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2.0*normaliw[1]*normaliw[2];
		Rref_i[6] = 0.0 - 2.0*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2.0*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2.0*normaliw[2]*normaliw[2];

		// Mirror particle velocity vi' = Ri * vi
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

					double vec_mijx = Particles->vel[j*3  ]-velMirrorXi;	
					double vec_mijy = Particles->vel[j*3+1]-velMirrorYi;	
					double vec_mijz = Particles->vel[j*3+2]-velMirrorZi;

					sum1 += vec_mijx*v0imj*wS/dstimj2;
					sum2 += vec_mijx*v1imj*wS/dstimj2;
					sum3 += vec_mijx*v2imj*wS/dstimj2;
					
					sum4 += vec_mijy*v0imj*wS/dstimj2;
					sum5 += vec_mijy*v1imj*wS/dstimj2;
					sum6 += vec_mijy*v2imj*wS/dstimj2;
					
					sum7 += vec_mijz*v0imj*wS/dstimj2;
					sum8 += vec_mijz*v1imj*wS/dstimj2;
					sum9 += vec_mijz*v2imj*wS/dstimj2;
					
					sum10 += Particles->press[j]*wS;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Rref_i * gradU
		// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/PSystem->pndGradientZero)
		aUxx = -PSystem->coeffPressGrad*(Rref_i[0]*sum1 + Rref_i[1]*sum4 + Rref_i[2]*sum7);
		aUxy = -PSystem->coeffPressGrad*(Rref_i[0]*sum2 + Rref_i[1]*sum5 + Rref_i[2]*sum8);
		aUxz = -PSystem->coeffPressGrad*(Rref_i[0]*sum3 + Rref_i[1]*sum6 + Rref_i[2]*sum9);

		aUyx = -PSystem->coeffPressGrad*(Rref_i[3]*sum1 + Rref_i[4]*sum4 + Rref_i[5]*sum7);
		aUyy = -PSystem->coeffPressGrad*(Rref_i[3]*sum2 + Rref_i[4]*sum5 + Rref_i[5]*sum8);
		aUyz = -PSystem->coeffPressGrad*(Rref_i[3]*sum3 + Rref_i[4]*sum6 + Rref_i[5]*sum9);

		aUzx = -PSystem->coeffPressGrad*(Rref_i[6]*sum1 + Rref_i[7]*sum4 + Rref_i[8]*sum7);
		aUzy = -PSystem->coeffPressGrad*(Rref_i[6]*sum2 + Rref_i[7]*sum5 + Rref_i[8]*sum8);
		aUzz = -PSystem->coeffPressGrad*(Rref_i[6]*sum3 + Rref_i[7]*sum6 + Rref_i[8]*sum9);

		// Rref_i * gradU * Rref_i
		Uxx = aUxx*Rref_i[0] + aUxy*Rref_i[3] + aUxz*Rref_i[6];
		Uxy = aUxx*Rref_i[1] + aUxy*Rref_i[4] + aUxz*Rref_i[7];
		Uxz = aUxx*Rref_i[2] + aUxy*Rref_i[5] + aUxz*Rref_i[8];

		Uyx = aUyx*Rref_i[0] + aUyy*Rref_i[3] + aUyz*Rref_i[6];
		Uyy = aUyx*Rref_i[1] + aUyy*Rref_i[4] + aUyz*Rref_i[7];
		Uyz = aUyx*Rref_i[2] + aUyy*Rref_i[5] + aUyz*Rref_i[8];

		Uzx = aUzx*Rref_i[0] + aUzy*Rref_i[3] + aUzz*Rref_i[6];
		Uzy = aUzx*Rref_i[1] + aUzy*Rref_i[4] + aUzz*Rref_i[7];
		Uzz = aUzx*Rref_i[2] + aUzy*Rref_i[5] + aUzz*Rref_i[8];

		// Addition of smoothed pressure for particles near mesh
		Particles->p_smooth[i] += sum10/PSystem->pndGradientZero;
		if(Particles->p_smooth[i] < PSystem->epsilonZero) Particles->p_smooth[i] = 0.0;

		// 0.5*(gradU + gradUt)
		Particles->S11[i] = 0.5*(Uxx + Uxx);
		Particles->S12[i] = 0.5*(Uxy + Uyx);
		Particles->S13[i] = 0.5*(Uxz + Uzx);
		Particles->S22[i] = 0.5*(Uyy + Uyy);
		Particles->S23[i] = 0.5*(Uyz + Uzy);
		Particles->S33[i] = 0.5*(Uzz + Uzz);

		//Particles->II[i] = 0.5*Uxx*Uxx + 0.5*Uyy*Uyy + 0.25*(Uxy + Uyx)*(Uxy + Uyx);

		// Addition of II for particles near mesh
		//Particles->II[i] += 0.5*(Particles->S11[i]*Particles->S11[i] + Particles->S12[i]*Particles->S12[i] + Particles->S13[i]*Particles->S13[i] + Particles->S12[i]*Particles->S12[i] + Particles->S22[i]*Particles->S22[i] + Particles->S23[i]*Particles->S23[i] + Particles->S13[i]*Particles->S13[i] + Particles->S23[i]*Particles->S23[i] + Particles->S33[i]*Particles->S33[i]);
		Particles->II[i] += 0.5*(Particles->S11[i]*Particles->S11[i] + 2.0*Particles->S12[i]*Particles->S12[i] + 2.0*Particles->S13[i]*Particles->S13[i] + Particles->S22[i]*Particles->S22[i] + 2.0*Particles->S23[i]*Particles->S23[i] + Particles->S33[i]*Particles->S33[i]);
//		Particles->II[i] = 0.5*(Particles->S11[i]*Particles->S11[i] + Particles->S12[i]*Particles->S12[i] + Particles->S13[i]*Particles->S13[i] + Particles->S12[i]*Particles->S12[i] + Particles->S22[i]*Particles->S22[i] + Particles->S23[i]*Particles->S23[i] + Particles->S13[i]*Particles->S13[i] + Particles->S23[i]*Particles->S23[i] + Particles->S33[i]*Particles->S33[i]);
		//Particles->II[i]= Particles->S11[i]*Particles->S22[i] +Particles->S22[i]*Particles->S33[i]+ Particles->S11[i]*Particles->S33[i] - Particles->S12[i]*Particles->S12[i] -Particles->S13[i]*Particles->S13[i]- Particles->S23[i]*Particles->S23[i] ;
		if(Particles->II[i] < PSystem->epsilonZero || Particles->II[i]*0 != 0) Particles->II[i] = 0.0;
		//II=fabs(Particles->S11[i]*Particles->S22[i]-Particles->S12[i]*Particles->S12[i]);
	}}
	
	// Newtonian viscosity
	if(PSystem->fluidType == viscType::NEWTONIAN)
	{
	// Loop only for particles near mesh
#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<Particles->numParticles; i++) {
		//if(Particles->particleType[i] == PSystem->fluid) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
			if(Particles->PTYPE[i] <= 1)Particles->MEU[i] = PSystem->KNM_VS1 * PSystem->DNS_FL1;
			if(Particles->PTYPE[i] != 1)Particles->MEU[i] = PSystem->KNM_VS2 * PSystem->DNS_FL2;
		}}

//		if(TURB>0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2.0*sqrt(Particles->II[i]);

//			if(NEUt[i] * 0 != 0)  NEUt[i] = 0.0;
//			if(NEUt[i]>1.0)     NEUt[i] = 1.0;
//		}
	}

	// Granular PSystem->fluid
	if(PSystem->fluidType == viscType::NON_NEWTONIAN)
	{
		// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
		// Loop only for particles near mesh
//#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<Particles->numParticles; i++) {
		//if(Particles->particleType[i] == PSystem->fluid) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
//			double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double vel2i = velXi*velXi + velYi*velYi + velZi*velZi;

			if(Particles->PTYPE[i] == 1) {
				Particles->MEU[i] = PSystem->KNM_VS1 * PSystem->DNS_FL1;
			}
			else if(Particles->PTYPE[i] == 2) {
				phi = (Particles->Cv[i] - 0.25)*PSystem->PHI_1/(1.0 - 0.25);
				phi2 = (Particles->Cv[i] - 0.25)*PSystem->PHI_2/(1.0 - 0.25);
				if(Particles->Cv[i] <= 0.25) { phi = PSystem->epsilonZero; phi2 = PSystem->epsilonZero; } // phi close to zero
				// if(Particles->PTYPE[i] <= 0) phi = PSystem->PHI_BED;

				// normal stress calculation (mehcanical pressure)
				Particles->p_rheo_new[i] = Particles->p_smooth[i];

				int kx = int( (posXi - Xmin)/aaL0 ) + 1;
				int ky = int( (posYi - Ymin)/aaL0 ) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress = (BL[kx][ky] - posZi + PSystem->partDist*0.5)*(PSystem->DNS_FL2)*gravityMod;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - posZi + PSystem->partDist*0.5)*(PSystem->DNS_FL2 - PSystem->DNS_FL1)*gravityMod - vel2i*(PSystem->DNS_FL2 - PSystem->DNS_FL1)*0.5;	// normal_stress= Gama.H

				if(Particles->p_smooth[i] < (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod) Particles->p_smooth[i] = (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod;
				if(PSystem->timeCurrent <= 1.0) normal_stress = (1.0 - PSystem->timeCurrent)*(Particles->p_smooth[i] - (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod) + PSystem->timeCurrent*normal_stress;

				// normal_stress = Particles->p_smooth[i] - (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod;
				// normal_stress = Particles->p_smooth[i];
				// normal_stress = normal_stress*0.61*1500/PSystem->DNS_FL2;
				if(normal_stress < 1.0 || Particles->Cv[i] < 0.5) normal_stress = 1.0;

				Particles->p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Particles->Inertia[i] = sqrt(Particles->II[i])*PSystem->DG/sqrt(normal_stress/PSystem->DNS_SDT);		// Free-fall (dry granular material)
				Particles->Inertia[i] = sqrt(Particles->II[i])*PSystem->DG/sqrt(normal_stress/(PSystem->DNS_FL1*PSystem->Cd));	// Grain inertia (submerged)
				//Particles->Inertia[i] = sqrt(Particles->II[i])*(PSystem->KNM_VS1*PSystem->DNS_FL1)/normal_stress ;	// Viscous regime

				// PSystem->VF_max PSystem->VF_min
				Particles->VF[i] = PSystem->VF_max - (PSystem->VF_max - PSystem->VF_min)*Particles->Inertia[i];
				if(Particles->VF[i] < PSystem->VF_min) Particles->VF[i] = PSystem->VF_min;
				Particles->RHO[i] = PSystem->DNS_SDT * Particles->VF[i] + (1.0-Particles->VF[i])*PSystem->DNS_FL1;
				phi = phi * Particles->VF[i] / PSystem->VF_max;

				double yield_stress = PSystem->cohes * cos(phi) + normal_stress * sin(phi);

				if(yield_stress < PSystem->epsilonZero) yield_stress = 0.0;

				double visc_max = (yield_stress*PSystem->mm*0.5 + PSystem->MEU0);

				if(Particles->II[i] > PSystem->epsilonZero)
					Particles->MEU_Y[i] = yield_stress*(1.0 - exp(-PSystem->mm*sqrt(Particles->II[i])))*0.5/sqrt(Particles->II[i]);
				else
					Particles->MEU_Y[i] = visc_max;

				// H-B rheology

				//meu_0 = PSystem->MEU0;

				// Non-linear Meu(I) rheology
				//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*PSystem->DG/(PSystem->I0*sqrt(normal_stress/PSystem->DNS_FL2)+sqrt(Particles->II[i])*PSystem->DG);			//free fall
				meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*PSystem->DG/(PSystem->I0*sqrt(normal_stress/(PSystem->DNS_FL1*PSystem->Cd))+sqrt(Particles->II[i])*PSystem->DG);		//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(PSystem->KNM_VS1*PSystem->DNS_FL1)/(PSystem->I0*normal_stress+sqrt(Particles->II[i])*(PSystem->KNM_VS1*PSystem->DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*PSystem->DG*sqrt(normal_stress*PSystem->DNS_FL2)/PSystem->I0;		//free fall
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*PSystem->DG*sqrt(normal_stress*PSystem->DNS_FL1*PSystem->Cd)/PSystem->I0;	//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*(PSystem->KNM_VS1*PSystem->DNS_FL1)/PSystem->I0;					//viscous

				if(Particles->II[i] <= PSystem->epsilonZero || (meu_0*0) != 0) meu_0 = PSystem->MEU0;

				visc_max = (yield_stress*PSystem->mm*0.5 + meu_0);

				// Herschel bulkley papanastasiou
				Particles->MEU[i] = Particles->MEU_Y[i] + PSystem->MEU0*pow(4.0*Particles->II[i], (PSystem->N - 1.0)*0.5);

				// MEU_Y rheological model
				//Particles->MEU[i] = Particles->MEU_Y[i] + meu_0;
				
				if(Particles->II[i] == 0 || Particles->MEU[i]>visc_max) Particles->MEU[i] = visc_max;
				// if(Particles->PTYPE[i] <= 0) Particles->MEU[i] = Particles->MEU[i]*Particles->Cv[i] + PSystem->DNS_FL1*PSystem->KNM_VS1*(1.0 - Particles->Cv[i]);
			}
			
			if(Particles->PTYPE[i] >= 2) {
				if(Particles->Cv[i] > 0.5) Particles->RHO[i] = PSystem->DNS_FL2;
				else Particles->RHO[i] = Particles->Cv[i]*PSystem->DNS_FL2 + (1.0 - Particles->Cv[i])*PSystem->DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if(stress_cal_method == 2)
//		{
//			for(i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0, sum8 = 0.0, sum9 = 0.0, sum10 = 0.0;
//				for(l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if(i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2.0 * Particles->MEU[i] * Particles->MEU[j] / (Particles->MEU[i] + Particles->MEU[j]);
//						if((NEUt[i] + NEUt[j])>0) meuij = meuij + 2.0 * NEUt[i] * Particles->RHO[i] * NEUt[j] * Particles->RHO[j] / (NEUt[i] * Particles->RHO[i] + NEUt[j] * Particles->RHO[j]);

//						sum1 = sum1 + meuij * (x_vel[j] - x_vel[i])*DX(i, j)*w / d / d;
//						sum2 = sum2 + meuij * (x_vel[j] - x_vel[i])*DY(i, j)*w / d / d;
//						sum3 = sum3 + meuij * (x_vel[j] - x_vel[i])*DZ(i, j)*w / d / d;

//						sum4 = sum4 + meuij * (y_vel[j] - y_vel[i])*DX(i, j)*w / d / d;
//						sum5 = sum5 + meuij * (y_vel[j] - y_vel[i])*DY(i, j)*w / d / d;
//						sum6 = sum6 + meuij * (y_vel[j] - y_vel[i])*DZ(i, j)*w / d / d;

//						sum7 = sum7 + meuij * (z_vel[j] - z_vel[i])*DX(i, j)*w / d / d;
//						sum8 = sum8 + meuij * (z_vel[j] - z_vel[i])*DY(i, j)*w / d / d;
//						sum9 = sum9 + meuij * (z_vel[j] - z_vel[i])*DZ(i, j)*w / d / d;
//					}
//				}

//				Tau_xx[i] = (PSystem->dim / n0) * 2.0 * sum1;
//				Tau_yy[i] = (PSystem->dim / n0) * 2.0 * sum5;
//				Tau_zz[i] = (PSystem->dim / n0) * 2.0 * sum9;

//				Tau_xy[i] = (PSystem->dim / n0)*(sum2 + sum4);
//				Tau_xz[i] = (PSystem->dim / n0)*(sum3 + sum7);
//				Tau_yz[i] = (PSystem->dim / n0)*(sum6 + sum8);
//			}
//		}

	} // if(PSystem->fluidType == viscType::NON_NEWTONIAN)

	//---------------------------------------------------------------

//	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS;// delete[]p_smooth;
//	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL;// p_smooth = NULL;

	delete[]BL; delete[]WL; delete[]PS;
	BL = NULL; WL = NULL; PS = NULL;

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// No-Slip condition. Viscosity interaction values
void MpsViscosity::calcWallNoSlipViscosityInteractionVal(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {

	double gravityMod = sqrt(PSystem->gravityX*PSystem->gravityX + PSystem->gravityY*PSystem->gravityY + PSystem->gravityZ*PSystem->gravityZ);

	//double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;//, grain_VF, *p_smooth;
	double phi = 0.0, phi2 = 0.0, meu_0, normal_stress;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
	// Xmin = 0.0 - PSystem->partDist*3.0; Xmax = 1.65 + PSystem->partDist*3.0;
	// Ymin = 0.0 - PSystem->partDist*3.0; Ymax = 0.15 + PSystem->partDist*3.0;
	// Zmin = 0.0 - PSystem->partDist*3.0; //Zmax = 0.7 + PSystem->partDist*30.0;
	// damErosion3D
	// Xmin = 0.0 - PSystem->partDist*3.0; Xmax = 2.00 + PSystem->partDist*3.0;
	// Ymin = 0.0 - PSystem->partDist*3.0; Ymax = 0.10 + PSystem->partDist*3.0;
	// Changed !!!

	Xmin = PSystem->domainMinX;
	Xmax = PSystem->domainMaxX;
	Ymin = PSystem->domainMinY;
	Ymax = Ymin;	Zmin = 0.0;
	if ((int)PSystem->dim == 3) {
		Ymax = PSystem->domainMaxY;
		Zmin = PSystem->domainMinZ;
	}

	// Search free-surface particles for each interval of aa = 2 particles in wall
	double aaL0 = 2.0*PSystem->partDist;
	int kx_max = int( (Xmax - Xmin)/aaL0 ) + 1;
	int ky_max = int( (Ymax - Ymin)/aaL0 ) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
/*
	S11 = new double[Particles->numParticles + 1];
	S22 = new double[Particles->numParticles + 1];
	S33 = new double[Particles->numParticles + 1];
	S12 = new double[Particles->numParticles + 1];
	S13 = new double[Particles->numParticles + 1];
	S23 = new double[Particles->numParticles + 1];
*/
	//p_smooth = new double[Particles->numParticles + 1];
	BL = new double*[kx_max + 1];  // bed level
	WL = new double*[kx_max + 1];  // water level
	PS = new double*[kx_max + 1];  // pressure sediment

#pragma omp parallel for
	for(int m = 1; m <= kx_max; m++)
	{
		BL[m] = new double[ky_max + 1];
		WL[m] = new double[ky_max + 1];
		PS[m] = new double[ky_max + 1];
	}

	// Determining the bed level
#pragma omp parallel for schedule(dynamic,64)
	for(int kx = 1; kx <= kx_max; kx++)
	{
		for(int ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->fluid) {
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		
			int kx = int( (posXi - Xmin)/aaL0 ) + 1;
			int ky = int( (posYi - Ymin)/aaL0 ) + 1;

			//if(posZi>BL[kx][ky] && Particles->Cv[i]>0.5) { BL[kx][ky] = posZi; PS[kx][ky] = pnew[i]; }
			if(posZi>BL[kx][ky] && Particles->Cv[i]>0.5) { BL[kx][ky] = posZi; PS[kx][ky] = Particles->press[i]; }
			if(posZi>WL[kx][ky] && Particles->PTYPE[i] == 1) { WL[kx][ky] = posZi; }
		}
	}

	// Strain rate calculation
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
	//if(Particles->particleType[i] == PSystem->fluid) {
	if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
		double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0, sum8 = 0.0, sum9 = 0.0, sum10 = 0.0;

//		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2];

		// Transformation matrix R_i = I
		double Rref_i[9], Rinv_i[9], normaliw[3], normalMod2;
	    // normal PSystem->fluid-wall particle = 0.5*(normal PSystem->fluid-mirror particle)
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
		// Mirror particle velocity vi' = Ri_inv * [vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}] 
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

					double vec_mijx = Particles->vel[j*3  ]-velMirrorXi;	
					double vec_mijy = Particles->vel[j*3+1]-velMirrorYi;	
					double vec_mijz = Particles->vel[j*3+2]-velMirrorZi;

					sum1 += vec_mijx*v0imj*wS/dstimj2;
					sum2 += vec_mijx*v1imj*wS/dstimj2;
					sum3 += vec_mijx*v2imj*wS/dstimj2;
					
					sum4 += vec_mijy*v0imj*wS/dstimj2;
					sum5 += vec_mijy*v1imj*wS/dstimj2;
					sum6 += vec_mijy*v2imj*wS/dstimj2;
					
					sum7 += vec_mijz*v0imj*wS/dstimj2;
					sum8 += vec_mijz*v1imj*wS/dstimj2;
					sum9 += vec_mijz*v2imj*wS/dstimj2;
					
					sum10 += Particles->press[j]*wS;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Rinv_i * gradU * Rref_i = - gradU * Rref_i
		// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/PSystem->pndGradientZero)
		Uxx = PSystem->coeffPressGrad*(sum1*Rref_i[0] + sum2*Rref_i[3] + sum3*Rref_i[6]);
		Uxy = PSystem->coeffPressGrad*(sum1*Rref_i[1] + sum2*Rref_i[4] + sum3*Rref_i[7]);
		Uxz = PSystem->coeffPressGrad*(sum1*Rref_i[2] + sum2*Rref_i[5] + sum3*Rref_i[8]);

		Uyx = PSystem->coeffPressGrad*(sum4*Rref_i[0] + sum5*Rref_i[3] + sum6*Rref_i[6]);
		Uyy = PSystem->coeffPressGrad*(sum4*Rref_i[1] + sum5*Rref_i[4] + sum6*Rref_i[7]);
		Uyz = PSystem->coeffPressGrad*(sum4*Rref_i[2] + sum5*Rref_i[5] + sum6*Rref_i[8]);

		Uzx = PSystem->coeffPressGrad*(sum7*Rref_i[0] + sum8*Rref_i[3] + sum9*Rref_i[6]);
		Uzy = PSystem->coeffPressGrad*(sum7*Rref_i[1] + sum8*Rref_i[4] + sum9*Rref_i[7]);
		Uzz = PSystem->coeffPressGrad*(sum7*Rref_i[2] + sum8*Rref_i[5] + sum9*Rref_i[8]);

		// Addition of smoothed pressure for particles near mesh
		Particles->p_smooth[i] += sum10/PSystem->pndGradientZero;
		if(Particles->p_smooth[i] < PSystem->epsilonZero) Particles->p_smooth[i] = 0.0;

		// - (Rref_i * gradU) - (Rref_i * gradU)t
		Particles->S11[i] = 0.5*(Uxx + Uxx);
		Particles->S12[i] = 0.5*(Uxy + Uyx);
		Particles->S13[i] = 0.5*(Uxz + Uzx);
		Particles->S22[i] = 0.5*(Uyy + Uyy);
		Particles->S23[i] = 0.5*(Uyz + Uzy);
		Particles->S33[i] = 0.5*(Uzz + Uzz);

		//Particles->II[i] = 0.5*Uxx*Uxx + 0.5*Uyy*Uyy + 0.25*(Uxy + Uyx)*(Uxy + Uyx);
		
		// Addition of II for particles near mesh
		//Particles->II[i] += 0.5*(Particles->S11[i]*Particles->S11[i] + Particles->S12[i]*Particles->S12[i] + Particles->S13[i]*Particles->S13[i] + Particles->S12[i]*Particles->S12[i] + Particles->S22[i]*Particles->S22[i] + Particles->S23[i]*Particles->S23[i] + Particles->S13[i]*Particles->S13[i] + Particles->S23[i]*Particles->S23[i] + Particles->S33[i]*Particles->S33[i]);
		Particles->II[i] += 0.5*(Particles->S11[i]*Particles->S11[i] + 2.0*Particles->S12[i]*Particles->S12[i] + 2.0*Particles->S13[i]*Particles->S13[i] + Particles->S22[i]*Particles->S22[i] + 2.0*Particles->S23[i]*Particles->S23[i] + Particles->S33[i]*Particles->S33[i]);
//		Particles->II[i] = 0.5*(Particles->S11[i]*Particles->S11[i] + Particles->S12[i]*Particles->S12[i] + Particles->S13[i]*Particles->S13[i] + Particles->S12[i]*Particles->S12[i] + Particles->S22[i]*Particles->S22[i] + Particles->S23[i]*Particles->S23[i] + Particles->S13[i]*Particles->S13[i] + Particles->S23[i]*Particles->S23[i] + Particles->S33[i]*Particles->S33[i]);
		//Particles->II[i]= Particles->S11[i]*Particles->S22[i] +Particles->S22[i]*Particles->S33[i]+ Particles->S11[i]*Particles->S33[i] - Particles->S12[i]*Particles->S12[i] -Particles->S13[i]*Particles->S13[i]- Particles->S23[i]*Particles->S23[i] ;
		if(Particles->II[i] < PSystem->epsilonZero || Particles->II[i]*0 != 0) Particles->II[i] = 0.0;
		//II=fabs(Particles->S11[i]*Particles->S22[i]-Particles->S12[i]*Particles->S12[i]);
	}}
	
	// Newtonian viscosity
	if(PSystem->fluidType == viscType::NEWTONIAN)
	{
	// Loop only for particles near mesh
#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<Particles->numParticles; i++) {
		//if(Particles->particleType[i] == PSystem->fluid) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
			if(Particles->PTYPE[i] <= 1)Particles->MEU[i] = PSystem->KNM_VS1 * PSystem->DNS_FL1;
			if(Particles->PTYPE[i] != 1)Particles->MEU[i] = PSystem->KNM_VS2 * PSystem->DNS_FL2;
		}}

//		if(TURB>0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2.0*sqrt(Particles->II[i]);

//			if(NEUt[i] * 0 != 0)  NEUt[i] = 0.0;
//			if(NEUt[i]>1.0)     NEUt[i] = 1.0;
//		}
	}

	// Granular PSystem->fluid
	if(PSystem->fluidType == viscType::NON_NEWTONIAN)
	{
		// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
		// Loop only for particles near mesh
//#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<Particles->numParticles; i++) {
		//if(Particles->particleType[i] == PSystem->fluid) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
//			double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double vel2i = velXi*velXi + velYi*velYi + velZi*velZi;

			if(Particles->PTYPE[i] == 1) 
				Particles->MEU[i] = PSystem->KNM_VS1 * PSystem->DNS_FL1;
			else if(Particles->PTYPE[i] == 2) {
				phi = (Particles->Cv[i] - 0.25)*PSystem->PHI_1/(1.0 - 0.25);
				phi2 = (Particles->Cv[i] - 0.25)*PSystem->PHI_2/(1.0 - 0.25);
				if(Particles->Cv[i] <= 0.25) { phi = PSystem->epsilonZero; phi2 = PSystem->epsilonZero; } // phi close to zero
				// if(Particles->PTYPE[i] <= 0) phi = PSystem->PHI_BED;

				// normal stress calculation (mechanical pressure)
				Particles->p_rheo_new[i] = Particles->p_smooth[i];

				int kx = int( (posXi - Xmin)/aaL0 ) + 1;
				int ky = int( (posYi - Ymin)/aaL0 ) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress = (BL[kx][ky] - posZi + PSystem->partDist*0.5)*(PSystem->DNS_FL2)*gravityMod;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - posZi + PSystem->partDist*0.5)*(PSystem->DNS_FL2 - PSystem->DNS_FL1)*gravityMod - vel2i*(PSystem->DNS_FL2 - PSystem->DNS_FL1)*0.5;	// normal_stress= Gama.H

				if(Particles->p_smooth[i] < (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod) Particles->p_smooth[i] = (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod;
				if(PSystem->timeCurrent <= 1.0) normal_stress = (1.0 - PSystem->timeCurrent)*(Particles->p_smooth[i] - (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod) + PSystem->timeCurrent*normal_stress;

				// normal_stress = Particles->p_smooth[i] - (WL[kx][ky] - posZi)*PSystem->DNS_FL1*gravityMod;
				// normal_stress = Particles->p_smooth[i];
				// normal_stress = normal_stress*0.61*1500/PSystem->DNS_FL2;
				if(normal_stress < 1.0 || Particles->Cv[i] < 0.5) normal_stress = 1.0;

				Particles->p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Particles->Inertia[i] = sqrt(Particles->II[i])*PSystem->DG/sqrt(normal_stress/PSystem->DNS_SDT);		// Free-fall (dry granular material)
				Particles->Inertia[i] = sqrt(Particles->II[i])*PSystem->DG/sqrt(normal_stress/(PSystem->DNS_FL1*PSystem->Cd));	// Grain inertia (submerged)
				//Particles->Inertia[i] = sqrt(Particles->II[i])*(PSystem->KNM_VS1*PSystem->DNS_FL1)/normal_stress ;	// Viscous regime

				// PSystem->VF_max PSystem->VF_min
				Particles->VF[i] = PSystem->VF_max - (PSystem->VF_max - PSystem->VF_min)*Particles->Inertia[i];
				if(Particles->VF[i] < PSystem->VF_min) Particles->VF[i] = PSystem->VF_min;
				Particles->RHO[i] = PSystem->DNS_SDT*Particles->VF[i] + (1.0-Particles->VF[i])*PSystem->DNS_FL1;
				phi = phi*Particles->VF[i]/PSystem->VF_max;

				double yield_stress = PSystem->cohes*cos(phi) + normal_stress*sin(phi);

				if(yield_stress < 0.0) yield_stress = 0.0;

				double visc_max = (yield_stress*PSystem->mm*0.5 + PSystem->MEU0);

				if(Particles->II[i] > PSystem->epsilonZero)
					Particles->MEU_Y[i] = yield_stress*(1.0 - exp(-PSystem->mm*sqrt(Particles->II[i])))*0.5/sqrt(Particles->II[i]);
				else
					Particles->MEU_Y[i] = visc_max;

				// H-B rheology

				//meu_0 = PSystem->MEU0;

				// Non-linear Meu(I) rheology
				//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*PSystem->DG/(PSystem->I0*sqrt(normal_stress/PSystem->DNS_FL2)+sqrt(Particles->II[i])*PSystem->DG);			//free fall
				meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*PSystem->DG/(PSystem->I0*sqrt(normal_stress/(PSystem->DNS_FL1*PSystem->Cd))+sqrt(Particles->II[i])*PSystem->DG);		//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(PSystem->KNM_VS1*PSystem->DNS_FL1)/(PSystem->I0*normal_stress+sqrt(Particles->II[i])*(PSystem->KNM_VS1*PSystem->DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*PSystem->DG*sqrt(normal_stress*PSystem->DNS_FL2)/PSystem->I0;		//free fall
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*PSystem->DG*sqrt(normal_stress*PSystem->DNS_FL1*PSystem->Cd)/PSystem->I0;	//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*(PSystem->KNM_VS1*PSystem->DNS_FL1)/PSystem->I0;					//viscous

				if(Particles->II[i] <= PSystem->epsilonZero || (meu_0*0) != 0) meu_0 = PSystem->MEU0;

				visc_max = (yield_stress*PSystem->mm*0.5 + meu_0);

				// Herschel bulkley papanastasiou
				Particles->MEU[i] = Particles->MEU_Y[i] + PSystem->MEU0*pow(4.0*Particles->II[i], (PSystem->N - 1.0)*0.5);

				// MEU_Y rheological model
				//Particles->MEU[i] = Particles->MEU_Y[i] + meu_0;
				
				if(Particles->II[i] == 0 || Particles->MEU[i]>visc_max) Particles->MEU[i] = visc_max;
				// if(Particles->PTYPE[i] <= 0) Particles->MEU[i] = Particles->MEU[i]*Particles->Cv[i] + PSystem->DNS_FL1*PSystem->KNM_VS1*(1.0 - Particles->Cv[i]);
			}
			
			if(Particles->PTYPE[i] >= 2) {
				if(Particles->Cv[i] > 0.5) Particles->RHO[i] = PSystem->DNS_FL2;
				else Particles->RHO[i] = Particles->Cv[i]*PSystem->DNS_FL2 + (1.0 - Particles->Cv[i])*PSystem->DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if(stress_cal_method == 2)
//		{
//			for(i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0, sum8 = 0.0, sum9 = 0.0, sum10 = 0.0;
//				for(l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if(i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2.0 * Particles->MEU[i] * Particles->MEU[j] / (Particles->MEU[i] + Particles->MEU[j]);
//						if((NEUt[i] + NEUt[j])>0) meuij = meuij + 2.0 * NEUt[i] * Particles->RHO[i] * NEUt[j] * Particles->RHO[j] / (NEUt[i] * Particles->RHO[i] + NEUt[j] * Particles->RHO[j]);

//						sum1 = sum1 + meuij * (x_vel[j] - x_vel[i])*DX(i, j)*w / d / d;
//						sum2 = sum2 + meuij * (x_vel[j] - x_vel[i])*DY(i, j)*w / d / d;
//						sum3 = sum3 + meuij * (x_vel[j] - x_vel[i])*DZ(i, j)*w / d / d;

//						sum4 = sum4 + meuij * (y_vel[j] - y_vel[i])*DX(i, j)*w / d / d;
//						sum5 = sum5 + meuij * (y_vel[j] - y_vel[i])*DY(i, j)*w / d / d;
//						sum6 = sum6 + meuij * (y_vel[j] - y_vel[i])*DZ(i, j)*w / d / d;

//						sum7 = sum7 + meuij * (z_vel[j] - z_vel[i])*DX(i, j)*w / d / d;
//						sum8 = sum8 + meuij * (z_vel[j] - z_vel[i])*DY(i, j)*w / d / d;
//						sum9 = sum9 + meuij * (z_vel[j] - z_vel[i])*DZ(i, j)*w / d / d;
//					}
//				}

//				Tau_xx[i] = (PSystem->dim / n0) * 2.0 * sum1;
//				Tau_yy[i] = (PSystem->dim / n0) * 2.0 * sum5;
//				Tau_zz[i] = (PSystem->dim / n0) * 2.0 * sum9;

//				Tau_xy[i] = (PSystem->dim / n0)*(sum2 + sum4);
//				Tau_xz[i] = (PSystem->dim / n0)*(sum3 + sum7);
//				Tau_yz[i] = (PSystem->dim / n0)*(sum6 + sum8);
//			}
//		}

	} // if(PSystem->fluidType == viscType::NON_NEWTONIAN)

	//---------------------------------------------------------------

//	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS; //delete[]p_smooth;
//	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL;// p_smooth = NULL;
	
	delete[]BL; delete[]WL; delete[]PS;
	BL = NULL; WL = NULL; PS = NULL;

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Free-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
void MpsViscosity::calcWallSlipViscosity(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	double VolumeForce = pow(PSystem->partDist,PSystem->dim);
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
	//if(Particles->particleType[i] == PSystem->fluid) {
	if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
		
		double meu_i = Particles->MEU[i];
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2

		// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
		double Rref_i[9], normaliw[3], normalMod2;
		// Normal PSystem->fluid-wall particle = 0.5*(normal PSystem->fluid-mirror particle)
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
				if(dstij2 < PSystem->reL2 && dstimj2 < PSystem->reL2 && dstij2 < dstimj2) {
				if(j != i) {
					double dst = sqrt(dstimj2);
					double wL = Particles->weight(dst, PSystem->reL, PSystem->weightType);
					double neu_ij;
					if((meu_i + Particles->MEU[j]) > PSystem->epsilonZero)
						neu_ij = 2.0 * meu_i * Particles->MEU[j] / (meu_i + Particles->MEU[j]);
					else
						neu_ij = 0.0;

//neu_ij = PSystem->KNM_VS2 * PSystem->DNS_FL2;

					if(Particles->PTYPE[i] == 1) neu_ij = neu_ij/PSystem->DNS_FL1;
					else neu_ij = neu_ij/PSystem->DNS_FL2;

					//if((NEUt[i] + NEUt[j])>0) neu_ij = neu_ij + (2.0 * NEUt[i] * Particles->RHO[j] * NEUt[j] * Particles->RHO[j] / (NEUt[i] * Particles->RHO[i] + NEUt[j] * Particles->RHO[j])) / Particles->RHO[i];

					// Original
//					accX +=(Particles->vel[j*3  ]-velMirrorXi)*w;
//					accY +=(Particles->vel[j*3+1]-velMirrorYi)*w;
//					accZ +=(Particles->vel[j*3+2]-velMirrorZi)*w;
					// Modified
					accX +=(Particles->vel[j*3  ]-velMirrorXi)*wL*neu_ij;
					accY +=(Particles->vel[j*3+1]-velMirrorYi)*wL*neu_ij;
					accZ +=(Particles->vel[j*3+2]-velMirrorZi)*wL*neu_ij;

					//accX +=(Velk[j*3  ]-velMirrorXi)*w;
					//accY +=(Velk[j*3+1]-velMirrorYi)*w;
					//accZ +=(Velk[j*3+2]-velMirrorZi)*w;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Add "i" contribution ("i" is a neighbor of "mirror i")
		double v0imi, v1imi, v2imi, dstimi2;
		Particles->sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
		
		if(dstimi2 < PSystem->reL2) {
			double dst = sqrt(dstimi2);
			double wL = Particles->weight(dst, PSystem->reL, PSystem->weightType);
			double neu_ij;
			if(meu_i > PSystem->epsilonZero)
				neu_ij = 2.0 * meu_i * meu_i / (meu_i + meu_i);
			else
				neu_ij = 0.0;

//neu_ij = PSystem->KNM_VS2 * PSystem->DNS_FL2;

			if(Particles->PTYPE[i] == 1) neu_ij = neu_ij/PSystem->DNS_FL1;
			else neu_ij = neu_ij/PSystem->DNS_FL2;

			// Original
//			accX +=(velXi-velMirrorXi)*w;
//			accY +=(velYi-velMirrorYi)*w;
//			accZ +=(velZi-velMirrorZi)*w;

			// Modified
			accX +=(velXi-velMirrorXi)*wL*neu_ij;
			accY +=(velYi-velMirrorYi)*wL*neu_ij;
			accZ +=(velZi-velMirrorZi)*wL*neu_ij;
		}

		// Wall laplacian Mitsume`s model
		// Correction of velocity
		// Original
//      acc[i*3  ] += (Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*PSystem->coeffViscosity;
//		acc[i*3+1] += (Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*PSystem->coeffViscosity;
//		acc[i*3+2] += (Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*PSystem->coeffViscosity;

		// coeffViscMultiphase = 2.0*PSystem->dim/(PSystem->pndLargeZero*PSystem->lambdaZero);
		// Modified
		Particles->acc[i*3  ] += (Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*PSystem->coeffViscMultiphase;
		Particles->acc[i*3+1] += (Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*PSystem->coeffViscMultiphase;
		Particles->acc[i*3+2] += (Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*PSystem->coeffViscMultiphase;

		
		// FSI
		// Force on wall
		Particles->forceWall[i*3  ] += - (Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*PSystem->coeffViscMultiphase*VolumeForce*Particles->RHO[i];
		Particles->forceWall[i*3+1] += - (Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*PSystem->coeffViscMultiphase*VolumeForce*Particles->RHO[i];
		Particles->forceWall[i*3+2] += - (Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*PSystem->coeffViscMultiphase*VolumeForce*Particles->RHO[i];
	}}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// No-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
void MpsViscosity::calcWallNoSlipViscosity(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	double VolumeForce = pow(PSystem->partDist,PSystem->dim);
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
	//if(Particles->particleType[i] == PSystem->fluid) {
	if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {

		double meu_i = Particles->MEU[i];
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2];

		// Inverse matrix Rinv_i = - I
		double Rinv_i[9], normaliw[3], normalMod2;
		// normal PSystem->fluid-wall particle = 0.5*(normal PSystem->fluid-mirror particle)
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
				if(dstij2 < PSystem->reL2 && dstimj2 < PSystem->reL2 && dstij2 < dstimj2) {
				if(j != i) {
					double dst = sqrt(dstimj2);
					double wL = Particles->weight(dst, PSystem->reL, PSystem->weightType);
					double neu_ij;
					if((meu_i + Particles->MEU[j]) > PSystem->epsilonZero)
						neu_ij = 2.0 * meu_i * Particles->MEU[j] / (meu_i + Particles->MEU[j]);
					else
						neu_ij = 0.0;

//neu_ij = PSystem->KNM_VS2 * PSystem->DNS_FL2;


					if(Particles->PTYPE[i] == 1) neu_ij = neu_ij/PSystem->DNS_FL1;
					else neu_ij = neu_ij/PSystem->DNS_FL2;
					
					//if((NEUt[i] + NEUt[j])>0) neu_ij = neu_ij + (2.0 * NEUt[i] * Particles->RHO[j] * NEUt[j] * Particles->RHO[j] / (NEUt[i] * Particles->RHO[i] + NEUt[j] * Particles->RHO[j])) / Particles->RHO[i];

					// Original
//					accX +=(Particles->vel[j*3  ]-velMirrorXi)*w;
//					accY +=(Particles->vel[j*3+1]-velMirrorYi)*w;
//					accZ +=(Particles->vel[j*3+2]-velMirrorZi)*w;
					// Modified
					accX +=(Particles->vel[j*3  ]-velMirrorXi)*wL*neu_ij;
					accY +=(Particles->vel[j*3+1]-velMirrorYi)*wL*neu_ij;
					accZ +=(Particles->vel[j*3+2]-velMirrorZi)*wL*neu_ij;
					
					//accX +=(Velk[j*3  ]-velMirrorXi)*w;
					//accY +=(Velk[j*3+1]-velMirrorYi)*w;
					//accZ +=(Velk[j*3+2]-velMirrorZi)*w;

					//if(i==2817) {
					//	Fwall[i*3  ] += 1;//AA[0];
					//	Fwall[i*3+1] += 1;//AA[1];
					//	Fwall[i*3+2] += 1;//AA[2];
					//	std::cout << j << " " << Velk[j*3] << " " << Velk[j*3+1] << " " << Velk[j*3+2] << " Pj " << Particles->press[j] << std::endl;
					//}
					//accX += (Rinv_i[0]*(Velk[j*3  ]-velMirrorXi)+ Rinv_i[1]*(Velk[j*3+1]-velMirrorYi) + Rinv_i[2]*(Velk[j*3+2]-velMirrorZi))*w;
					//accY += (Rinv_i[3]*(Velk[j*3  ]-velMirrorXi)+ Rinv_i[4]*(Velk[j*3+1]-velMirrorYi) + Rinv_i[5]*(Velk[j*3+2]-velMirrorZi))*w;
					//accZ += (Rinv_i[6]*(Velk[j*3  ]-velMirrorXi)+ Rinv_i[7]*(Velk[j*3+1]-velMirrorYi) + Rinv_i[8]*(Velk[j*3+2]-velMirrorZi))*w;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Add "i" contribution ("i" is a neighbor of "mirror i")
		double v0imi, v1imi, v2imi, dstimi2;
		Particles->sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
		
		if(dstimi2 < PSystem->reL2) {
			double dst = sqrt(dstimi2);
			double wL = Particles->weight(dst, PSystem->reL, PSystem->weightType);
			double neu_ij;
			if(meu_i > PSystem->epsilonZero)
				neu_ij = 2.0 * meu_i * meu_i / (meu_i + meu_i);
			else
				neu_ij = 0.0;

//neu_ij = PSystem->KNM_VS2 * PSystem->DNS_FL2;

			if(Particles->PTYPE[i] == 1) neu_ij = neu_ij/PSystem->DNS_FL1;
			else neu_ij = neu_ij/PSystem->DNS_FL2;

			// Original
//			accX +=(velXi-velMirrorXi)*w;
//			accY +=(velYi-velMirrorYi)*w;
//			accZ +=(velZi-velMirrorZi)*w;

			// Modified
			accX +=(velXi-velMirrorXi)*wL*neu_ij;
			accY +=(velYi-velMirrorYi)*wL*neu_ij;
			accZ +=(velZi-velMirrorZi)*wL*neu_ij;

			//accX += (Rinv_i[0]*(velXi-velMirrorXi)+ Rinv_i[1]*(velYi-velMirrorYi) + Rinv_i[2]*(velZi-velMirrorZi))*w;
			//accY += (Rinv_i[3]*(velXi-velMirrorXi)+ Rinv_i[4]*(velYi-velMirrorYi) + Rinv_i[5]*(velZi-velMirrorZi))*w;
			//accZ += (Rinv_i[6]*(velXi-velMirrorXi)+ Rinv_i[7]*(velYi-velMirrorYi) + Rinv_i[8]*(velZi-velMirrorZi))*w;
		}

		//if(i==6) {
		//	printf("Vx:%lf Vy:%lf Vz:%lf Vx:%lf Vy:%lf Vz:%lf\n",velXi,velMirrorXi,velYi,velMirrorYi,velZi,velMirrorZi);
		//	printf("Accx:%lf Bccy:%lf Bccz:%lf \n", acc[i*3],acc[i*3+1],acc[i*3+2]);
			//printf("Bccx:%lf Bccy:%lf Bccz:%lf \n", acc[i*3],acc[i*3+1],acc[i*3+2]);
		//}
		// Wall laplacian Mitsume`s model
		// Correction of velocity
		// coeffViscMultiphase = 2.0*PSystem->dim/(PSystem->pndLargeZero*PSystem->lambdaZero);
		// Original
//     	acc[i*3  ] += (Rinv_i[0]*accX + Rinv_i[1]*accY + Rinv_i[2]*accZ)*PSystem->coeffViscosity;
//		acc[i*3+1] += (Rinv_i[3]*accX + Rinv_i[4]*accY + Rinv_i[5]*accZ)*PSystem->coeffViscosity;
//		acc[i*3+2] += (Rinv_i[6]*accX + Rinv_i[7]*accY + Rinv_i[8]*accZ)*PSystem->coeffViscosity;
		// Modified
		Particles->acc[i*3  ] += (Rinv_i[0]*accX + Rinv_i[1]*accY + Rinv_i[2]*accZ)*PSystem->coeffViscMultiphase;
		Particles->acc[i*3+1] += (Rinv_i[3]*accX + Rinv_i[4]*accY + Rinv_i[5]*accZ)*PSystem->coeffViscMultiphase;
		Particles->acc[i*3+2] += (Rinv_i[6]*accX + Rinv_i[7]*accY + Rinv_i[8]*accZ)*PSystem->coeffViscMultiphase;
		//Acv[i*3  ] = (Rinv_i[0]*accX + Rinv_i[1]*accY + Rinv_i[2]*accZ)*PSystem->coeffViscosity;
		//Acv[i*3+1] = (Rinv_i[3]*accX + Rinv_i[4]*accY + Rinv_i[5]*accZ)*PSystem->coeffViscosity;
		//Acv[i*3+2] = (Rinv_i[6]*accX + Rinv_i[7]*accY + Rinv_i[8]*accZ)*PSystem->coeffViscosity;

		//double AA[3];
		
		//AA[0] = (Rinv_i[0]*accX + Rinv_i[1]*accY + Rinv_i[2]*accZ)*PSystem->coeffViscosity;
		//AA[1] = (Rinv_i[3]*accX + Rinv_i[4]*accY + Rinv_i[5]*accZ)*PSystem->coeffViscosity;
		//AA[2] = (Rinv_i[6]*accX + Rinv_i[7]*accY + Rinv_i[8]*accZ)*PSystem->coeffViscosity;

		//wallParticleForce1[i*3  ] = AA[0];
		//wallParticleForce1[i*3+1] = AA[1];
		//wallParticleForce1[i*3+2] = AA[2];

		//if(i==2817) {
			//Fwall[i*3  ] = AA[0];
			//Fwall[i*3+1] = AA[1];
			//Fwall[i*3+2] = AA[2];
			//std::cout << "t: " << PSystem->timeCurrent << std::endl;
			//std::cout << "Fwall " << AA[0] << " " << AA[1] << " " << AA[2] << std::endl;
			//std::cout << "Veli " << velXi << " " << velYi << " " << velZi << std::endl;
			//std::cout << "Posmi " << posMirrorXi << " " << posMirrorYi << " " << posMirrorZi << std::endl;
			//std::cout << "pndi " << pndi[i] << " Pi " << Particles->press[i] << std::endl;
			//printf("Accx:%lf Accy:%lf Accz:%lf \n", AA[0],AA[1],AA[2]);
		//}
		//if(i==6) {
			//printf("Accx:%lf Accy:%lf Accz:%lf \n", acc[i*3],acc[i*3+1],acc[i*3+2]);
		//	printf("Time:%e\n", PSystem->timeCurrent);
		//	printf("Xi:%e %e %e Xm:%e %e %e\n", posXi,posYi,posZi,posMirrorXi,posMirrorYi,posMirrorZi);
		//	printf("Vi:%e %e %e Vm:%e %e %e\n", velXi,velYi,velZi,velMirrorXi,velMirrorYi,velMirrorZi);
		//	printf("acc:%e %e %e\n", AA[0],AA[1],AA[2]);
		//}


		// FSI
		// Force on wall
		Particles->forceWall[i*3  ] += - (Rinv_i[0]*accX + Rinv_i[1]*accY + Rinv_i[2]*accZ)*PSystem->coeffViscMultiphase*VolumeForce*Particles->RHO[i];
		Particles->forceWall[i*3+1] += - (Rinv_i[3]*accX + Rinv_i[4]*accY + Rinv_i[5]*accZ)*PSystem->coeffViscMultiphase*VolumeForce*Particles->RHO[i];
		Particles->forceWall[i*3+2] += - (Rinv_i[6]*accX + Rinv_i[7]*accY + Rinv_i[8]*accZ)*PSystem->coeffViscMultiphase*VolumeForce*Particles->RHO[i];
	}}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}