// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <iostream>		///< cerr
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "MpsPressure.h"

using namespace std;

// Constructor declaration
MpsPressure::MpsPressure()
{
}
// Destructor declaration
MpsPressure::~MpsPressure()
{
}

// Allocates memory for pressure and source term Vector
void MpsPressure::allocateMemory(MpsParticle *Particles) {
	Particles->pressurePPE = Eigen::VectorXd::Zero(Particles->numParticles);
	Particles->sourceTerm = Eigen::VectorXd::Zero(Particles->numParticles);
}


// Select how to compute pressure
void MpsPressure::calcPress(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets){

	if(PSystem->mpsType == calcPressType::EXPLICIT) {
		calcPressEMPS(PSystem, Particles);
	}
	else if(PSystem->mpsType == calcPressType::WEAKLY) {
		calcPressWCMPS(PSystem, Particles);
	}
	else if(PSystem->mpsType == calcPressType::IMPLICIT_PND)
	{
		solvePressurePoissonPnd(PSystem, Particles, Buckets);
	}
	else if(PSystem->mpsType == calcPressType::IMPLICIT_PND_DIVU)
	{
		calcVelDivergence(PSystem, Particles, Buckets);
		if(PSystem->wallType == boundaryWallType::POLYGON) {
			if(PSystem->slipCondition == slipBC::FREE_SLIP) {
				calcWallSlipVelDivergence(PSystem, Particles, Buckets); // Free-Slip condition
			}
			else if(PSystem->slipCondition == slipBC::NO_SLIP) {
				calcWallNoSlipVelDivergence(PSystem, Particles, Buckets); // No-Slip condition
			}
		}
		solvePressurePoissonPndDivU(PSystem, Particles, Buckets);
	}

}


// Compute pressure EMPS (mpsType = calcPressType::EXPLICIT)
void MpsPressure::calcPressEMPS(MpsParticleSystem *PSystem, MpsParticle *Particles) {
// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		double mi;
		if(Particles->PTYPE[i] == 1) mi = PSystem->DNS_FL1;
		else mi = PSystem->DNS_FL2;
		//if(Particles->particleType[i]==fluid)
		//	mi = Dns[partType::FLUID];
		//else
		//	mi = Dns[partType::WALL];

		double pressure = 0.0;
		if(Particles->particleBC[i] == PSystem->inner) {
			pressure = (Particles->pndi[i] - PSystem->pndSmallZero) * PSystem->coeffPressEMPS * mi;
		}

//		if(Particles->pndSmall[i] < PSystem->pndThreshold*PSystem->pndSmallZero && numNeigh[i] < PSystem->neighThreshold*PSystem->numNeighZero)
//			Particles->particleBC[i] = PSystem->surface;
//		else
//		{
//			Particles->particleBC[i] = PSystem->inner;
//			pressure = (Particles->pndi[i] - PSystem->pndSmallZero) * PSystem->coeffPressEMPS * mi;
//		}
//		if(PSystem->wallType == boundaryWallType::POLYGON) {
//			if(Particles->pndi[i] > PSystem->pndThreshold*PSystem->pndSmallZero) {
//				pressure = -mi*PSystem->gravityZ*(0.3-posZi);
//			}
//		}
//		else if(PSystem->wallType == boundaryWallType::PARTICLE) {
//			if(Particles->pndi[i] > PSystem->pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero) {
//				pressure = -mi*PSystem->gravityZ*(0.3-posZi);
//			}
//		}

		if(pressure < 0.0) {
			pressure = 0.0;
		}
		Particles->press[i] = pressure;
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Compute pressure WCMPS (mpsType = calcPressType::WEAKLY)
void MpsPressure::calcPressWCMPS(MpsParticleSystem *PSystem, MpsParticle *Particles) {
// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		double mi;
		if(Particles->PTYPE[i] == 1) mi = PSystem->DNS_FL1;
		else mi = PSystem->DNS_FL2;
		//if(Particles->particleType[i]==PSystem->fluid)
		//	mi = Dns[partType::FLUID];
		//else
		//	mi = Dns[partType::WALL];

		double pressure = 0.0;
		if(Particles->particleBC[i] == PSystem->inner) {
		// if(Particles->particleBC[i] == PSystem->inner && Particles->particleType[i] == PSystem->fluid)
			pressure = (mi*PSystem->coeffPressWCMPS/PSystem->gamma)*(pow(Particles->pndi[i]/PSystem->pndSmallZero,PSystem->gamma)-1);
		}
	
//		if(Particles->pndSmall[i] < PSystem->pndThreshold*PSystem->pndSmallZero && numNeigh[i] < PSystem->neighThreshold*PSystem->numNeighZero)
//			Particles->particleBC[i] = PSystem->surface;
//		else
//		{
//			Particles->particleBC[i] = PSystem->inner;
//			pressure = (mi*PSystem->coeffPressWCMPS/PSystem->gamma)*(pow(Particles->pndi[i]/PSystem->pndSmallZero,PSystem->gamma)-1);
		//pressure = (ni - n0) * PSystem->coeffPressEMPS * mi;
		//pressure = (Particles->pndi[i] - PSystem->pndSmallZero) * PSystem->coeffPressEMPS * mi;
//		}
	
//			if(PSystem->wallType == boundaryWallType::POLYGON) {
//				if(Particles->pndi[i] > PSystem->pndThreshold*PSystem->pndSmallZero){
//					pressure = -mi*PSystem->gravityZ*(0.20 - 0.5*partDist -posZi);
//					pressure = -mi*PSystem->gravityZ*(0.18 - 0.5*partDist -posZi); // lat
//				}
//			}
//		else if(PSystem->wallType == boundaryWallType::PARTICLE) 
//				if(Particles->pndi[i] > PSystem->pndThreshold*PSystem->pndSmallZero || numNeigh[i] > PSystem->neighThreshold*PSystem->numNeighZero) {
//					pressure = -mi*PSystem->gravityZ*(0.20 - 0.5*partDist -posZi);
//					pressure = -mi*PSystem->gravityZ*(0.18 - 0.5*partDist -posZi); // lat
//				}
//			}
		
		if(pressure < 0.0) {
			pressure = 0.0;
		}
		Particles->press[i] = pressure;
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Solve linear system solver PPE (mpsType = calcPressType::IMPLICIT_PND)
// Perform conjugate gradient method on symmetry matrix A to solve Ax=b
// matA			symmetric (sparse) matrix
// sourceTerm	vector
void MpsPressure::solvePressurePoissonPnd(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {

	using T = Eigen::Triplet<double>;
	double lap_r = PSystem->reL*PSystem->invPartDist;
	int n_size = (int)(pow(lap_r * 2, PSystem->dim)); // maximum number of neighbors
	Eigen::SparseMatrix<double> matA(Particles->numParticles, Particles->numParticles); // declares a column-major sparse matrix type of double
	Particles->sourceTerm.resize(Particles->numParticles); // Resizing a dynamic-size matrix
	Particles->sourceTerm.setZero(); // Right hand side-vector set to zero
	vector<T> coeffs(Particles->numParticles * n_size); // list of non-zeros coefficients

// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleBC[i] == PSystem->other || Particles->particleBC[i] == PSystem->surface) {
			coeffs.push_back(T(i, i, 1.0));
			continue;
		}

		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double ni = Particles->pndSmall[i];
		double sum = 0.0;

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
					// PSystem->coeffPPE = 2.0*PSystem->dim/(PSystem->pndLargeZero*PSystem->lambdaZero)
					double mat_ij = wL*PSystem->coeffPPE;
					if (Particles->particleType[j] == PSystem->dummyWall) {
						double pgh = Particles->RHO[j]*(v0ij*PSystem->gravityX + v1ij*PSystem->gravityY + v2ij*PSystem->gravityZ)*wL;
						Particles->sourceTerm(i) -= PSystem->coeffPPE*pgh;
					}
					else
					{
						sum -= mat_ij;
						if (Particles->particleBC[j] == PSystem->inner) {
							coeffs.push_back(T(i, j, mat_ij));
						}
					}
				}}
				
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		double density;
		if(Particles->PTYPE[i] == 1) density = PSystem->DNS_FL1;
		else density = PSystem->DNS_FL2;

		// Increase diagonal
		sum -= PSystem->alphaCompressibility*density/(PSystem->timeStep*PSystem->timeStep);

		//double Cdiag = 1.0;
		//double beta = 0.9;
		//double DI2 = Cdiag*beta*PSystem->coeffPPE*PSystem->relaxPND*pndWallContribution[i]*density;
		//sum -= DI2;

		coeffs.push_back(T(i, i, sum));

		//PSystem->coeffPPESource = PSystem->relaxPND/(PSystem->timeStep*PSystem->timeStep*PSystem->pndSmallZero)
		Particles->sourceTerm(i) += - PSystem->coeffPPESource*density*(ni - PSystem->pndSmallZero);

		// 2019 - Enhancement of stabilization of MPS to arbitrary geometries with a generic wall boundary condition
		//double pndc = 0.0;
		//if(ni > 0)
		//	pndc = (ni - pndWallContribution[i])/ni;
		//Particles->sourceTerm(i) += pndc*((1.0-PSystem->relaxPND)*(Particles->pndki[i] - ni) + PSystem->relaxPND*(PSystem->pndSmallZero - pndski[i]))*ddt/PSystem->pndSmallZero;
		//Particles->sourceTerm(i) += - PSystem->relaxPND*ddt*(pndski[i] - PSystem->pndSmallZero)/PSystem->pndSmallZero;

		//double riw[3], riwSqrt;
		// Normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		//riw[0] = 0.5*(posXi - posMirrorXi); riw[1] = 0.5*(posYi - posMirrorYi); riw[2] = 0.5*(posZi - posMirrorZi);
		//riwSqrt = sqrt(riw[0]*riw[0] + riw[1]*riw[1] + riw[2]*riw[2]);
		//double ST1 = - PSystem->coeffPPESource*density*(ni - PSystem->pndSmallZero);
		//double ST2 = 0.0;
		//if(riwSqrt < 0.5*partDist)
		//	ST2 = - Cdiag*(1.0-beta)*2.0*partDist/PSystem->lambdaZero*(0.5*partDist - riwSqrt)/(PSystem->timeStep*PSystem->timeStep);
		//Particles->sourceTerm(i) += ST1 + ST2;
	}

	// Finished setup matrix
	matA.setFromTriplets(coeffs.begin(), coeffs.end());
	// Solve PPE
	if(PSystem->solverType == solvPressType::CG)
		solveConjugateGradient(matA, Particles);
	else if (PSystem->solverType == solvPressType::BICGSTAB)
		solveBiConjugateGradientStabilized(matA, Particles);
	// Set zero to negative pressures
	setZeroOnNegativePressure(Particles);

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Solve linear system solver PPE (mpsType = calcPressType::IMPLICIT_PND_DIVU)
void MpsPressure::solvePressurePoissonPndDivU(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {

	using T = Eigen::Triplet<double>;
	double lap_r = PSystem->reL*PSystem->invPartDist;
	int n_size = (int)(pow(lap_r * 2, PSystem->dim)); // maximum number of neighbors
	Eigen::SparseMatrix<double> matA(Particles->numParticles, Particles->numParticles); // declares a column-major sparse matrix type of double
	Particles->sourceTerm.resize(Particles->numParticles); // Resizing a dynamic-size matrix
	Particles->sourceTerm.setZero(); // Right hand side-vector set to zero
	vector<T> coeffs(Particles->numParticles * n_size); // list of non-zeros coefficients

// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleBC[i] == PSystem->other || Particles->particleBC[i] == PSystem->surface) {
			coeffs.push_back(T(i, i, 1.0));
			continue;
		}

		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		// double ni = Particles->pndSmall[i];
		double sum = 0.0;

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
					// PSystem->coeffPPE = 2.0*PSystem->dim/(PSystem->pndLargeZero*PSystem->lambdaZero)
					double mat_ij = wL*PSystem->coeffPPE;
					if (Particles->particleType[j] == PSystem->dummyWall) {
						double pgh = Particles->RHO[j]*(v0ij*PSystem->gravityX + v1ij*PSystem->gravityY + v2ij*PSystem->gravityZ)*wL;
						Particles->sourceTerm(i) -= PSystem->coeffPPE*pgh;
					}
					else
					{
						sum -= mat_ij;
						if (Particles->particleBC[j] == PSystem->inner) {
							coeffs.push_back(T(i, j, mat_ij));
						}
					}
				}}
				
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		double density;
		if(Particles->PTYPE[i] == 1) density = PSystem->DNS_FL1;
		else density = PSystem->DNS_FL2;

		// Increase diagonal
		sum -= PSystem->alphaCompressibility*density/(PSystem->timeStep*PSystem->timeStep);

		//double Cdiag = 1.0;
		//double beta = 0.9;
		//double DI2 = Cdiag*beta*PSystem->coeffPPE*PSystem->relaxPND*pndWallContribution[i]*density;
		//sum -= DI2;

		coeffs.push_back(T(i, i, sum));

		//PSystem->coeffPPESource = PSystem->relaxPND/(PSystem->timeStep*PSystem->timeStep*PSystem->pndSmallZero)
		//Particles->sourceTerm(i) += - PSystem->coeffPPESource*density*(ni - PSystem->pndSmallZero) + (1.0-PSystem->relaxPND)*density*Particles->velDivergence[i]/PSystem->timeStep;
		Particles->sourceTerm(i) += - PSystem->coeffPPESource*density*(Particles->pndki[i] - PSystem->pndSmallZero) + (1.0-PSystem->relaxPND)*density*Particles->velDivergence[i]/PSystem->timeStep;
		//Particles->sourceTerm(i) += - 4.0*density*(ni - PSystem->pndSmallZero)/(partDist*partDist*PSystem->pndSmallZero) 
		//				+ 2.0*density*Particles->velDivergence[i]*PSystem->invPartDist;

		// Sun et al., 2015. Modified MPS method for the 2D fluid structure interaction problem with free surface
		////double dtPhysical = partDist/20.0;
		//double dtPhysical = PSystem->timeStep;
		//double a1 = fabs(ni - PSystem->pndSmallZero)/PSystem->pndSmallZero;
		//if ((PSystem->pndSmallZero-ni)*Particles->velDivergence[i] > PSystem->epsilonZero)
		//{
		//	a1 += dtPhysical*fabs(Particles->velDivergence[i]);
		//}
		////double a2 = fabs((ni - PSystem->pndSmallZero)/PSystem->pndSmallZero);
		//Particles->sourceTerm(i) += - a1*density/(dtPhysical*dtPhysical)*(ni - PSystem->pndSmallZero)/PSystem->pndSmallZero 
		//	+ density*Particles->velDivergence[i]/dtPhysical;

		// 2019 - Enhancement of stabilization of MPS to arbitrary geometries with a generic wall boundary condition
		//double pndc = 0.0;
		//if(ni > 0)
		//	pndc = (ni - pndWallContribution[i])/ni;
		//Particles->sourceTerm(i) += pndc*((1.0-PSystem->relaxPND)*(Particles->pndki[i] - ni) + PSystem->relaxPND*(PSystem->pndSmallZero - pndski[i]))*ddt/PSystem->pndSmallZero;
		//Particles->sourceTerm(i) += - PSystem->relaxPND*ddt*(pndski[i] - PSystem->pndSmallZero)/PSystem->pndSmallZero;

		//double riw[3], riwSqrt;
		// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		//riw[0] = 0.5*(posXi - posMirrorXi); riw[1] = 0.5*(posYi - posMirrorYi); riw[2] = 0.5*(posZi - posMirrorZi);
		//riwSqrt = sqrt(riw[0]*riw[0] + riw[1]*riw[1] + riw[2]*riw[2]);
		//double ST1 = - PSystem->coeffPPESource*density*(ni - PSystem->pndSmallZero);
		//double ST2 = 0.0;
		//if(riwSqrt < 0.5*partDist)
		//	ST2 = - Cdiag*(1.0-beta)*2.0*partDist/PSystem->lambdaZero*(0.5*partDist - riwSqrt)/(PSystem->timeStep*PSystem->timeStep);
		//Particles->sourceTerm(i) += ST1 + ST2;
	}

	// Finished setup matrix
	matA.setFromTriplets(coeffs.begin(), coeffs.end());
	// Solve PPE
	if(PSystem->solverType == solvPressType::CG)
		solveConjugateGradient(matA, Particles);
	else if (PSystem->solverType == solvPressType::BICGSTAB)
		solveBiConjugateGradientStabilized(matA, Particles);
	// Set zero to negative pressures
	setZeroOnNegativePressure(Particles);

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Solve linear system using Conjugate Gradient (solverType = solvPressType::CG)
void MpsPressure::solveConjugateGradient(Eigen::SparseMatrix<double> p_mat, MpsParticle *Particles) {
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
	//cg.setTolerance(1.0e-9);
	cg.compute(p_mat);
	if (cg.info() != Eigen::ComputationInfo::Success) {
		cerr << "Error: Failed decompostion." << endl;
	}
	//Particles->pressurePPE = cg.solve(Particles->sourceTerm);
	Particles->pressurePPE = cg.solveWithGuess(Particles->sourceTerm, Particles->pressurePPE);

#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		Particles->press[i] = Particles->pressurePPE(i);
	}

	if (cg.info() != Eigen::ComputationInfo::Success) {
		cerr << "Error: Failed solving." << endl;
	}
	solverIter = cg.iterations();
	solverError = cg.error();
}

// Solve linear system using Bi Conjugate Gradient Stabilized (solverType = solvPressType::BICGSTAB)
void MpsPressure::solveBiConjugateGradientStabilized(Eigen::SparseMatrix<double> p_mat, MpsParticle *Particles) {
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> bicg;
	//bicg.setTolerance(1.0e-9);
	bicg.compute(p_mat);
	if (bicg.info() != Eigen::ComputationInfo::Success) {
		cerr << "Error: Failed decompostion." << endl;
	}

	Particles->pressurePPE = bicg.solveWithGuess(Particles->sourceTerm, Particles->pressurePPE);

#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		Particles->press[i] = Particles->pressurePPE(i);
	}

	if (bicg.info() != Eigen::ComputationInfo::Success) {
		cerr << "Error: Failed solving." << endl;
	}
	solverIter = bicg.iterations();
	solverError = bicg.error();
}

// Set negative pressure to zero
void MpsPressure::setZeroOnNegativePressure(MpsParticle *Particles){
#pragma omp parallel for
	for(int i=0; i<Particles->numParticles; i++) {
		if (Particles->press[i] < 0) Particles->press[i] = 0.0;
	}
}

// Divergence of velocity
void MpsPressure::calcVelDivergence(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		Particles->velDivergence[i] = 0.0;
		if(Particles->particleType[i] == PSystem->fluid) {
		double DivV = 0.0;
		double ni = Particles->pndi[i];
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
//						if(Particles->particleType[i] == PSystem->fluid && Particles->particleType[j] == PSystem->fluid) {
							double vijx = Particles->vel[j*3  ]-velXi;
							double vijy = Particles->vel[j*3+1]-velYi;
							double vijz = Particles->vel[j*3+2]-velZi;
							double dst = sqrt(dstij2);
							double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
							
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
//						}
					}
				}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		
		Particles->velDivergence[i] = DivV;
	}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Divergence of velocity (Polygon wall) - Free-slip
void MpsPressure::calcWallSlipVelDivergence(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
//	if(Particles->particleType[i] == PSystem->fluid) {
		double ni = Particles->pndi[i];
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true && ni > PSystem->epsilonZero) {
			double DivV = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			
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

			Particles->velDivergence[i] += DivV;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}


// Divergence of velocity (Polygon wall) - No-slip
void MpsPressure::calcWallNoSlipVelDivergence(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
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
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double velXi = Particles->vel[i*3  ];	double velYi = Particles->vel[i*3+1];	double velZi = Particles->vel[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];

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

			Particles->velDivergence[i] += DivV;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}


// Extrapolate pressure to wall and dummy particles
void MpsPressure::extrapolatePressParticlesWallDummy(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
		if(Particles->particleType[i] == PSystem->dummyWall || (PSystem->mpsType == calcPressType::WEAKLY && Particles->particleType[i] == PSystem->wall)) {
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double ni = 0.0;
			double pressure = 0.0;
			
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
	//				if(j != i) {
						if(j != i && Particles->particleType[j] == PSystem->fluid) {
							double dst = sqrt(dstij2);
							double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
							ni += wS;
							pressure += (Particles->press[j] - Particles->RHO[j]*(PSystem->gravityX*v0ij+PSystem->gravityY*v1ij+PSystem->gravityZ*v2ij))*wS;
							//pressure += (Particles->press[j] + Particles->RHO[j]*(PSystem->gravityX*v0ij+PSystem->gravityY*v1ij+PSystem->gravityZ*v2ij))*wS;
						}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			if(pressure < 0.0){
				pressure = 0.0;
				continue;
			}
			if(ni > 0)
				Particles->press[i] = pressure/ni;
			else
				Particles->press[i] = pressure;
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Extrapolate pressure to inner particles near polygon walls
void MpsPressure::extrapolatePressParticlesNearPolygonWall(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
	//if(Particles->particleType[i] == PSystem->fluid && Particles->press[i] == 0 /*&& Particles->particleNearWall[i] == true*/) {
		if(Particles->particleType[i] == PSystem->fluid && Particles->press[i] == 0 && Particles->particleNearWall[i] == true) {
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			double pressure = 0.0;
			double sumWij = 0.0;
			int nTotal = 0;
			int nFree = 0;
			
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
					if(dstij2 < PSystem->reS2 && dstij2 < dstimj2) {
	//				if(dstij2 < 1.2*partDist) {
						if(j != i) {
							double dst = sqrt(dstij2);
							double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
		//					double wS = WEI_WEND(dst, 1.2*partDist);
							sumWij += wS;
		//					pressure += Particles->press[j]*wS;
							pressure += (Particles->press[j] - Particles->RHO[j]*(PSystem->gravityX*v0ij+PSystem->gravityY*v1ij+PSystem->gravityZ*v2ij))*wS;
							nTotal += 1;
							if(Particles->particleBC[j] == PSystem->surface)
								nFree += 1;
						}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			
			if(nTotal > 0)
				Particles->numNeighborsSurfaceParticles[i] = double(nFree)/nTotal;
			else
				Particles->numNeighborsSurfaceParticles[i] = 1.0;

	//		if(Particles->numNeighborsSurfaceParticles[i]<=0.5){
			if(pressure > 0) {
		  		if(sumWij > PSystem->epsilonZero)
					Particles->press[i] = pressure/sumWij;
				else
					Particles->press[i] = pressure;
			}
	//		}
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}




// Prediction of pressure gradient
void MpsPressure::predictionPressGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++){
//		if(Particles->particleType[i] == PSystem->fluid) {
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double Pi = Particles->press[i];			double ni = Particles->pndi[i];		double pressMin = Pi;
		
		int ix, iy, iz;
		Buckets->bucketCoordinates(ix, iy, iz, posXi, posYi, posZi, PSystem);
		int minZ = (iz-1)*((int)(PSystem->dim-2.0)); int maxZ = (iz+1)*((int)(PSystem->dim-2.0));
		if(PSystem->gradientType == 0 || PSystem->gradientType == 2) {
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
							if(pressMin > Particles->press[j]) pressMin = Particles->press[j];
						}
					}
					j = Particles->nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
		}
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
						double wS = Particles->weightGradient(dst, PSystem->reS, PSystem->weightType);
						if(PSystem->gradientType == 0)
							wS *= (Particles->press[j] - pressMin)/dstij2;
						else if(PSystem->gradientType == 1)
							wS *= (Particles->press[j] + Pi)/dstij2;
						else if(PSystem->gradientType == 2)
							wS *= (Particles->press[j] + Pi - 2.0*pressMin)/dstij2;
						else if(PSystem->gradientType == 3) {
							double nj = Particles->pndi[j];
							if(ni > PSystem->epsilonZero && nj > PSystem->epsilonZero)
								wS *= (ni*Particles->press[j]/nj + nj*Pi/ni)/dstij2;
						}
						accX += v0ij*wS;	accY += v1ij*wS;	accZ += v2ij*wS;
					}
				}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/noGrad)
		// Original
//		Particles->acc[i*3  ]+=(1.0-PSystem->relaxPress)*accX*invDns[partType::FLUID]*PSystem->coeffPressGrad;
//		Particles->acc[i*3+1]+=(1.0-PSystem->relaxPress)*accY*invDns[partType::FLUID]*PSystem->coeffPressGrad;
//		Particles->acc[i*3+2]+=(1.0-PSystem->relaxPress)*accZ*invDns[partType::FLUID]*PSystem->coeffPressGrad;
		// Modified
		Particles->acc[i*3  ]+=(1.0-PSystem->relaxPress)*accX*PSystem->coeffPressGrad/Particles->RHO[i];
		Particles->acc[i*3+1]+=(1.0-PSystem->relaxPress)*accY*PSystem->coeffPressGrad/Particles->RHO[i];
		Particles->acc[i*3+2]+=(1.0-PSystem->relaxPress)*accZ*PSystem->coeffPressGrad/Particles->RHO[i];
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Prediction of pressure gradient (Polygon wall)
void MpsPressure::predictionWallPressGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	// Maximum velocity is the minimum of the computed and expected maximum velocities
	double maxVelocity = min(PSystem->velMax, PSystem->expectMaxVelocity);
	double velMax2 = maxVelocity*maxVelocity;
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {
		//Particles->particleNearWall[i]=true; // Only to show particles near polygon
		if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
			double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
			double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
			double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
			double Pi = Particles->press[i];			double ni = Particles->pndi[i];		double pressMin = Pi;
			
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
		    	normaliw[0] = 0.0;
		    	normaliw[1] = 0.0;
		    	normaliw[2] = 0.0;
		    }
		    //  Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
		    Rref_i[0] = 1.0 - 2.0*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2.0*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2.0*normaliw[0]*normaliw[2];
			Rref_i[3] = 0.0 - 2.0*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2.0*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2.0*normaliw[1]*normaliw[2];
			Rref_i[6] = 0.0 - 2.0*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2.0*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2.0*normaliw[2]*normaliw[2];
			// Taylor pressure Pj
			double Rai[3];
			Rai[0] = Rref_i[0]*Particles->accStar[i*3] + Rref_i[1]*Particles->accStar[i*3+1] + Rref_i[2]*Particles->accStar[i*3+2];
			Rai[1] = Rref_i[3]*Particles->accStar[i*3] + Rref_i[4]*Particles->accStar[i*3+1] + Rref_i[5]*Particles->accStar[i*3+2];
			Rai[2] = Rref_i[6]*Particles->accStar[i*3] + Rref_i[7]*Particles->accStar[i*3+1] + Rref_i[8]*Particles->accStar[i*3+2];

			int ix, iy, iz;
			Buckets->bucketCoordinates(ix, iy, iz, posXi, posYi, posZi, PSystem);
			int minZ = (iz-1)*((int)(PSystem->dim-2.0)); int maxZ = (iz+1)*((int)(PSystem->dim-2.0));
			if(PSystem->gradientType == 0 || PSystem->gradientType == 2) {
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
						if(dstij2 < PSystem->reS2 && dstij2 < dstimj2) {
						if(j != i) {
							if(pressMin > Particles->press[j]) pressMin = Particles->press[j];
						}}
						j = Particles->nextParticleInSameBucket[j];
						if(j == -1) break;
					}
				}}}
			}
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
							
							// Taylor pressure Pj
							double Pj;
							// Pj = Pi + Particles->RHO[i]*(Rai[0]*v0 + Rai[1]*v1 + Rai[2]*v2);
							Pj = Particles->press[j];

							if(PSystem->gradientType == 0)
								wS *= (Pj - pressMin)/dstimj2;//(Particles->press[j] - pressMin)/dstimj2;
							else if(PSystem->gradientType == 1)
								wS *= (Pj + Pi)/dstimj2;//(Particles->press[j] + Pi)/dstimj2;
							else if(PSystem->gradientType == 2)
								wS *= (Pj + Pi - 2.0*pressMin)/dstimj2;//(Particles->press[j] + Pi - 2.0*pressMin)/dstimj2;
							else if(PSystem->gradientType == 3) {
								double nj = Particles->pndi[j];
								if(ni > PSystem->epsilonZero && nj > PSystem->epsilonZero)
									wS *= (ni*Pj/nj + nj*Pi/ni)/dstimj2;
							}
						accX += v0imj*wS;	accY += v1imj*wS;	accZ += v2imj*wS;
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

				// Taylor pressure Pj
				double Pj;
				// Pj = Pi + Particles->RHO[i]*(Rai[0]*v0 + Rai[1]*v1 + Rai[2]*v2);
				Pj = Pi;

				if(PSystem->gradientType == 0)
					wS *= (Pj - pressMin)/dstimi2;//(Pi - pressMin)/dstimi2
				else if(PSystem->gradientType == 1)
					wS *= (Pj + Pi)/dstimi2;//(Pi + Pi)/dstimi2
				else if(PSystem->gradientType == 2)
					wS *= (Pj + Pi - 2.0*pressMin)/dstimi2;//(Pi + Pi - 2.0*pressMin)/dstimi2
				else if(PSystem->gradientType == 3) {
					double nj = Particles->pndi[i];
					if(ni > PSystem->epsilonZero && nj > PSystem->epsilonZero)
						wS *= (ni*Pj/nj + nj*Pi/ni)/dstimi2;
				}
				accX += v0imi*wS;	accY += v1imi*wS;	accZ += v2imi*wS;
		  	}

			// Repulsive force
			double rpsForce[3];
			rpsForce[0]=rpsForce[1]=rpsForce[2] = 0.0;

			if(normaliwSqrt < PSystem->reRepulsiveForce && normaliwSqrt > PSystem->epsilonZero) {
				if(PSystem->repulsiveForceType == repForceType::HARADA)
					repulsiveForceHarada(rpsForce, normaliw, normaliwSqrt, i, PSystem, Particles);
				else if(PSystem->repulsiveForceType == repForceType::MITSUME)
					repulsiveForceMitsume(rpsForce, normaliw, normaliwSqrt, i, PSystem, Particles);
				else if(PSystem->repulsiveForceType == repForceType::LENNARD_JONES)
					repulsiveForceLennardJones(rpsForce, normaliw, normaliwSqrt, velMax2, i, PSystem, Particles);
				else
					repulsiveForceMonaghanKajtar(rpsForce, normaliw, normaliwSqrt, velMax2, i, PSystem, Particles);
			}
			
			// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/noGrad)
			// Original
	//		Particles->acc[i*3  ] += ((1.0-PSystem->relaxPress)*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*PSystem->coeffPressGrad - rpsForce[0])*invDns[partType::FLUID];
	//		Particles->acc[i*3+1] += ((1.0-PSystem->relaxPress)*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*PSystem->coeffPressGrad - rpsForce[1])*invDns[partType::FLUID];
	//		Particles->acc[i*3+2] += ((1.0-PSystem->relaxPress)*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*PSystem->coeffPressGrad - rpsForce[2])*invDns[partType::FLUID];
			// Modified
			Particles->acc[i*3  ] += ((1.0-PSystem->relaxPress)*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*PSystem->coeffPressGrad - rpsForce[0])/Particles->RHO[i];
			Particles->acc[i*3+1] += ((1.0-PSystem->relaxPress)*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*PSystem->coeffPressGrad - rpsForce[1])/Particles->RHO[i];
			Particles->acc[i*3+2] += ((1.0-PSystem->relaxPress)*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*PSystem->coeffPressGrad - rpsForce[2])/Particles->RHO[i];

			//Fwall[i*3  ] =  (Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*invDns[partType::FLUID]*PSystem->coeffPressGrad - rpsForce[0]*invDns[partType::FLUID];
			//Fwall[i*3+1] =  (Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*invDns[partType::FLUID]*PSystem->coeffPressGrad - rpsForce[1]*invDns[partType::FLUID];
			//Fwall[i*3+2] =  (Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*invDns[partType::FLUID]*PSystem->coeffPressGrad - rpsForce[2]*invDns[partType::FLUID];
		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Acceleration due to pressure gradient
void MpsPressure::calcPressGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
//	if(Particles->particleType[i] == PSystem->fluid) {
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double pressMin = Particles->press[i];
		double Pi = Particles->press[i];
		double ni = Particles->pndi[i];
		double MC[9];
		MC[0] = Particles->correcMatrixRow1[i*3];	MC[1] = Particles->correcMatrixRow1[i*3+1];	MC[2] = Particles->correcMatrixRow1[i*3+2];
		MC[3] = Particles->correcMatrixRow2[i*3];	MC[4] = Particles->correcMatrixRow2[i*3+1];	MC[5] = Particles->correcMatrixRow2[i*3+2];
		MC[6] = Particles->correcMatrixRow3[i*3];	MC[7] = Particles->correcMatrixRow3[i*3+1];	MC[8] = Particles->correcMatrixRow3[i*3+2];
		int ix, iy, iz;
		Buckets->bucketCoordinates(ix, iy, iz, posXi, posYi, posZi, PSystem);
		int minZ = (iz-1)*((int)(PSystem->dim-2.0)); int maxZ = (iz+1)*((int)(PSystem->dim-2.0));
		if(PSystem->gradientType == 0 || PSystem->gradientType == 2) {
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
					if(pressMin > Particles->press[j]) pressMin = Particles->press[j];
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}}
		
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
					double wS = Particles->weightGradient(dst, PSystem->reS, PSystem->weightType);
					if(PSystem->gradientType == 0)
						wS *= (Particles->press[j] - pressMin)/dstij2;
					else if(PSystem->gradientType == 1)
						wS *= (Particles->press[j] + Pi)/dstij2;
					else if(PSystem->gradientType == 2)
						wS *= (Particles->press[j] + Pi - 2.0*pressMin)/dstij2;
					else if(PSystem->gradientType == 3) {
						double nj = Particles->pndi[j];
						if(ni > PSystem->epsilonZero && nj > PSystem->epsilonZero)
							wS *= (ni*Particles->press[j]/nj + nj*Pi/ni)/dstij2;
					}
					if(PSystem->gradientCorrection == false) {
						accX += v0ij*wS;	accY += v1ij*wS;	accZ += v2ij*wS;
					}
					else {
						accX += (v0ij*MC[0] + v1ij*MC[1] + v2ij*MC[2])*wS;
						accY += (v0ij*MC[3] + v1ij*MC[4] + v2ij*MC[5])*wS;
						accZ += (v0ij*MC[6] + v1ij*MC[7] + v2ij*MC[8])*wS;
					}
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/noGrad)
		// Original
//		Particles->acc[i*3  ]=PSystem->relaxPress*accX*invDns[partType::FLUID]*PSystem->coeffPressGrad;
//		Particles->acc[i*3+1]=PSystem->relaxPress*accY*invDns[partType::FLUID]*PSystem->coeffPressGrad;
//		Particles->acc[i*3+2]=PSystem->relaxPress*accZ*invDns[partType::FLUID]*PSystem->coeffPressGrad;
		// Modified
		Particles->acc[i*3  ]=PSystem->relaxPress*accX*PSystem->coeffPressGrad/Particles->RHO[i];
		Particles->acc[i*3+1]=PSystem->relaxPress*accY*PSystem->coeffPressGrad/Particles->RHO[i];
		Particles->acc[i*3+2]=PSystem->relaxPress*accZ*PSystem->coeffPressGrad/Particles->RHO[i];
		/*if(PSystem->gradientCorrection == false) {
			Particles->acc[i*3  ]=PSystem->relaxPress*accX*PSystem->coeffPressGrad/Particles->RHO[i];
			Particles->acc[i*3+1]=PSystem->relaxPress*accY*PSystem->coeffPressGrad/Particles->RHO[i];
			Particles->acc[i*3+2]=PSystem->relaxPress*accZ*PSystem->coeffPressGrad/Particles->RHO[i];
		}
		else {
		//	if(Particles->correcMatrixRow1[1*3] > 1.0) {
		//		printf("\n X %e %e %e ", Particles->correcMatrixRow1[i*3  ], Particles->correcMatrixRow1[i*3+1], Particles->correcMatrixRow1[i*3+2]);
		//		printf("\n Y %e %e %e ", Particles->correcMatrixRow2[i*3  ], Particles->correcMatrixRow2[i*3+1], Particles->correcMatrixRow2[i*3+2]);
		//		printf("\n Z %e %e %e \n", Particles->correcMatrixRow3[i*3  ], Particles->correcMatrixRow3[i*3+1], Particles->correcMatrixRow3[i*3+2]);
			//}
			Particles->acc[i*3  ]=(PSystem->relaxPress*PSystem->coeffPressGrad/Particles->RHO[i])*(accX*Particles->correcMatrixRow1[i*3] + accY*Particles->correcMatrixRow1[i*3+1] + accZ*Particles->correcMatrixRow1[i*3+2]);
			Particles->acc[i*3+1]=(PSystem->relaxPress*PSystem->coeffPressGrad/Particles->RHO[i])*(accX*Particles->correcMatrixRow2[i*3] + accY*Particles->correcMatrixRow2[i*3+1] + accZ*Particles->correcMatrixRow2[i*3+2]);
			Particles->acc[i*3+2]=(PSystem->relaxPress*PSystem->coeffPressGrad/Particles->RHO[i])*(accX*Particles->correcMatrixRow3[i*3] + accY*Particles->correcMatrixRow3[i*3+1] + accZ*Particles->correcMatrixRow3[i*3+2]);
		}*/
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Acceleration due to pressure gradient (Polygon wall)
void MpsPressure::calcWallPressGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
	//int nPartNearMesh = partNearMesh.size();
	double VolumeForce = pow(PSystem->partDist,PSystem->dim);
	// Maximum velocity is the minimum of the computed and expected maximum velocities
	double maxVelocity = min(PSystem->velMax, PSystem->expectMaxVelocity);
	double velMax2 = maxVelocity*maxVelocity;
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<Particles->numParticles; i++) {

	//Particles->particleNearWall[i]=true; // Only to show particles near polygon
	
	//if(Particles->particleType[i] == PSystem->fluid) {
	if(Particles->particleType[i] == PSystem->fluid && Particles->particleNearWall[i] == true) {
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		double pressMin = Particles->press[i];
		double Pi = Particles->press[i];
		double ni = Particles->pndi[i];
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

		// Taylor pressure Pj
		double Rai[3];
		Rai[0] = Rref_i[0]*Particles->accStar[i*3] + Rref_i[1]*Particles->accStar[i*3+1] + Rref_i[2]*Particles->accStar[i*3+2];
		Rai[1] = Rref_i[3]*Particles->accStar[i*3] + Rref_i[4]*Particles->accStar[i*3+1] + Rref_i[5]*Particles->accStar[i*3+2];
		Rai[2] = Rref_i[6]*Particles->accStar[i*3] + Rref_i[7]*Particles->accStar[i*3+1] + Rref_i[8]*Particles->accStar[i*3+2];

		// if(i == 16107)
		// {
		// 	printf("\ni:%5d timeCurrent: %lf / Rref_i: ", i, PSystem->timeCurrent);
		// 	for(int rr=0; rr<9; rr++)
		// 		printf("%lf ", Rref_i[rr]);
		// 	printf("\ni:%5d timeCurrent: %lf / Rai: ", i, PSystem->timeCurrent);
		// 	for(int rr=0; rr<3; rr++)
		// 		printf("%lf ", Rai[rr]);
		// 	printf("\ni:%5d timeCurrent: %lf / ai: %lf, %lf, %lf", i, PSystem->timeCurrent, Particles->accStar[i*3], Particles->accStar[i*3+1], Particles->accStar[i*3+2]);
			
		// }
			
		int ix, iy, iz;
		Buckets->bucketCoordinates(ix, iy, iz, posXi, posYi, posZi, PSystem);
		int minZ = (iz-1)*((int)(PSystem->dim-2.0)); int maxZ = (iz+1)*((int)(PSystem->dim-2.0));
		if(PSystem->gradientType == 0 || PSystem->gradientType == 2) {
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
				if(dstij2 < PSystem->reS2 && dstij2 < dstimj2) {
				if(j != i) {
					if(pressMin > Particles->press[j]) pressMin = Particles->press[j];
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}}
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

					// Taylor pressure Pj
					double Pj;
					// Pj = Pi + Particles->RHO[i]*(Rai[0]*v0imj + Rai[1]*v1imj + Rai[2]*v2imj);
					Pj = Particles->press[j];

					// if(i == 16107)
					// 	printf("\ni:%5d j:%5d timeCurrent: %lf / Pj: %lf / Pi: %lf / Zj: %lf / Zi: %lf", i, j, PSystem->timeCurrent, Pj, Pi, Particles->pos[j*3+2], posMirrorZi);

					if(PSystem->gradientType == 0)
						wS *= (Pj - pressMin)/dstimj2;//(Particles->press[j] - pressMin)/dstimj2
					else if(PSystem->gradientType == 1)
						wS *= (Pj + Pi)/dstimj2;//(Particles->press[j] + Particles->press[i])/dstimj2
					else if(PSystem->gradientType == 2)
						wS *= (Pj + Pi - 2.0*pressMin)/dstimj2;//(Particles->press[j] + Particles->press[i] - 2.0*pressMin)/dstimj2
					else if(PSystem->gradientType == 3) {
						double nj = Particles->pndi[j];
						if(ni > PSystem->epsilonZero && nj > PSystem->epsilonZero)
							wS *= (ni*Pj/nj + nj*Pi/ni)/dstimj2;
					}
					accX += v0imj*wS;	accY += v1imj*wS;	accZ += v2imj*wS;
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

			// Taylor pressure Pj
			double Pj;
			// Pj = Pi + Particles->RHO[i]*(Rai[0]*v0imi + Rai[1]*v1imi + Rai[2]*v2imi);
			Pj = Pi;

			//if(i == 16107)
			//	printf("\ni:%5d timeCurrent: %lf / Pj: %lf / Pi: %lf / Zj: %lf / Zi: %lf", i, PSystem->timeCurrent, Pj, Pi, posZi, posMirrorZi);

			if(PSystem->gradientType == 0)
				wS *= (Pj - pressMin)/dstimi2;//(Particles->press[i] - pressMin)/dstimi2
			else if(PSystem->gradientType == 1)
				wS *= (Pj + Pi)/dstimi2;//(Particles->press[i] + Particles->press[i])/dstimi2;
			else if(PSystem->gradientType == 2)
				wS *= (Pj + Pi - 2.0*pressMin)/dstimi2;//(Particles->press[i] + Particles->press[i] - 2.0*pressMin)/dstimi2;
			else if(PSystem->gradientType == 3) {
				double nj = Particles->pndi[i];
				if(ni > PSystem->epsilonZero && nj > PSystem->epsilonZero)
					wS *= (ni*Pj/nj + nj*Pi/ni)/dstimi2;
			}
			accX += v0imi*wS;	accY += v1imi*wS;	accZ += v2imi*wS;
	  	}

		// Repulsive force
		double rpsForce[3];
		rpsForce[0]=rpsForce[1]=rpsForce[2] = 0.0;

		if(normaliwSqrt < PSystem->reRepulsiveForce && normaliwSqrt > PSystem->epsilonZero) {
			if(PSystem->repulsiveForceType == repForceType::HARADA)
				repulsiveForceHarada(rpsForce, normaliw, normaliwSqrt, i, PSystem, Particles);
			else if(PSystem->repulsiveForceType == repForceType::MITSUME)
				repulsiveForceMitsume(rpsForce, normaliw, normaliwSqrt, i, PSystem, Particles);
			else if(PSystem->repulsiveForceType == repForceType::LENNARD_JONES)
				repulsiveForceLennardJones(rpsForce, normaliw, normaliwSqrt, velMax2, i, PSystem, Particles);
			else
				repulsiveForceMonaghanKajtar(rpsForce, normaliw, normaliwSqrt, velMax2, i, PSystem, Particles);
		}

		// PSystem->coeffPressGrad is a negative cte (-PSystem->dim/noGrad)
		// Original
//		Particles->acc[i*3  ] += (PSystem->relaxPress*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*PSystem->coeffPressGrad - rpsForce[0])*invDns[partType::FLUID];
//		Particles->acc[i*3+1] += (PSystem->relaxPress*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*PSystem->coeffPressGrad - rpsForce[1])*invDns[partType::FLUID];
//		Particles->acc[i*3+2] += (PSystem->relaxPress*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*PSystem->coeffPressGrad - rpsForce[2])*invDns[partType::FLUID];
		// Modified
		Particles->acc[i*3  ] += (PSystem->relaxPress*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*PSystem->coeffPressGrad - rpsForce[0])/Particles->RHO[i];
		Particles->acc[i*3+1] += (PSystem->relaxPress*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*PSystem->coeffPressGrad - rpsForce[1])/Particles->RHO[i];
		Particles->acc[i*3+2] += (PSystem->relaxPress*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*PSystem->coeffPressGrad - rpsForce[2])/Particles->RHO[i];

		// FSI
		// Force on wall
		Particles->forceWall[i*3  ] += - (PSystem->relaxPress*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*PSystem->coeffPressGrad - rpsForce[0])*VolumeForce;
		Particles->forceWall[i*3+1] += - (PSystem->relaxPress*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*PSystem->coeffPressGrad - rpsForce[1])*VolumeForce;
		Particles->forceWall[i*3+2] += - (PSystem->relaxPress*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*PSystem->coeffPressGrad - rpsForce[2])*VolumeForce;
	}}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Parallel analysis system for free-surface flow using MPS method with explicitly represented polygon wall boundary model
// https://doi.org/10.1007/s40571-019-00269-6
void MpsPressure::repulsiveForceHarada(double *force, const double *normal, const double normalSqrt, const int i, MpsParticleSystem *PSystem, MpsParticle *Particles) {
	double wijRep = Particles->RHO[i]/(PSystem->timeStep*PSystem->timeStep)*(PSystem->reRepulsiveForce-normalSqrt);
	force[0] = - wijRep*normal[0];
	force[1] = - wijRep*normal[1];
	force[2] = - wijRep*normal[2];

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Parallel analysis system for free-surface flow using MPS method with explicitly represented polygon wall boundary model
// https://doi.org/10.1007/s40571-019-00269-6
void MpsPressure::repulsiveForceMitsume(double *force, const double *normal, const double normalSqrt, const int i, MpsParticleSystem *PSystem, MpsParticle *Particles) {
	double wijRep = PSystem->repForceCoefMitsume*Particles->weightGradient(normalSqrt, PSystem->reRepulsiveForce, PSystem->weightType);
	force[0] = - wijRep*normal[0];
	force[1] = - wijRep*normal[1];
	force[2] = - wijRep*normal[2];

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Simulating Free Surface Flows with SPH
// https://doi.org/10.1006/jcph.1994.1034
void MpsPressure::repulsiveForceLennardJones(double *force, const double *normal, const double normalSqrt, const double vMax2, const int i,MpsParticleSystem *PSystem,  MpsParticle *Particles) {
	
	double R1 = (PSystem->reRepulsiveForce/normalSqrt)*(PSystem->reRepulsiveForce/normalSqrt);
	double R2 = R1*R1;
	double wijRep = (PSystem->repForceCoefLennardJones*vMax2/normalSqrt)*(R2-R1)*Particles->RHO[i];
	force[0] = - wijRep*normal[0];
	force[1] = - wijRep*normal[1];
	force[2] = - wijRep*normal[2];

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// SPH particle boundary forces for arbitrary boundaries 
// https://doi.org/10.1016/j.cpc.2009.05.008
void MpsPressure::repulsiveForceMonaghanKajtar(double *force, const double *normal, const double normalSqrt, const double vMax2, const int i, MpsParticleSystem *PSystem, MpsParticle *Particles) {
	double W1 = (1.0+3.0*0.5*normalSqrt/(PSystem->reRepulsiveForce));
	double W2 = (1.0-normalSqrt/(PSystem->reRepulsiveForce))*(1.0-normalSqrt/(PSystem->reRepulsiveForce))*(1.0-normalSqrt/(PSystem->reRepulsiveForce));
	double wijRep = (PSystem->repForceCoefMonaghanKajtar*vMax2/(normalSqrt - 0.0*PSystem->partDist))*(1.0/8.0)*(W1)*(W2)*Particles->RHO[i];
	force[0] = - wijRep*normal[0];
	force[1] = - wijRep*normal[1];
	force[2] = - wijRep*normal[2];

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

