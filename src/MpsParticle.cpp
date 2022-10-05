// Copyright (c) 2021 Rubens AMARO
// Distributed under the MIT License.

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "MpsParticle.h"

using namespace std;

// Constructor declaration
MpsParticle::MpsParticle()
{
}
// Destructor declaration
MpsParticle::~MpsParticle()
{
}

// Return the square distance between thwo particles "i" and "j"
void MpsParticle::sqrDistBetweenParticles(const int j, 
	const double rxi, const double ryi, const double rzi,
	double &rx, double &ry, double &rz, double &rij2) {
	rx = pos[j*3  ] - rxi;
	ry = pos[j*3+1] - ryi;
	rz = pos[j*3+2] - rzi;

	rij2 = rx*rx+ry*ry+rz*rz;
}
// Return the square distance between thwo particles "i" and "j" considering Periodic BC
void MpsParticle::sqrDistBetweenParticles(const int j, 
	const double rxi, const double ryi, const double rzi,
	double &rx, double &ry, double &rz, double &rij2, 
	const double plx, const double ply, const double plz) {

	rx = (pos[j*3  ] + plx) - rxi;
	ry = (pos[j*3+1] + ply) - ryi;
	rz = (pos[j*3+2] + plz) - rzi;

	rij2 = rx*rx+ry*ry+rz*rz;
}

// Get Periodic lenghts
void MpsParticle::getPeriodicLengths(const int jb, double &perlx, double &perly, 
	double &perlz, MpsParticleSystem *PSystem) {
	int bPBC = bucketPeriodicBC[jb];
	perlx = bPBC*PSystem->periodicLength[0];
	perly = bPBC*PSystem->periodicLength[1];
	perlz = bPBC*PSystem->periodicLength[2];
}

// Weight function
double MpsParticle::weight(const double dst, const double re, const int wijType) {
	switch (wijType) {
		case 0:
			return re/dst - 1.0;
		case 1:
			return re/dst + dst/re - 2.0;
		case 2:
			return re/dst - dst/re;
		case 3:
			return (1.0-dst/re)*(1.0-dst/re)*(1.0-dst/re);
		case 4:
			return (1.0-dst/re)*(1.0-dst/re);
		default:
			return re/dst - 1.0;
	}
}

// Weight function for gradient
double MpsParticle::weightGradient(const double dst, const double re, const int wijType) {
	switch (wijType) {
		case 0:
			return re/dst - 1.0;
		case 1:
			return re/dst + dst/re - 2.0;
		case 2:
			return re/dst - dst/re;
		case 3:
			return (1.0-dst/re)*(1.0-dst/re)*(1.0-dst/re);
		case 4:
			return (1.0-dst/re)*(1.0-dst/re);
		default:
			return re/dst - 1.0;
	}
}

// Derivate of weight function
double MpsParticle::delWeight(const double dst, const double re, const int wijType) {
	switch (wijType) {
		case 0:
			return -re/(dst*dst);
		case 1:
			return -re/(dst*dst) + 1.0/re;
		case 2:
			return -re/(dst*dst) - 1.0/re;
		case 3:
			return -3.0/re*(1.0-dst/re)*(1.0-dst/re);
		case 4:
			return -2.0/re*(1.0-dst/re);
		default:
			return -re/(dst*dst);
	}
}

////////////////////////////////////////////////////////////
// Functions called only at the initial instant (t=0)
////////////////////////////////////////////////////////////
// Set parameters
void MpsParticle::setParameters(MpsParticleSystem *PSystem) {

	// Extend domain
	PSystem->domainMinX -= PSystem->partDist*3.0;
	PSystem->domainMinY -= PSystem->partDist*3.0;
	PSystem->domainMinZ -= PSystem->partDist*3.0;
	PSystem->domainMaxX += PSystem->partDist*3.0;
	PSystem->domainMaxY += PSystem->partDist*3.0;
	PSystem->domainMaxZ += PSystem->partDist*3.0;
	if((int)PSystem->dim == 2) {	
		PSystem->domainMinZ = 0.0;
		PSystem->domainMaxZ = 0.0;
	}

	// Number of meshs
	PSystem->numOfRigidMesh = 0;	PSystem->numOfDeformableMesh = 0;	PSystem->numOfForcedMesh = 0;
	if(PSystem->wallType == 1) PSystem->numOfRigidMesh = 1;
	if(PSystem->femOn == true) PSystem->numOfDeformableMesh = 1;
	if(PSystem->forcedOn == true) PSystem->numOfForcedMesh = 1;
	PSystem->numOfMeshs = PSystem->numOfRigidMesh + PSystem->numOfDeformableMesh + PSystem->numOfForcedMesh;

	// Normalized input variables are multiplied by the particle distance
	PSystem->reS = PSystem->partDist*PSystem->reS;								///< Influence radius small
	PSystem->reL = PSystem->partDist*PSystem->reL;								///< Influence radius large
	PSystem->reRepulsiveForce = PSystem->partDist*PSystem->reRepulsiveForce;	///< Influence radius for repulsive force

	PSystem->reS2 = PSystem->reS*PSystem->reS;									///< Influence radius small to square
	PSystem->reL2 = PSystem->reL*PSystem->reL;									///< Influence radius large to square
	PSystem->EPS_RE = PSystem->EPS_RE*PSystem->reS2/4.0;

	// Computes Initial particle number density and Initial number of neighbors for a fully filled compact support
	PSystem->pndSmallZero = PSystem->pndLargeZero = PSystem->pndGradientZero = PSystem->lambdaZero = PSystem->numNeighZero = 0.0;
	int lmin = ceil(PSystem->reL/PSystem->partDist) + 1;
	int lmax = ceil(PSystem->reL/PSystem->partDist) + 2;
	int flag2D = 0;
	int flag3D = 1;
	if((int)PSystem->dim == 2) {
		flag2D = 1;
		flag3D = 0;
	}
	for(int ix= -lmin; ix<lmax; ix++) {
	for(int iy= -lmin; iy<lmax; iy++) {
	for(int iz= -lmin*flag3D; iz<lmax*flag3D+flag2D; iz++) {
		double x = PSystem->partDist* (double)ix;
		double y = PSystem->partDist* (double)iy;
		double z = PSystem->partDist* (double)iz;
		double dst2 = x*x+y*y+z*z;
		if(dst2 <= PSystem->reL2) {
			if(dst2 <= PSystem->epsilonZero) continue; 						// equals to zero
			double dst = sqrt(dst2);
			PSystem->pndLargeZero += weight(dst, PSystem->reL, PSystem->weightType);	///< Initial particle number density (large)
			PSystem->lambdaZero += dst2 * weight(dst, PSystem->reL, PSystem->weightType);
			PSystem->numNeighZero += 1;										///< Initial number of neighbors
			if(dst2 <= PSystem->reS2) {
				PSystem->pndSmallZero += weight(dst, PSystem->reS, PSystem->weightType);			///< Initial particle number density (small)
				PSystem->pndGradientZero += weightGradient(dst, PSystem->reS, PSystem->weightType);	///< Initial particle number density (gradient operator)
			}
		}
	}}}
	PSystem->lambdaZero = PSystem->lambdaZero/PSystem->pndLargeZero;							///< Coefficient λ of Laplacian model
	PSystem->coeffViscosity = 2.0*PSystem->KNM_VS1*PSystem->dim/(PSystem->pndLargeZero*PSystem->lambdaZero);	///< Coefficient used to calculate viscosity term
	PSystem->coeffViscMultiphase = 2.0*PSystem->dim/(PSystem->pndLargeZero*PSystem->lambdaZero);///< Coefficient used to calculate viscosity term Multiphase
	PSystem->coeffPressEMPS = PSystem->soundSpeed*PSystem->soundSpeed/PSystem->pndSmallZero;	///< Coefficient used to calculate pressure E-MPS
	PSystem->coeffPressGrad = -PSystem->dim/PSystem->pndGradientZero;							///< Coefficient used to calculate pressure gradient term
	PSystem->coeffPressWCMPS = PSystem->soundSpeed*PSystem->soundSpeed;							///< Coefficient used to calculate pressure WC-MPS
	PSystem->coeffShifting1 = PSystem->dri*PSystem->partDist/PSystem->pndSmallZero;				///< Coefficient used to adjust velocity type 1
	PSystem->coeffShifting2 = PSystem->coefA*PSystem->partDist*PSystem->partDist*PSystem->cflNumber*PSystem->machNumber;	///< Coefficient used to adjust velocity type 2
	PSystem->coeffPPE = 2.0*PSystem->dim/(PSystem->pndLargeZero*PSystem->lambdaZero);			///< Coefficient used to PPE
	PSystem->coeffPPESource = PSystem->relaxPND/(PSystem->timeStep*PSystem->timeStep*PSystem->pndSmallZero);	///< Coefficient used to PPE source term
	Dns[partType::FLUID]=PSystem->densityFluid;			Dns[partType::WALL]=PSystem->densityWall;
	invDns[partType::FLUID]=1.0/PSystem->densityFluid;	invDns[partType::WALL]=1.0/PSystem->densityWall;
	PSystem->invPartDist = 1.0/PSystem->partDist;
	PSystem->distCollisionLimit = PSystem->partDist*PSystem->distLimitRatio;	///< A distance that does not allow further access between particles
	PSystem->distCollisionLimit2 = PSystem->distCollisionLimit*PSystem->distCollisionLimit;
	PSystem->restitutionCollision = 1.0 + PSystem->collisionRatio;
	PSystem->numOfIterations = 0;										///< Number of iterations
	PSystem->fileNumber = 0;											///< File number
	PSystem->timeCurrent = 0.0;											///< Simulation time
	PSystem->velMax = 0.0;												///< Maximum flow velocity
	PSystem->CFLcurrent = PSystem->cflNumber;							///< Current Courant number
	PSystem->betaPnd = PSystem->pndThreshold*PSystem->pndSmallZero;		///< Surface cte PND
	PSystem->betaNeigh = PSystem->neighThreshold*PSystem->numNeighZero;	///< Surface cte Neighbors
	PSystem->delta2 = PSystem->npcdThreshold*PSystem->npcdThreshold*PSystem->partDist*PSystem->partDist;	///< Surface cte NPCD 
	PSystem->thetaArc = PSystem->thetaThreshold/180.0*3.14159265;				///< Surface cte theta ARC
	PSystem->hThreshold2 = 1.33*1.33*PSystem->partDist*PSystem->partDist;		///< Surface cte radius ARC
	PSystem->dstThreshold2 = 2.0*PSystem->hThreshold2;							///< Surface cte radius ARC
	PSystem->normThreshold2 = PSystem->normThreshold*PSystem->normThreshold;	///< Surface cte Normal
	
	//cout << "lo: " << partDist << " m, dt: " << timeStep << " s, PND0Small: " << pndSmallZero << " PND0Large: ";
	//cout << pndLargeZero << " PND0Grad: " << pndGradientZero << " lambda: " << lambdaZero << std::endl;
	//cout << "bPnd: " << betaPnd << "betaNeigh: " << betaNeigh << endl;
}

////////////////////////////////////////////////////////////
// Functions called during the simulation (main loop)
////////////////////////////////////////////////////////////

// Verify if particle is out of domain
void MpsParticle::checkParticleOutDomain(MpsParticleSystem *PSystem) {

#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		for(int b=0; b<PSystem->numBC; b++) {
			if(bucketTypeBC[b] == domainBC::PERIODIC) {
				// Periodic in X
				if(periodicDirection[b*3]) {
					if(pos[i*3  ]>PSystem->physDomMaxX) {
						pos[i*3  ] -= PSystem->periodicLength[0];
					}
					if(pos[i*3  ]<PSystem->physDomMinX) {
						pos[i*3  ] += PSystem->periodicLength[0];
					}
				}
				// Periodic in Y
				if(periodicDirection[b*3+1]) {
					if(pos[i*3+1]>PSystem->physDomMaxY) {
						pos[i*3+1] -= PSystem->periodicLength[1];
					}
					if(pos[i*3+1]<PSystem->physDomMinY) {
						pos[i*3+1] += PSystem->periodicLength[1];
					}
				}
				// Periodic in Z
				if(periodicDirection[b*3+2]) {
					if(pos[i*3+2]>PSystem->physDomMaxZ) {
						pos[i*3+2] -= PSystem->periodicLength[2];
					}
					if(pos[i*3+2]<PSystem->physDomMinZ) {
						pos[i*3+2] += PSystem->periodicLength[2];
					}
				}
			}
		}
	}

	// Auxliar variables
	int newNumRealParticles = numParticles;
	int newNumRealIOParticles = numRealAndIOParticles;

	double limMinX = PSystem->domainMinX;
	double limMaxX = PSystem->domainMaxX;
	double limMinY = PSystem->domainMinY;
	double limMaxY = PSystem->domainMaxY;
	double limMinZ = PSystem->domainMinZ;
	double limMaxZ = PSystem->domainMaxZ;
	// limitTypeBC = 0: Border particle positions
	// limitTypeBC = 1: Domain limits min and max
	////if(limitTypeBC == 0)
	////{
		limMinX += PSystem->bucketSide; limMaxX -= PSystem->bucketSide;
		limMinY += PSystem->bucketSide; limMaxY -= PSystem->bucketSide;
		limMinZ += PSystem->bucketSide; limMaxZ -= PSystem->bucketSide;
	////}
	for(int i=0; i<numParticles; i++) {
		if(	particleType[i] == PSystem->fluid && i < newNumRealParticles &&
			(pos[i*3  ]>limMaxX || pos[i*3  ]<limMinX ||
			pos[i*3+1]>limMaxY || pos[i*3+1]<limMinY ||
			(PSystem->dim == 3 && (pos[i*3+2]>limMaxZ || pos[i*3+2]<limMinZ)))) {
			
			// ID of last Real particle
			int iLastRealParticle = newNumRealParticles - 1;
			
			// Move the data from "Last Real Particle" to i-th ghost particle and 
			// update some data of "Last Real Particle" as ghost
			moveDataLastRealPartToGhostPart(i, iLastRealParticle, PSystem);

			// Decrease number of Real particles
			newNumRealParticles--;

			// Verify if the simulation has InOutflow bondary conditions
			if(PSystem->inOutflowOn == true && PSystem->numInOutflowPlane > 0) {
				
				// ID of last Real+IO particle
				int iLastRealIOParticle = newNumRealIOParticles - 1;

				// Move the data from "Last Real+IO Particle" to "Last Real Particle" and 
				// update some data of "Last Real+IO Particle" as ghost
				moveDataLastRealIOPartToLastRealPart(iLastRealParticle, iLastRealIOParticle, PSystem);
				
				// Decrease number of Real+IO particles
				newNumRealIOParticles--;
			}
			
		}
		
		// Update number of Real and Real+IO particles
		numParticles = newNumRealParticles;
		numRealAndIOParticles = newNumRealIOParticles;
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif
}

// Moves the data from "Last Real Particle" to i-th ghost particle and updates some data of "Last Particle" as ghost
void MpsParticle::moveDataLastRealPartToGhostPart(const int i, const int iLastReal, MpsParticleSystem *PSystem) {
	// Scalars
	particleType[i]=particleType[iLastReal];
	particleBC[i]=particleBC[iLastReal];
	numNeigh[i]=numNeigh[iLastReal];
	press[i]=press[iLastReal];
	pressAverage[i]=pressAverage[iLastReal];
	pndi[i]=pndi[iLastReal];
	pndki[i]=pndki[iLastReal];
	pndski[i]=pndski[iLastReal];
	pndSmall[i]=pndSmall[iLastReal];
	npcdDeviation2[i]=npcdDeviation2[iLastReal];
	concentration[i]=concentration[iLastReal];
	velDivergence[i]=velDivergence[iLastReal];
	diffusiveTerm[i]=diffusiveTerm[iLastReal];
	nearMeshType[i]=nearMeshType[iLastReal];
	particleNearWall[i]=particleNearWall[iLastReal];
	numNeighWallContribution[i]=numNeighWallContribution[iLastReal];
	pndWallContribution[i]=pndWallContribution[iLastReal];
	deviationDotPolygonNormal[i]=deviationDotPolygonNormal[iLastReal];
	numNeighborsSurfaceParticles[i]=numNeighborsSurfaceParticles[iLastReal];
	distParticleWall2[i]=distParticleWall2[iLastReal];
	PTYPE[i]=PTYPE[iLastReal];
	RHO[i]=RHO[iLastReal];
	MEU[i]=MEU[iLastReal];
	elementID[i]=elementID[iLastReal];

	if(PSystem->fluidType == viscType::NON_NEWTONIAN) {
		Cv[i]=Cv[iLastReal];
		II[i]=II[iLastReal];
		MEU_Y[i]=MEU_Y[iLastReal];
		Inertia[i]=Inertia[iLastReal];
		pnew[i]=pnew[iLastReal];
		p_rheo_new[i]=p_rheo_new[iLastReal];
		p_smooth[i]=p_smooth[iLastReal];
		VF[i]=VF[iLastReal];
		S12[i]=S12[iLastReal];
		S13[i]=S13[iLastReal];
		S23[i]=S23[iLastReal];
		S11[i]=S11[iLastReal];
		S22[i]=S22[iLastReal];
		S33[i]=S33[iLastReal];
	}

	if(PSystem->mpsType == calcPressType::IMPLICIT_PND || 
		PSystem->mpsType == calcPressType::IMPLICIT_PND_DIVU) {
		pressurePPE(i)=pressurePPE(iLastReal);
		sourceTerm(i)=sourceTerm(iLastReal);
	}
	
	if(PSystem->inOutflowOn == true && PSystem->numInOutflowPlane > 0) {
		signDist[i]=signDist[iLastReal];
		isInIORegion[i]=isInIORegion[iLastReal];
	}

	// Vectors
	for (int j = 0; j < 3; j++)
	{
		acc[i*3+j]=acc[iLastReal*3+j];
		accStar[i*3+j]=accStar[iLastReal*3+j];
		pos[i*3+j]=pos[iLastReal*3+j];
		vel[i*3+j]=vel[iLastReal*3+j];
		npcdDeviation[i*3+j]=npcdDeviation[iLastReal*3+j];
		gradConcentration[i*3+j]=gradConcentration[iLastReal*3+j];
		correcMatrixRow1[i*3+j]=correcMatrixRow1[iLastReal*3+j];
		correcMatrixRow2[i*3+j]=correcMatrixRow2[iLastReal*3+j];
		correcMatrixRow3[i*3+j]=correcMatrixRow3[iLastReal*3+j];
		normal[i*3+j]=normal[iLastReal*3+j];
		dvelCollision[i*3+j]=dvelCollision[iLastReal*3+j];
		particleAtWallPos[i*3+j]=particleAtWallPos[iLastReal*3+j];
		mirrorParticlePos[i*3+j]=mirrorParticlePos[iLastReal*3+j];
		wallParticleForce1[i*3+j]=wallParticleForce1[iLastReal*3+j];
		wallParticleForce2[i*3+j]=wallParticleForce2[iLastReal*3+j];
		polygonNormal[i*3+j]=polygonNormal[iLastReal*3+j];
		forceWall[i*3+j]=forceWall[iLastReal*3+j];
		//Posk[i*3+j]=Posk[iLastReal*3+j];
		//Velk[i*3+j]=Velk[iLastReal*3+j];
		//Acv[i*3+j]=Acv[iLastReal*3+j];
	}

	// Update some data of "Last Particle" that now represents the ghost
	particleType[iLastReal]=PSystem->ghost;
	particleBC[iLastReal]=PSystem->other;
	particleNearWall[iLastReal]=false;
	nearMeshType[iLastReal]=meshType::FIXED;
	distParticleWall2[iLastReal]=10e8*PSystem->partDist;
	// Set zero to velocity and press of lastParticle
	for (int j = 0; j < 3; j++) {
		//pos[iLastReal*3+j]=0.0;
		//Posk[iLastReal*3+j]=0.0;
		vel[iLastReal*3+j] = 0.0;
	}
	press[iLastReal] = 0.0;
	// Set maximum position to lastParticle
	pos[iLastReal*3  ] = PSystem->domainMaxX - PSystem->partDist;
	pos[iLastReal*3+1] = PSystem->domainMaxY - PSystem->partDist;
	if (PSystem->dim == 2) {
		pos[iLastReal*3+2] = 0.0;
	}
	else {
		pos[iLastReal*3+2] = PSystem->domainMaxZ - PSystem->partDist;
	}

	if(PSystem->inOutflowOn == true && PSystem->numInOutflowPlane > 0) {
		signDist[iLastReal]=10e8*PSystem->partDist;
		isInIORegion[iLastReal]=false;
	}
}

// Move the data from "Last Real+IO Particle" to "Last Real Particle" and update some data of "Last Real+IO Particle" as ghost
void MpsParticle::moveDataLastRealIOPartToLastRealPart(const int iLastReal, const int iLastRealIO, MpsParticleSystem *PSystem) {
	// Scalars
	particleType[iLastReal]=particleType[iLastRealIO];
	particleBC[iLastReal]=particleBC[iLastRealIO];
	numNeigh[iLastReal]=numNeigh[iLastRealIO];
	press[iLastReal]=press[iLastRealIO];
	pressAverage[iLastReal]=pressAverage[iLastRealIO];
	pndi[iLastReal]=pndi[iLastRealIO];
	pndki[iLastReal]=pndki[iLastRealIO];
	pndski[iLastReal]=pndski[iLastRealIO];
	pndSmall[iLastReal]=pndSmall[iLastRealIO];
	npcdDeviation2[iLastReal]=npcdDeviation2[iLastRealIO];
	concentration[iLastReal]=concentration[iLastRealIO];
	velDivergence[iLastReal]=velDivergence[iLastRealIO];
	diffusiveTerm[iLastReal]=diffusiveTerm[iLastRealIO];
	nearMeshType[iLastReal]=nearMeshType[iLastRealIO];
	particleNearWall[iLastReal]=particleNearWall[iLastRealIO];
	numNeighWallContribution[iLastReal]=numNeighWallContribution[iLastRealIO];
	pndWallContribution[iLastReal]=pndWallContribution[iLastRealIO];
	deviationDotPolygonNormal[iLastReal]=deviationDotPolygonNormal[iLastRealIO];
	numNeighborsSurfaceParticles[iLastReal]=numNeighborsSurfaceParticles[iLastRealIO];
	distParticleWall2[iLastReal]=distParticleWall2[iLastRealIO];
	PTYPE[iLastReal]=PTYPE[iLastRealIO];
	RHO[iLastReal]=RHO[iLastRealIO];
	MEU[iLastReal]=MEU[iLastRealIO];
	elementID[iLastReal]=elementID[iLastRealIO];

	if(PSystem->fluidType == viscType::NON_NEWTONIAN) {
		Cv[iLastReal]=Cv[iLastRealIO];
		II[iLastReal]=II[iLastRealIO];
		MEU_Y[iLastReal]=MEU_Y[iLastRealIO];
		Inertia[iLastReal]=Inertia[iLastRealIO];
		pnew[iLastReal]=pnew[iLastRealIO];
		p_rheo_new[iLastReal]=p_rheo_new[iLastRealIO];
		p_smooth[iLastReal]=p_smooth[iLastRealIO];
		VF[iLastReal]=VF[iLastRealIO];
		S12[iLastReal]=S12[iLastRealIO];
		S13[iLastReal]=S13[iLastRealIO];
		S23[iLastReal]=S23[iLastRealIO];
		S11[iLastReal]=S11[iLastRealIO];
		S22[iLastReal]=S22[iLastRealIO];
		S33[iLastReal]=S33[iLastRealIO];
	}

	// if(PSystem->mpsType == calcPressType::IMPLICIT_PND || 
	// 	PSystem->mpsType == calcPressType::IMPLICIT_PND_DIVU) {
	// 	pressurePPE(iLastReal)=pressurePPE(iLastRealIO);
	// 	sourceTerm(iLastReal)=sourceTerm(iLastRealIO);
	// }
	
	if(PSystem->inOutflowOn == true && PSystem->numInOutflowPlane > 0) {
		signDist[iLastReal]=signDist[iLastRealIO];
		isInIORegion[iLastReal]=isInIORegion[iLastRealIO];
	}

	// Vectors
	for (int j = 0; j < 3; j++)
	{
		acc[iLastReal*3+j]=acc[iLastRealIO*3+j];
		accStar[iLastReal*3+j]=accStar[iLastRealIO*3+j];
		pos[iLastReal*3+j]=pos[iLastRealIO*3+j];
		vel[iLastReal*3+j]=vel[iLastRealIO*3+j];
		npcdDeviation[iLastReal*3+j]=npcdDeviation[iLastRealIO*3+j];
		gradConcentration[iLastReal*3+j]=gradConcentration[iLastRealIO*3+j];
		correcMatrixRow1[iLastReal*3+j]=correcMatrixRow1[iLastRealIO*3+j];
		correcMatrixRow2[iLastReal*3+j]=correcMatrixRow2[iLastRealIO*3+j];
		correcMatrixRow3[iLastReal*3+j]=correcMatrixRow3[iLastRealIO*3+j];
		normal[iLastReal*3+j]=normal[iLastRealIO*3+j];
		dvelCollision[iLastReal*3+j]=dvelCollision[iLastRealIO*3+j];
		particleAtWallPos[iLastReal*3+j]=particleAtWallPos[iLastRealIO*3+j];
		mirrorParticlePos[iLastReal*3+j]=mirrorParticlePos[iLastRealIO*3+j];
		wallParticleForce1[iLastReal*3+j]=wallParticleForce1[iLastRealIO*3+j];
		wallParticleForce2[iLastReal*3+j]=wallParticleForce2[iLastRealIO*3+j];
		polygonNormal[iLastReal*3+j]=polygonNormal[iLastRealIO*3+j];
		forceWall[iLastReal*3+j]=forceWall[iLastRealIO*3+j];
		//Posk[iLastReal*3+j]=Posk[iLastRealIO*3+j];
		//Velk[iLastReal*3+j]=Velk[iLastRealIO*3+j];
		//Acv[iLastReal*3+j]=Acv[iLastRealIO*3+j];
	}

	// Update some data of "Last Particle" that now represents the ghost
	particleType[iLastRealIO]=PSystem->ghost;
	particleBC[iLastRealIO]=PSystem->other;
	particleNearWall[iLastRealIO]=false;
	nearMeshType[iLastRealIO]=meshType::FIXED;
	distParticleWall2[iLastRealIO]=10e8*PSystem->partDist;
	// Set zero to velocity and press of lastParticle
	for (int j = 0; j < 3; j++) {
		//pos[iLastRealIO*3+j]=0.0;
		//Posk[iLastRealIO*3+j]=0.0;
		vel[iLastRealIO*3+j] = 0.0;
	}
	press[iLastRealIO] = 0.0;
	// Set maximum position to lastParticle
	pos[iLastRealIO*3  ] = PSystem->domainMaxX - PSystem->partDist;
	pos[iLastRealIO*3+1] = PSystem->domainMaxY - PSystem->partDist;
	if (PSystem->dim == 2) {
		pos[iLastRealIO*3+2] = 0.0;
	}
	else {
		pos[iLastRealIO*3+2] = PSystem->domainMaxZ - PSystem->partDist;
	}

	if(PSystem->inOutflowOn == true && PSystem->numInOutflowPlane > 0) {
		signDist[iLastRealIO]=10e8*PSystem->partDist;
		isInIORegion[iLastRealIO]=false;
	}
}

// Set force on wall to zero
//void MpsParticle::WallZeroForce_omp(int nNodes, int nSolids, solid_fem * &solid) {
void MpsParticle::setWallForceZero(const int nNodes, double *nodeforceX, double *nodeforceY, double *nodeforceZ, MpsParticleSystem *PSystem) {
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == PSystem->fluid) {
			forceWall[i*3]=forceWall[i*3+1]=forceWall[i*3+2]=0.0;
		}
	}
#pragma omp parallel for
	for(int nn=0;nn<nNodes;nn++)
		nodeforceX[nn]=nodeforceY[nn]=nodeforceZ[nn]=0.0;
/*	for(int ss=0; ss<nSolids; ss++) {
//#pragma omp parallel for
		for(int ns=0; ns<solid[ss].nNodes; ns++) {
			solid[ss].node[ns].forceX = 0.0;
			solid[ss].node[ns].forceY = 0.0;
			solid[ss].node[ns].forceZ = 0.0;
		}}*/

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif
}

// // Force on wall due fluid particles - FSI
// void MpsParticle::forceParticlesToWall(mesh mesh, solid_fem * &solid) {
// 	//int nPartNearMesh = partNearMesh.size();

// 	double resForce_x, resForce_y, resForce_z, pForce_x, pForce_y, pForce_z;
// 	double AreaForce = pow(partDist,dim-1.0);
// 	resForce_x=0.0;	resForce_y=0.0;	resForce_z=0.0;
// 	pForce_x=0.0;	pForce_y=0.0;	pForce_z=0.0;

// 	//printf(" Mesh %d \n", nPartNearMesh);
// 	// Loop only for particles near mesh
// #pragma omp parallel for reduction(+: resForce_x, resForce_y, resForce_z, pForce_x, pForce_y, pForce_z)
// 	//for(int im=0;im<nPartNearMesh;im++) {
// 	//int i = partNearMesh[im];
// 	for(int i=0; i<numParticles; i++) {

// 	//if(particleType[i] == fluid && nearMeshType[i] == meshType::DEFORMABLE) {
// 	if(particleType[i] == fluid && particleNearWall[i] == true && nearMeshType[i] == meshType::DEFORMABLE) {
// 		//printf("\n%5d th timeCurrent: %lf / ENTROU !!! i: %d", numOfIterations, timeCurrent, i);
		
// 		resForce_x += forceWall[i*3]; resForce_y += forceWall[i*3+1]; resForce_z += forceWall[i*3+2];

// 		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
// 		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];

// 		double  normaliw[3], normalMod2;
// 		// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
// 		normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
// 		normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

// 		//if(i==16173)
// 		//   	printf("\ni:%5d th timeCurrent: %lf / ri: %lf %lf %lf / rim: %lf %lf %lf / N: %lf", i, timeCurrent, posXi, posYi, posZi, posMirrorXi, posMirrorYi, posMirrorZi, normalMod2);

// 		if(normalMod2 <= 0.26*partDist*partDist) {
// 			if(normalMod2 > epsilonZero) {
// 				double normalMod = sqrt(normalMod2);
// 				normaliw[0] = normaliw[0]/normalMod;
// 				normaliw[1] = normaliw[1]/normalMod;
// 				normaliw[2] = normaliw[2]/normalMod;
// 		}
// 		else {
// 			normaliw[0] = 0;
// 			normaliw[1] = 0;
// 			normaliw[2] = 0;
// 		}

// 			pForce_x += press[i]*normaliw[0]*AreaForce; pForce_y += press[i]*normaliw[1]*AreaForce; pForce_z += press[i]*normaliw[2]*AreaForce;
// 		}



// 		double FWG[3], FWL[3];	// Global and local force
// 		double V1[3],V2[3],V3[3],PW[3],RM[9];

// 		// Mesh element ID
// 		int elemID = elementID[i];

// 		// Triangle vertices
// 		int node1 = elemNode1id[elemID];
// 		int node2 = elemNode2id[elemID];
// 		int node3 = elemNode3id[elemID];

// 		int ss = 1;
// //		int node1 = solid[ss].element[elemID].node1ID;
// //		int node2 = solid[ss].element[elemID].node2ID;
// //		int node3 = solid[ss].element[elemID].node3ID;
// /*
// 		if(i==65)
// 		{
// 			std::cout << "NN: " << node1 << " X: " << nodeX[node1] << " Y: " << nodeY[node1] << " Z: " << nodeZ[node1] << std::endl;
// 			std::cout << "NN: " << node2 << " X: " << nodeX[node2] << " Y: " << nodeY[node2] << " Z: " << nodeZ[node2] << std::endl;
// 			std::cout << "NN: " << node3 << " X: " << nodeX[node3] << " Y: " << nodeY[node3] << " Z: " << nodeZ[node3] << std::endl;
// 		}
// 	*/
// 		V1[0] = nodeX[node1]; V1[1] = nodeY[node1]; V1[2] = nodeZ[node1];
// 		V2[0] = nodeX[node2]; V2[1] = nodeY[node2]; V2[2] = nodeZ[node2];
// 		V3[0] = nodeX[node3]; V3[1] = nodeY[node3]; V3[2] = nodeZ[node3];

// //		V1[0] = solid[ss].node[node1].x; V1[1] = solid[ss].node[node1].y; V1[2] = solid[ss].node[node1].z;
// //		V2[0] = solid[ss].node[node2].x; V2[1] = solid[ss].node[node2].y; V2[2] = solid[ss].node[node2].z;
// //		V3[0] = solid[ss].node[node3].x; V3[1] = solid[ss].node[node3].y; V3[2] = solid[ss].node[node3].z;
// /*
// 		if(i==65)
// 		{
// 			std::cout << "Elem: " << elemID << " N1: " << node1 << " N2: " << node2 << " N3: " << node3 << std::endl;
// 			std::cout << "N1 X: " << V1[0] << " Y: " << V1[1] << " Z: " << V1[2] << std::endl;
// 			std::cout << "N2 X: " << V2[0] << " Y: " << V2[1] << " Z: " << V2[2] << std::endl;
// 			std::cout << "N3 X: " << V3[0] << " Y: " << V3[1] << " Z: " << V3[2] << std::endl;
			
// 		}
// */
// 		// Transformation matrix
//		transformMatrix(V1, V2, V3, RM, PSystem->epsilonZero);

// 		if(i==65)
// 		{
// 			std::cout << "RM: " << RM[0] << " Y: " << RM[1] << " Z: " << RM[2] << std::endl;
// 			std::cout << "RM: " << RM[3] << " Y: " << RM[4] << " Z: " << RM[5] << std::endl;
// 			std::cout << "RM: " << RM[6] << " Y: " << RM[7] << " Z: " << RM[8] << std::endl;
			
// 		}

// 		// Wall point
// 		PW[0] = particleAtWallPos[i*3  ];	PW[1] = particleAtWallPos[i*3+1];	PW[2] = particleAtWallPos[i*3+2];
// /*
// 		if(i==65)
// 		{
// 			std::cout << "PW: " << PW[0] << " Y: " << PW[1] << " Z: " << PW[2] << std::endl;
// 		}
// */
// 		FWL[0] = forceWall[i*3]; FWL[1] = forceWall[i*3+1]; FWL[2] = forceWall[i*3+2];

// 		// Arbitrary plane to XY plane
// 		transform3Dto2D(V1, RM);
// 		transform3Dto2D(V2, RM);
// 		transform3Dto2D(V3, RM);
// 		transform3Dto2D(PW, RM);
// 		transform3Dto2D(FWL, RM);
// /*
// 		if(i==65)
// 		{
// 			std::cout << "PW: " << PW[0] << " Y: " << PW[1] << " Z: " << PW[2] << std::endl;
// 			std::cout << "FWL: " << FWL[0] << " Y: " << FWL[1] << " Z: " << FWL[2] << std::endl;
// 		}
// */
// 		// Element area 0.5*[(x2*y3 - x3*y2) + (y2 - y3)*x1 + (x3 - x2)*y1]
// 		double Ae = fabs(0.5*((V2[0]*V3[1] - V3[0]*V2[1]) + (V2[1] - V3[1])*V1[0] + (V3[0] - V2[0])*V1[1]));

// 		// Shape functions
// 		// N1 = Ar1/Ae, N2 = Ar2/Ae, N3 = Ar3/Ae
// 		// Ar1 = 0.5*[(x2*y3 - x3*y2) + (y2 - y3)*x + (x3 - x2)*y]
// 		// Ar2 = 0.5*[(x3*y1 - x1*y3) + (y3 - y1)*x + (x1 - x3)*y]
// 		// Ar3 = 0.5*[(x1*y2 - x2*y1) + (y1 - y2)*x + (x2 - x1)*y] = 1 - Ar2 - Ar1
// 		double Ar1 = fabs(0.5*((V2[0]*V3[1] - V3[0]*V2[1]) + (V2[1] - V3[1])*PW[0] + (V3[0] - V2[0])*PW[1]));
// 		double Ar2 = fabs(0.5*((V3[0]*V1[1] - V1[0]*V3[1]) + (V3[1] - V1[1])*PW[0] + (V1[0] - V3[0])*PW[1]));
// 		double Ar3 = Ae - Ar1 - Ar2;

// /*
// 		if(i==65)
// 		{
// 			std::cout << "Elem: " << elemID << " N1: " << node1 << " N2: " << node2 << " N3: " << node3 << std::endl;
// 			std::cout << "N1 X: " << V1[0] << " Y: " << V1[1] << " Z: " << V1[2] << std::endl;
// 			std::cout << "N2 X: " << V2[0] << " Y: " << V2[1] << " Z: " << V2[2] << std::endl;
// 			std::cout << "N3 X: " << V3[0] << " Y: " << V3[1] << " Z: " << V3[2] << std::endl;
// 			std::cout << "Ae: " << Ae << " Ar1: " << Ar1 << " Ar2: " << Ar2 << " Ar3: " << Ar3 << std::endl;
// 			std::cout << "t: " << timeCurrent << " fiX: " << forceWall[i*3] << " fiY: " << forceWall[i*3+1] << " fiZ: " << forceWall[i*3+2] << std::endl;
// 		}
// */
		
// 		// Nodal forces
// 		// Node 1
// 		FWG[0] = (Ar1/Ae)*FWL[0]; FWG[1] = (Ar1/Ae)*FWL[1]; FWG[2] = (Ar1/Ae)*FWL[2];
// 		// XY plane to Arbitrary plane
// 		transform2Dto3D(FWG, RM);
// 		nodeFx[node1] += FWG[0]; nodeFy[node1] += FWG[1]; nodeFz[node1] += FWG[2];
// //		solid[ss].node[node1].forceX += FWG[0]; solid[ss].node[node1].forceY += FWG[1]; solid[ss].node[node1].forceZ += FWG[2];
// /*		
// 		if(i==65)
// 		{
// 			std::cout << "t1: " << timeCurrent << " fWX: " << FWG[0] << " fWY: " << FWG[1] << " fWZ: " << FWG[2] << std::endl;
// 		}
// */
// 		// Node 2
// 		FWG[0] = (Ar2/Ae)*FWL[0]; FWG[1] = (Ar2/Ae)*FWL[1]; FWG[2] = (Ar2/Ae)*FWL[2];
// 		// XY plane to Arbitrary plane
// 		transform2Dto3D(FWG, RM);
// 		nodeFx[node2] += FWG[0]; nodeFy[node2] += FWG[1]; nodeFz[node2] += FWG[2];
// //		solid[ss].node[node2].forceX += FWG[0]; solid[ss].node[node2].forceY += FWG[1]; solid[ss].node[node2].forceZ += FWG[2];
// /*
// 		if(i==65)
// 		{
// 			std::cout << "t2: " << timeCurrent << " fWX: " << FWG[0] << " fWY: " << FWG[1] << " fWZ: " << FWG[2] << std::endl;
// 		}
// */
// 		// Node 3
// 		FWG[0] = (Ar3/Ae)*FWL[0]; FWG[1] = (Ar3/Ae)*FWL[1]; FWG[2] = (Ar3/Ae)*FWL[2];
// 		// XY plane to Arbitrary plane
// 		transform2Dto3D(FWG, RM);
// 		nodeFx[node3] += FWG[0]; nodeFy[node3] += FWG[1]; nodeFz[node3] += FWG[2];
// //		solid[ss].node[node3].forceX += FWG[0]; solid[ss].node[node3].forceY += FWG[1]; solid[ss].node[node3].forceZ += FWG[2];
// /*
// 		if(i==65)
// 		{
// 			std::cout << "t3: " << timeCurrent << " fWX: " << FWG[0] << " fWY: " << FWG[1] << " fWZ: " << FWG[2] << std::endl;
// 		}
// 		*/
// 	}}


// 	if(numOfIterations%1 == 0) {
// 		printf("\n%5d th timeCurrent: %lf / Fx: %lf / Fy: %lf / Fz: %lf", numOfIterations, timeCurrent, resForce_x, resForce_y, resForce_z);
// 		printf("\n%5d th timeCurrent: %lf / pFx: %lf / pFy: %lf / pFz: %lf", numOfIterations, timeCurrent, pForce_x, pForce_y, pForce_z);
// 	}
	
// 	// Open the File to write
// 	forceTxtFile = fopen(OUT_FORCE, "a");
// 	if(forceTxtFile == NULL) perror ("Error opening force txt file");
// 	fprintf(forceTxtFile,"\n%lf\t%lf\t%lf\t%lf",timeCurrent, resForce_x,resForce_y,resForce_z);
// 	// Close force file
// 	fclose(forceTxtFile);
// }