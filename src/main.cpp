/**
 * @defgroup   MAIN main
 *
 * @brief      This file is part of PolyMPS, a open-soruce CFD code based on MPS
 * Copyright (C) 2022  Rubens Augusto Amaro Junior <rubens.amaro@usp.br>
 * 
 * @date       2022
 */

/** @file
 * @brief The PolyMPS starting point, also known as the host.
 */

/** @mainpage PolyMPS
 * PolyMPS is a C++ code for numerical modeling of free-surface flow using the mesh-less moving particle 
 * semi-implicit (MPS) method from two different ways, namely, the incompressible and weakly-compressible
 * approaches. Besides the conventional solid wall boundaries represented by discrete layers of wall and dummy
 * (ghost) particles, boundary walls can be modeled by polygonal mesh of triangles, therefore providing flexibility
 * to model the surface of three-dimensional complex-shaped bodies. The present code has a lightweight and
 * efficient CPU implementation using OpenMP for research and practical applications.
 *
 * <hr>
 *
 */

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <limits>

#include "PolygonMesh.h"
#include "MpsParticle.h"
#include "MpsParticleSystem.h"
#include "MpsViscosity.h"
#include "MpsInputOutput.h"
#include "MpsParticleCollision.h"
#include "MpsShifting.h"
#include "MpsPressure.h"
#include "MpsPndNeigh.h"
#include "MpsBucket.h"
#include "MpsBoundaryCondition.h"
#include "MpsVectorMatrix.h"
#include "MpsParticleVelPos.h"

void initMesh(MpsParticleSystem* partSyst, MpsParticle* part, PolygonMesh* mesh, 
	MpsPndNeigh* pndNeigh, MpsInputOutput* io, MpsBucket* buck);
void mainLoopOfSimulation(MpsParticleSystem* partSyst, MpsParticle* part, PolygonMesh* mesh, MpsInputOutput* io, 
	MpsBucket* buck, MpsPndNeigh* partPndNeig, MpsPressure* partPress, MpsShifting* partShift, MpsViscosity* partVisc,
	MpsParticleCollision* partColl, MpsBoundaryCondition* boundCond, MpsVectorMatrix* vectMatr, MpsParticleVelPos* partVelPos);


// Global variables
std::vector<int> partNearMesh; 						///< Vector with ID of particles near to mesh
// Forced rigid wall
double *nodeFRWX, *nodeFRWY, *nodeFRWZ; 			///< Nodal positions

/**
 * @brief      The source code begins execution
 *
 * @param[in]  argc  An integer that contains the count of arguments that follow in argv. 
 * The argc parameter is always greater than or equal to 1.
 * @param      argv  An array of null-terminated strings representing command-line arguments 
 * 					entered by the user of the program. By convention, argv[0] is the command with which the 
 *					program is invoked. argv[1] is the first command-line argument. The last argument from the 
 * 					command line is argv[argc - 1], and argv[argc] is always NULL.
 *
 * @return     return 0 in the main function means that the program executed successfully
 */
int main( int argc, char** argv) {

	printf("Start PolyMPS.\n");

	Eigen::initParallel();

	//////////////////////////////
	////// Initializations ///////
	//////////////////////////////
	MpsParticle *particles = nullptr;					///< MpsParticle engine
	PolygonMesh *solidMesh = nullptr;					///< PolygonMesh engine
	MpsParticleSystem *particleSystem = nullptr;		///< MpsParticleSystem engine
	MpsViscosity *particleViscosity = nullptr;			///< MpsViscosity engine
	MpsInputOutput *inputOutput = nullptr;				///< MpsInputOutput engine
	MpsParticleCollision *particleCollision = nullptr;	///< MpsParticleCollision engine
	MpsShifting *particleShifting = nullptr;			///< MpsShifting engine
	MpsPressure *particlePress = nullptr;				///< MpsPressure engine
	MpsPndNeigh *particlePndNeigh = nullptr;			///< MpsPndNeigh engine
	MpsBucket *buckets = nullptr;						///< MpsBucket engine
	MpsBoundaryCondition *boundaryConditions = nullptr;	///< MpsBoundaryCondition engine
	MpsVectorMatrix *vectorMatrix = nullptr;			///< MpsVectorMatrix engine
	MpsParticleVelPos *particleVelPos = nullptr;		///< MpsParticleVelPos engine

	// Creates MpsParticleSystem class
	particleSystem = new MpsParticleSystem();
	// Creates MpsParticle class
	particles = new MpsParticle();
	// Creates MpsInputOutput class
	inputOutput = new MpsInputOutput();
	// Reads and allocates memory for particles data
	inputOutput->readInputFile(particleSystem, particles);
	
	// Write header of output txt files (force and pressure)
	//inputOutput->writeHeaderTxtFiles(particles); (NOT WORKING !!!)
	
	// Setting parameters
	particles->setParameters(particleSystem);
	// Creates MpsBuckets class
	buckets = new MpsBucket();
	// Allocation of memory for buckets data
	buckets->allocateMemory(particleSystem, particles);
	// Creates MpsBoundaryCondition class
	boundaryConditions = new MpsBoundaryCondition();
	// Set Periodic Boundary Condition of the bucket
	boundaryConditions->setBucketBC(particleSystem, particles, buckets);
	// Update particle ID's in buckets
	buckets->updateParticlesID(particleSystem, particles);
	// Verify if particle is out of domain
	particles->checkParticleOutDomain(particleSystem);

	printf("Particle class initialized.\n");

	// Creates MpsPndNeigh class
	particlePndNeigh = new MpsPndNeigh();

	if(particleSystem->wallType == boundaryWallType::POLYGON) {
		// Creates an array of PolygonMesh class
		solidMesh = new PolygonMesh[particleSystem->numOfMeshs];
		initMesh(particleSystem, particles, solidMesh, particlePndNeigh, inputOutput, buckets);

		printf("Mesh class initialized.\n");
	}

	// Creates MpsViscosity class
	particleViscosity = new MpsViscosity();
	// Creates MpsParticleCollision class
	particleCollision = new MpsParticleCollision();
	// Creates MpsShifting class
	particleShifting = new MpsShifting();
	// Creates MpsPressure class
	particlePress = new MpsPressure();
	if(particleSystem->mpsType == calcPressType::IMPLICIT_PND || particleSystem->mpsType == calcPressType::IMPLICIT_PND_DIVU) {
		particlePress->allocateMemory(particles);
	}
	// Creates MpsVectorMatrix class
	vectorMatrix = new MpsVectorMatrix();
	// Creates MpsParticleVelPos class
	particleVelPos = new MpsParticleVelPos();

	printf("Initial Step... ");
	// Updates variables at 0th step
	
	// Initial PND
	particlePndNeigh->setInitialPndNumberOfNeighNPCD(particleSystem, particles, buckets);
	if(particleSystem->wallType == boundaryWallType::POLYGON) {
		// Contribution to mean PND due polygon wall
		particlePndNeigh->meanWallPnd(particleSystem, particles, buckets);
	}
	// Mean of PND
	particlePndNeigh->meanPnd(particleSystem, particles, buckets);
	// Mean fluid neighbor PND
	particlePndNeigh->meanNeighFluidPnd(particleSystem, particles, buckets);
	
	// Update type of particle
	if(particleSystem->freeSurfType == calcBCType::PND_ARC) {
		// Compute fluid particles normal vector
		particleShifting->calcNormalParticles(particleSystem, particles, buckets);
		if(particleSystem->wallType == boundaryWallType::POLYGON) {
			// Contribution to normal vector due polygon wall
			particleShifting->calcWallNormalParticles(particleSystem, particles, buckets);
		}
	}

	boundaryConditions->updateParticleBC(particleSystem, particles, buckets);
	// Computes pressure
	particlePress->calcPress(particleSystem, particles, buckets);
	// Writes header for vtu files
	inputOutput->writePvd(particleSystem, particles);
	// Deletes all files inside the simulation folder
	inputOutput->deleteDirectoryFiles();
	// Writes VTK file of buckets
	inputOutput->writeBuckets(particleSystem, particles);

	printf("OK\n");

	printf("\nStart Main Loop of Simulation\n\n");
	///////////////////////////
	//////// Main loop ////////
	///////////////////////////
	mainLoopOfSimulation(particleSystem, particles, solidMesh, inputOutput, 
		buckets, particlePndNeigh, particlePress, particleShifting, particleViscosity, 
		particleCollision, boundaryConditions, vectorMatrix, particleVelPos);
	
	// Deallocate memory block
	free(nodeFRWX); free(nodeFRWY); free(nodeFRWZ);
	delete particles;
	delete[] solidMesh;
	delete particleSystem;
	delete particleViscosity;
	delete inputOutput;
	delete particleCollision;
	delete particleShifting;
	delete particlePress;
	delete particlePndNeigh;
	delete buckets;
	delete boundaryConditions;
	delete vectorMatrix;
	delete particleVelPos;

	printf("End PolyMPS.\n");
	return 0;
}

// Main Loop
// Class MpsParticle: Particles (mps)
// Class PolygonMesh: Polygons from the point of view of particles (mps)
void mainLoopOfSimulation(MpsParticleSystem* partSyst, MpsParticle* part, PolygonMesh* mesh, MpsInputOutput* io, 
	MpsBucket* buck, MpsPndNeigh* partPndNeig, MpsPressure* partPress, MpsShifting* partShift, MpsViscosity* partVisc,
	MpsParticleCollision* partColl, MpsBoundaryCondition* boundCond, MpsVectorMatrix* vectMatr, MpsParticleVelPos* partVelPos) {

	// string -> char
	char *output_folder_char = nullptr;
	io->stringToChar(output_folder_char);

	// Break if simulation reaches the final time
	while(true) {

		// Display simulation informations at each 100 iterations
		io->displayInfo(partSyst, part, partPress, 100);
		if(partSyst->numOfIterations%partSyst->iterOutput == 0) {
			// Write output particle files
			io->writeOutputFiles(partSyst, part); ///< Write output particle files
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				// Write output mesh files
				if(partSyst->femOn == true) {
					// Write deformable mesh (STL files)
					mesh[meshType::DEFORMABLE].writePolygonMeshFile(meshType::DEFORMABLE, output_folder_char, partSyst->fileNumber);
				}
				if(partSyst->forcedOn == true) {
					// Write forced rigid mesh (STL files)
					mesh[meshType::FORCED].writePolygonMeshFile(meshType::FORCED, output_folder_char, partSyst->fileNumber);
				}
			}
			partSyst->fileNumber++; // Integer number
		}

		///////////////////////////
		/// Numerical simulation //
		///////////////////////////
		
		// Update particle ID's in buckets
		buck->updateParticlesID(partSyst, part);

		// Non newtonian calculation
		if(partSyst->fluidType == viscType::NON_NEWTONIAN) {
			partVisc->calcVolumeFraction(partSyst, part, buck); ///< Volume of fraction if phase II in the mixture
			partVisc->calcViscosityInteractionVal(partSyst, part, buck); ///< Viscosity interaction values for "real" fluid particles
		}

		// Calculation of acceleration due laplacian of velocity and gravity
		partVisc->calcViscosity(partSyst, part, buck);

		// Add gravity to acceleration
		part->addGravity(partSyst, part);

		// Add acceleration due pressure gradient (Prediction)
		if(partSyst->relaxPress < 1.0) {
			partPress->predictionPressGradient(partSyst, part, buck);
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				partPress->predictionWallPressGradient(partSyst, part, buck); ///< Pressure gradient due to the polygon wall
			}
		}
		
		// Prediction of velocities and positions. Set some variables to zero or inf
		partVelPos->predictionVelocityPosition(partSyst, part);
		
		// Verify if particle is out of the domain
		part->checkParticleOutDomain(partSyst);
		
		// Check collision between particles
		if(partSyst->collisionType == colType::PC) {
			partColl->checkParticleCollisions(partSyst, part, buck);
		}
		else {;
			partColl->checkDynamicParticleCollisions(partSyst, part, buck);
		}

		// Contributions due polygon wall
		if(partSyst->wallType == boundaryWallType::POLYGON) {
			// Fluid particles: Calculation of PND due wall and number of neighboors
			// Positions of wall and mirror particles
			for(int me = 0; me < partSyst->numOfMeshs; me++) {
				mesh[me].closestPointPNDBoundaryAABB(partSyst->reS2, partSyst->reL2, part->numParticles, partSyst->weightType, part->particleType, 
					partSyst->fluid, me, meshType::FIXED, meshType::DEFORMABLE, meshType::FORCED, part->pos, part->particleAtWallPos, part->mirrorParticlePos, 
					part->distParticleWall2, part->elementID, part->nearMeshType, part->polygonNormal);
				// mesh[me].closestPointPNDBoundaryAABB(partSyst->reS2, partSyst->reL2, part->numParticles, partSyst->weightType, part->particleType, 
				// 		partSyst->fluid, part->pos, part->particleAtWallPos, part->mirrorParticlePos, part->distParticleWall2, 
				// 		part->pndWallContribution, part->numNeighWallContribution, elementID, partNearMesh);
			}
			partNearMesh.clear();
			//for(int me = 0; me < partSyst->numOfMeshs; me++) {
			//	mesh[me].updateParticlesNearPolygonMesh(partSyst->reS2, partSyst->reL2, part->numParticles, partSyst->weightType, part->particleType, 
			//	partSyst->fluid, part->distParticleWall2, part->pndWallContribution, part->numNeighWallContribution, partNearMesh);
			//}
			// Only call once since the distance riw2 have the values from all the meshes
			mesh[0].updateParticlesNearPolygonMesh(partSyst->reS2, partSyst->reL2, part->numParticles, partSyst->weightType, part->particleType, partSyst->fluid, 
				part->distParticleWall2, part->pndWallContribution, part->numNeighWallContribution, partNearMesh, part->particleNearWall);

			partPndNeig->calcWallNPCD(partSyst, part, buck);	///< NPCD PND due to the polygon wall
		}

		// Compute correction matrix
		if(partSyst->gradientCorrection == true) {
			vectMatr->correctionMatrix(partSyst, part, buck);
		}

		// PND, number of neighbors and NPCD calculation
		partPndNeig->calcPndnNeighNPCD(partSyst, part, buck);
		// Diffusion term
		if(partSyst->pndType == calcPNDType::DIFFUSIVE) {
			partPndNeig->calcPndDiffusiveTerm(partSyst, part, buck);
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				// Add contribution of PND due polygon wall
				if(partSyst->slipCondition == slipBC::FREE_SLIP) {
					partPndNeig->calcWallSlipPndDiffusiveTerm(partSyst, part, buck);
				}
				else if (partSyst->slipCondition == slipBC::NO_SLIP) {
					partPndNeig->calcWallNoSlipPndDiffusiveTerm(partSyst, part, buck);
				}
			}
			// Update PND
			partPndNeig->updatePnd(part);
			// Mean PND at wall and dummy
			partPndNeig->meanPndParticlesWallDummySurface(partSyst, part, buck);
		}
		// Mean PND
		else if(partSyst->pndType == calcPNDType::MEAN_SUM_WIJ) {
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				partPndNeig->meanWallPnd(partSyst, part, buck);	///< Contribution due to the polygon wall
			}
			partPndNeig->meanPnd(partSyst, part, buck);
		}
		// Mean fluid neighbor PND
		///partPndNeig->meanNeighFluidPnd(part, part, buck);
		
		// Update type of particle
		if(partSyst->freeSurfType == calcBCType::PND_ARC) {
			// Compute fluid particles normal vector
			partShift->calcNormalParticles(partSyst, part, buck);
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				// Contribution to normal vector due polygon wall
				partShift->calcWallNormalParticles(partSyst, part, buck);
			}
		}
		boundCond->updateParticleBC(partSyst, part, buck);

		// Pressure calculation
		partPress->calcPress(partSyst, part, buck);

		// Extrapolate pressure to wall and dummy particles
		if(partSyst->wallType == boundaryWallType::PARTICLE) {
			partPress->extrapolatePressParticlesWallDummy(partSyst, part, buck);
		}
		if(partSyst->wallType == boundaryWallType::POLYGON) {
			// Pressure at inner particles near corners
			// Extrapolate pressure to inner particles near polygon walls (Not working !!)
			///partPress->extrapolatePressParticlesNearPolygonWall(partSyst, part, buck);
		}

		// Compute correction matrix
		if(partSyst->gradientCorrection == true) {
			vectMatr->correctionMatrix(partSyst, part, buck);
		}

		// Calculation of acceleration due pressure gradient
		partPress->calcPressGradient(partSyst, part, buck);	///< Add acceleration due pressure gradient 
		if(partSyst->wallType == boundaryWallType::POLYGON) {
			partPress->calcWallPressGradient(partSyst, part, buck);	///< Add acceleration from pressure gradient due to the polygon wall
		}

		// Add acceleration due laplacian of viscosity on wall
		///////////////// Non-Newtonian flow /////////////////
		if(partSyst->fluidType == viscType::NON_NEWTONIAN) {
			partVisc->calcVolumeFraction(partSyst, part, buck); ///< Calculation of the volume of fraction if phase II in the mixture
			partVisc->calcViscosityInteractionVal(partSyst, part, buck); ///< Viscosity interaction values for "real" fluid particles
			
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				if(partSyst->slipCondition == slipBC::FREE_SLIP) {
					partVisc->calcWallSlipViscosityInteractionVal(partSyst, part, buck); ///< Free-slip condition. Viscosity interaction values
				}
				else if(partSyst->slipCondition == slipBC::NO_SLIP) {
					partVisc->calcWallNoSlipViscosityInteractionVal(partSyst, part, buck); ///< No-slip condition. Viscosity interaction values
				}
			}
		}
		///////////////// Newtonian flow /////////////////////
		if(partSyst->wallType == boundaryWallType::POLYGON) {
			if(partSyst->slipCondition == slipBC::FREE_SLIP) {
				partVisc->calcWallSlipViscosity(partSyst, part, buck); ///< Free-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
			}
			else if(partSyst->slipCondition == slipBC::NO_SLIP) {
				partVisc->calcWallNoSlipViscosity(partSyst, part, buck); ///< No-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
			}
		}

		// Update forced solid wall modeled by polygons
		if(partSyst->forcedOn == true) {
			// Update forced mesh
			mesh[meshType::FORCED].updateForcedPolygonMesh(nodeFRWX, nodeFRWY, nodeFRWZ, partSyst->uniformVelWall, partSyst->timeStep, partSyst->timeCurrent);
		}
		
		// Correction of velocities and positions
		partVelPos->correctionVelocityPosition(partSyst, part);
		
		// Verify if particle is out of the domain
		part->checkParticleOutDomain(partSyst);
		
		// Shifting techniques
		// Adjust velocity
		if(partSyst->shiftingType == 1) {
			//partShift->correctionMatrix(part);
			//partShift->calcNormalParticles(partSyst, part, buck);
			partShift->calcShifting(partSyst, part, buck);
			//partShift->calcWallShifting(partSyst, part, buck);
		}
		// Adjust position
		else if(partSyst->shiftingType == 2) {
			//partShift->calcNormalConcentration(part);
			partShift->calcConcAndConcGradient(partSyst, part, buck);
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				partShift->calcWallConcAndConcGradient(partSyst, part, buck);
			}
		}

		// Wall and dummy velocity
		if(partSyst->wallType == boundaryWallType::PARTICLE) {
			partVelPos->updateVelocityParticlesWallDummy(partSyst, part, buck);
		}
		
		// Pressure calculation
		// if(part->mpsType == calcPressType::EXPLICIT) {
		// 	part->calcPressEMPSandParticleBC();
		// }
		// else if(part->mpsType == calcPressType::WEAKLY) {
		// 	part->calcPressWCMPSandParticleBC();
		// }
		// for(int i=0; i<part->numParticles; i++) {
		// 	part->pressAverage[i] += part->press[i];
		// }
		
		// Verify if particle is out of the domain
		part->checkParticleOutDomain(partSyst);

		// Update iteration and time
		partSyst->numOfIterations++;
		partSyst->timeCurrent += partSyst->timeStep;

		// Break if simulation reachs the final time
		if(partSyst->timeCurrent > partSyst->timeSimulation ) {
			delete[] output_folder_char;
			output_folder_char = NULL;
			break;
		}
	}
}

void initMesh(MpsParticleSystem* partSyst, MpsParticle* part, PolygonMesh* mesh,
	MpsPndNeigh* pndNeigh, MpsInputOutput* io, MpsBucket* buck){
	
	printf("Reading STL File... \n");
	// It is necessary to read the msesh in the following order (maximum of 3 meshs)
	// 1st: Rigid mesh file
	// 2nd: Deformable mesh file
	// 3rd: Forced mesh file
	for(int me = 0; me < partSyst->numOfMeshs; me++) {
		mesh[me].initPolygonMesh(part->numParticles);
		if(me == 0) {
			mesh[me].readPolygonMeshFile(io->meshRigidFilename);
		}
		if(partSyst->femOn == true) {
			if(me == 1) {
				mesh[me].readPolygonMeshFile(io->meshDeformableFilename);
			}
		}
		if(partSyst->forcedOn == true) {
			if(me == 2) {
				mesh[me].readPolygonMeshFile(io->meshForcedFilename);
			}
		}
	}

	// Libigl remove duplicates
	std::cout << " Number of meshes: " << partSyst->numOfMeshs << std::endl;
	for(int me = 0; me < partSyst->numOfMeshs; me++) {
		std::cout << " after removing duplicates, mesh containts " << mesh[me].NV.rows() << " vertices and " << mesh[me].NF.rows() << " faces" << std::endl;
	}
	printf("...OK\n");

	// Creation of the wall weight (Zij) and number of neighboors functions (numNeighWall)
	// std::cout << " reS: " << partSyst->reS << " reL: " << partSyst->reL << std::endl;
	for(int me = 0; me < partSyst->numOfMeshs; me++) {
		mesh[me].initWijnNeigh((int)partSyst->dim, partSyst->weightType, partSyst->partDist, partSyst->reL, partSyst->reS);
	}

	//////////////////////////////
	//////// Forced mesh /////////
	//////////////////////////////
	if(partSyst->forcedOn == true) {
	
		nodeFRWX = (double*)malloc(sizeof(double)*mesh[meshType::FORCED].NV.rows());	// Node position
		nodeFRWY = (double*)malloc(sizeof(double)*mesh[meshType::FORCED].NV.rows());	// Node position
		nodeFRWZ = (double*)malloc(sizeof(double)*mesh[meshType::FORCED].NV.rows());	// Node position

		for(int nn=0;nn<mesh[meshType::FORCED].NV.rows();nn++) {
			nodeFRWX[nn] = mesh[meshType::FORCED].NV(nn,0);
			nodeFRWY[nn] = mesh[meshType::FORCED].NV(nn,1);
			nodeFRWZ[nn] = mesh[meshType::FORCED].NV(nn,2);
			// std::cout << "N: " << nn << " X: " << nodeX[nn] << " Y: " << nodeY[nn] << " Z: " << nodeZ[nn] << std::endl;
		}

		// Hard-code. Need to be implemented using json input variable
		partSyst->uniformVelWall[0]=partSyst->uniformVelWall[1]=partSyst->uniformVelWall[2]=0.0;

	}

	//////////////////////////////
	////// Polygon wall mesh /////
	//////////////////////////////
	if(partSyst->wallType == boundaryWallType::POLYGON) {
		// Positions of wall and mirror particles
		for(int me = 0; me < partSyst->numOfMeshs; me++) {
			mesh[me].closestPointPNDBoundaryAABB(partSyst->reS2, partSyst->reL2, part->numParticles, partSyst->weightType, 
				part->particleType, partSyst->fluid, me, meshType::FIXED, meshType::DEFORMABLE, 
				meshType::FORCED, part->pos, part->particleAtWallPos, part->mirrorParticlePos, 
				part->distParticleWall2, part->elementID, part->nearMeshType, part->polygonNormal);
			// mesh[me].closestPointPNDBoundaryAABB(reS2, reL2, nP, WGT_TYP, 
			// Typ, FLD, GST, Pos, wallPos, mirrorPos, riw2, niw, numNeighw, elementID, partNearMesh);
		}
		partNearMesh.clear();
		//for(int me = 0; me < partSyst->numOfMeshs; me++) {
			// mesh[me].updateParticlesNearPolygonMesh(partSyst->reS2, partSyst->reL2, part->numParticles, partSyst->weightType, 
					// Typ, FLD, GST, riw2, niw, numNeighw, partNearMesh);
		// }
		// Only call once since the distance riw2 have the values from all the meshes (Rigid, Deformable and Forced)
		mesh[0].updateParticlesNearPolygonMesh(partSyst->reS2, partSyst->reL2, part->numParticles, partSyst->weightType, 
			part->particleType, partSyst->fluid, part->distParticleWall2, part->pndWallContribution,
			part->numNeighWallContribution, partNearMesh, part->particleNearWall);
		
		// NPCD PND due polygon wall
		pndNeigh->calcWallNPCD(partSyst, part, buck);
	}	
}