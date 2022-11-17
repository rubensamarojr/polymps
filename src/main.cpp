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
	MpsParticleCollision* partColl, MpsBoundaryCondition* boundCond, MpsVectorMatrix* vectMatr, MpsParticleVelPos* partVelPos,
	MpsInflowOutflow* inOutflow);


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

	MpsInflowOutflow *inflowOutflow = nullptr;			///< MpsInflowOutflow engine

	// Creates MpsParticleSystem class
	particleSystem = new MpsParticleSystem();
	// Creates MpsParticle class
	particles = new MpsParticle();
	// Creates MpsInputOutput class
	inputOutput = new MpsInputOutput();
	// Reads, allocates memory and sets initial values for particles data
	try {
		inputOutput->readInputFile(particleSystem, particles);
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what() << std::endl;
		return -1;
	}
	
	// Write header of output txt files (force and pressure)
	// try {
	// 	inputOutput->writeHeaderTxtFiles(particles); // NOT WORKING !!!
	// }
	// catch (std::invalid_argument& e) {
	// 	std::cerr << e.what() << std::endl;
	// 	return -1;
	// }
	
	// Setting parameters
	particles->setParameters(particleSystem);
	// Creates MpsBuckets class
	buckets = new MpsBucket();
	// Allocation of memory for buckets data
	buckets->allocateMemory(particleSystem, particles);
	// Creates MpsBoundaryCondition class
	boundaryConditions = new MpsBoundaryCondition();
	// Set Periodic Boundary Condition of the bucket
	boundaryConditions->setBucketBC(particleSystem, particles);
	// Update particle ID's in buckets
	buckets->updateParticlesID(particleSystem, particles);
	// Set number of ghost particles to zero (in the current step)
	particles->numGhostParticles = 0;
	// Verify if particle is out of the domain, assign as Ghost and updates the current number of ghost particles
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

	// Creates MpsInflowOutflow class
	if(particleSystem->inOutflowOn == true && particleSystem->numInOutflowPlane > 0) {
		inflowOutflow = new MpsInflowOutflow[particleSystem->numInOutflowPlane];
		for (int ii = 0; ii < particleSystem->numInOutflowPlane; ii++) {
			inflowOutflow[ii].setPlaneFromPointNormal(particleSystem, ii);
		}
	}
	

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
	particlePress->calcPress(particleSystem, particles, buckets, inflowOutflow);
	// Writes header for vtu files
	inputOutput->writePvd(particleSystem);
	// Deletes all files inside the simulation folder
	inputOutput->deleteDirectoryFiles();
	// Writes VTK file of buckets
	inputOutput->writeBuckets(particleSystem);

	// Writes VTK file of initial Inflow/Outflow plans
	if(particleSystem->inOutflowOn == true && particleSystem->numInOutflowPlane > 0) {
		try {
			inputOutput->writeInOutFlowPlan(particleSystem, inflowOutflow);
		}
		catch (std::invalid_argument& e) {
			std::cerr << e.what() << std::endl;
			return -1;
		}
	}

	printf("OK\n");

	printf("\nStart Main Loop of Simulation\n\n");
	///////////////////////////
	//////// Main loop ////////
	///////////////////////////
	mainLoopOfSimulation(particleSystem, particles, solidMesh, inputOutput, 
		buckets, particlePndNeigh, particlePress, particleShifting, particleViscosity, 
		particleCollision, boundaryConditions, vectorMatrix, particleVelPos, 
		inflowOutflow);
	
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

	delete[] inflowOutflow;

	printf("\nEnd PolyMPS.\n");
	return 0;
}

// Main Loop
// Class MpsParticle: Particles (mps)
// Class PolygonMesh: Polygons from the point of view of particles (mps)
void mainLoopOfSimulation(MpsParticleSystem* partSyst, MpsParticle* part, PolygonMesh* mesh, MpsInputOutput* io, 
	MpsBucket* buck, MpsPndNeigh* partPndNeig, MpsPressure* partPress, MpsShifting* partShift, MpsViscosity* partVisc,
	MpsParticleCollision* partColl, MpsBoundaryCondition* boundCond, MpsVectorMatrix* vectMatr, MpsParticleVelPos* partVelPos,
	MpsInflowOutflow* inOutflow) {

	// string -> char
	char *output_folder_char = nullptr;
	io->stringToChar(output_folder_char);

	// Break if simulation reaches the final time
	while(true) {


		// // Display simulation informations at each 100 iterations
		// io->displayInfo(partSyst, part, partPress, 100);

		// // Write output vtk and stl files
		// if(partSyst->numOfIterations%partSyst->iterOutput == 0) {
		// 	// Write output particle files
		// 	io->writeOutputFiles(partSyst, part); ///< Write output particle files
		// 	if(partSyst->wallType == boundaryWallType::POLYGON) {
		// 		// Write output mesh files
		// 		if(partSyst->femOn == true) {
		// 			// Write deformable mesh (STL files)
		// 			mesh[meshType::DEFORMABLE].writePolygonMeshFile(meshType::DEFORMABLE, output_folder_char, partSyst->fileNumber);
		// 		}
		// 		if(partSyst->forcedOn == true) {
		// 			// Write forced rigid mesh (STL files)
		// 			mesh[meshType::FORCED].writePolygonMeshFile(meshType::FORCED, output_folder_char, partSyst->fileNumber);
		// 		}
		// 	}
		// 	partSyst->fileNumber++; // Integer number
		// }
		
		
		///////////////////////////
		/// Numerical simulation //
		///////////////////////////
		// Update particle ID's in buckets
		buck->updateParticlesID(partSyst, part);



		// Set number of ghost particles to zero (in the current step)
		part->numGhostParticles = 0;
		
		// Verify overlaped IO-IO or IO-Real particles, assign overlaped IO particles as Ghost 
		// and updates the current number of ghost particles
		if (partSyst->inOutflowOn == true && partSyst->numInOutflowPlane > 0) {
			if(part->numIOParticles > 0) {
				// Here only call once, (InOutflow plan of id = 0)
				inOutflow[0].checkOverlapedIOParticles(partSyst, part, buck);
			}
		}
		

		// Display simulation informations at each 100 iterations
		io->displayInfo(partSyst, part, partPress, 100);

		// Write output vtk and stl files
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


		// Non newtonian calculation
		///////////////// Non-Newtonian flow /////////////////
		if(partSyst->fluidType == viscType::NON_NEWTONIAN) {
			partVisc->calcVolumeFraction(partSyst, part, buck); ///< Volume of fraction if phase II in the mixture
			partVisc->calcViscosityInteractionVal(partSyst, part, buck); ///< Viscosity interaction values for "real" fluid particles
		}
		
		// Calculation of acceleration due laplacian of velocity and gravity
		partVisc->calcViscosityGravity(partSyst, part, buck);

		// Add acceleration due pressure gradient (Prediction)
		if(partSyst->relaxPress < 1.0) {
			partPress->predictionPressGradient(partSyst, part, buck);
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				partPress->predictionWallPressGradient(partSyst, part, buck); ///< Pressure gradient due to the polygon wall
			}
		}

		// Update velocity and positions. Set some variables to zero or inf
		partVelPos->updateVelocityPosition1st(partSyst, part);

		// Verify if particle is out of the domain, assign as Ghost and updates the current number of ghost particles
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
				mesh[me].closestPointPNDBoundaryAABB(part->numParticlesZero, part->particleType, partSyst->fluid, me,
					meshType::FIXED, meshType::DEFORMABLE, meshType::FORCED, part->pos, part->particleAtWallPos, part->mirrorParticlePos, 
					part->distParticleWall2, part->particleID, part->elementID, part->nearMeshType, part->polygonNormal);
				// mesh[me].closestPointPNDBoundaryAABB(partSyst->reS2, partSyst->reL2, part->numParticlesZero, partSyst->weightType, part->particleType, 
				// 		partSyst->fluid, part->pos, part->particleAtWallPos, part->mirrorParticlePos, part->distParticleWall2, 
				// 		part->pndWallContribution, part->numNeighWallContribution, elementID, partNearMesh);
			}
			partNearMesh.clear();
			//for(int me = 0; me < partSyst->numOfMeshs; me++) {
			//	mesh[me].updateParticlesNearPolygonMesh(partSyst->reS2, partSyst->reL2, part->numParticlesZero, partSyst->weightType, part->particleType, 
			//	partSyst->fluid, part->distParticleWall2, part->pndWallContribution, part->numNeighWallContribution, partNearMesh);
			//}
			// Only call once since the distance riw2 have the values from all the meshes
			mesh[0].updateParticlesNearPolygonMesh(partSyst->reS2, partSyst->reL2, part->numParticlesZero, part->particleType, partSyst->fluid, 
				part->distParticleWall2, part->pndWallContribution, part->numNeighWallContribution, part->particleID, part->particleNearWall, partNearMesh);

			partPndNeig->calcWallNPCD(partSyst, part, buck);	///< NPCD PND due to the polygon wall
		}

		// Compute correction matrix
		if(partSyst->gradientCorrection == true || partSyst->divergenceCorrection == true) {
			vectMatr->correctionMatrix(partSyst, part, buck);
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				vectMatr->correctionMatrixWall(partSyst, part, buck);
			}
			vectMatr->updateCorrectionMatrix(partSyst, part);
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
		partPress->calcPress(partSyst, part, buck, inOutflow);
		
		// Compute correction matrix
		if(partSyst->gradientCorrection == true || partSyst->divergenceCorrection == true) {
			vectMatr->correctionMatrix(partSyst, part, buck);
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				vectMatr->correctionMatrixWall(partSyst, part, buck);
			}
			vectMatr->updateCorrectionMatrix(partSyst, part);
		}
		
		// Minimum and Maximum pressure calculation
		if(partSyst->gradientType == 0 || partSyst->gradientType == 2 || partSyst->gradientType == 4) {
			partPress->calcPressMinMax(partSyst, part, buck);
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
		
		if(partSyst->forcedOn == true) {
			// Update forced mesh
			mesh[meshType::FORCED].updateForcedPolygonMesh(nodeFRWX, nodeFRWY, nodeFRWZ, partSyst->uniformVelWall, partSyst->timeStep, partSyst->timeCurrent);
		}


		// // Here only for InOutflow plan of id = 0
		// if(partSyst->inOutflowOn == true && partSyst->numInOutflowPlane > 0) {
		// 	// Impose motion to particles
		// 	inOutflow[0].imposeMotionParticles(partSyst, part);
		// }


		// Update velocity and positions
		partVelPos->updateVelocityPosition2nd(partSyst, part);
		
		// Verify if particle is out of the domain, assign as Ghost and updates the current number of ghost particles
		part->checkParticleOutDomain(partSyst);
		
		// Shifting techniques
		// Adjust velocity
		if(partSyst->shiftingType == 1) {
			//partShift->correctionMatrix(part);
			//partShift->calcNormalParticles(partSyst, part, buck);
			partShift->calcShifting(partSyst, part, buck);
			partShift->calcWallShifting(partSyst, part, buck);
			partShift->updateVelocity(partSyst, part);

			// partShift->calcVelGradient(partSyst, part, buck);
			// if(partSyst->wallType == boundaryWallType::POLYGON) {
			// 	if(partSyst->slipCondition == slipBC::FREE_SLIP) {
			// 		partShift->calcWallSlipVelGradient(partSyst, part, buck); // Free-Slip condition
			// 	}
			// 	else if(partSyst->slipCondition == slipBC::NO_SLIP) {
			// 		partShift->calcWallNoSlipVelGradient(partSyst, part, buck); // No-Slip condition
			// 	}
			// }
			// partShift->interpolateVelocity(partSyst, part);
		}
		// Adjust position
		else if(partSyst->shiftingType == 2) {
			// partShift->calcNormalConcentration(part);
			partShift->calcConcentration(partSyst, part, buck);
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				partShift->calcWallConcentration(partSyst, part, buck);
			}
			partShift->calcConcentrationGradient(partSyst, part, buck);
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				partShift->calcWallConcentrationGradient(partSyst, part, buck);
			}
			partShift->updatePosition(partSyst, part);

			partShift->calcVelGradient(partSyst, part, buck);
			if(partSyst->wallType == boundaryWallType::POLYGON) {
				if(partSyst->slipCondition == slipBC::FREE_SLIP) {
					partShift->calcWallSlipVelGradient(partSyst, part, buck); // Free-Slip condition
				}
				else if(partSyst->slipCondition == slipBC::NO_SLIP) {
					partShift->calcWallNoSlipVelGradient(partSyst, part, buck); // No-Slip condition
				}
			}
			partShift->interpolateVelocity(partSyst, part);
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
		// for(int ip=0; ip<part->numParticlesZero; ip++) {
		//	int ii = Particles->particleID[ip];
		// 	part->pressAverage[ii] += part->press[ii];
		// }
		


		// InOutflow bondary conditions
		if(partSyst->inOutflowOn == true && partSyst->numInOutflowPlane > 0) {
			// Set initial values to inOutflow variables
			inOutflow[0].setInOutflowVariables(partSyst, part);

			// Check particles in the Inflow/Outflow region, create real+IO particles, or delete real particles
			for(int ioID = 0; ioID < partSyst->numInOutflowPlane; ioID++) {
				// Check particles in the Inflow/Outflow region, create real and IO particles, 
				// or delete real particles
				inOutflow[ioID].checkCreateDeleteParticlesInOutflow(partSyst, part);

				// inOutflow[0].checkIOParticlesInOutflow(partSyst, part);
				
				// Verify if real particles were created
				if(part->realParticleCreated == true) {
					// Swap the data between Real particles in the array part->numRealAndIOParticles 
					// and the array part->numParticles
					inOutflow[ioID].swapIdRealAndIOParticlesInOutflow(partSyst, part);
				}
			}

			// // Move ghost particles to the last positions of the array
			// if(part->numDeletedParticles > 0) {
			// 	part->moveGhostToLastPosArray(partSyst);
			// }

			// Updates the number of ghost particles in the current step
			part->numGhostParticles += part->numDeletedParticles;

			// for(int ioID = 0; ioID < partSyst->numInOutflowPlane; ioID++) {
			// 	// Check particles in the Inflow/Outflow region
			// 	// inOutflow[ioID].checkRealParticlesInOutflow(partSyst, part);
			// 	// inOutflow[0].checkIOParticlesInOutflow(partSyst, part);
			// 	inOutflow[ioID].swapIdRealAndIOParticlesInOutflow(partSyst, part);
			// }
		}


				
		// Verify if particle is out of the domain, assign as Ghost and updates the current number of ghost particles
		part->checkParticleOutDomain(partSyst);

		// If the number of ghost particles is positive (in the current step), 
		// then move (swaps) ghost particles to the last positions of the array
		if(part->numGhostParticles > 0) {
			part->moveGhostToLastPosArray(partSyst);
		}

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
		mesh[me].initPolygonMesh(part->numParticlesZero);
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
			mesh[me].closestPointPNDBoundaryAABB(part->numParticlesZero, part->particleType, partSyst->fluid, me, meshType::FIXED,
				meshType::DEFORMABLE, meshType::FORCED, part->pos, part->particleAtWallPos, part->mirrorParticlePos,
				part->distParticleWall2, part->particleID, part->elementID, part->nearMeshType, part->polygonNormal);
			// mesh[me].closestPointPNDBoundaryAABB(reS2, reL2, nP, WGT_TYP, 
			// Typ, FLD, GST, Pos, wallPos, mirrorPos, riw2, niw, numNeighw, elementID, partNearMesh);
		}
		partNearMesh.clear();
		//for(int me = 0; me < partSyst->numOfMeshs; me++) {
			// mesh[me].updateParticlesNearPolygonMesh(partSyst->reS2, partSyst->reL2, part->numParticlesZero, partSyst->weightType, 
					// Typ, FLD, GST, riw2, niw, numNeighw, partNearMesh);
		// }
		// Only call once since the distance riw2 have the values from all the meshes (Rigid, Deformable and Forced)
		mesh[0].updateParticlesNearPolygonMesh(partSyst->reS2, partSyst->reL2, part->numParticlesZero,
			part->particleType, partSyst->fluid, part->distParticleWall2, part->pndWallContribution,
			part->numNeighWallContribution, part->particleID, part->particleNearWall, partNearMesh);
		
		// NPCD PND due polygon wall
		pndNeigh->calcWallNPCD(partSyst, part, buck);
	}	
}