#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include "PolygonMesh.h"
#include "MpsParticle.h"

// Grid file
//#define IN_FILE "input/dambreak_fluid.prof"
//#define IN_FILE "input/brumadinho_fluid_lo10p00.prof"
//#define IN_FILE "input/brumadinho_fluid_lo05p00_desl.prof"
//#define IN_FILE "input/dam1610_3D_fluid_lo0p010_mps.grid"
//#define IN_FILE "input/hydro_3D_lo0p010_fluid_mps.grid"

// Static rigid wall
//#define IN_MESH_0 "input/dam1610.stl"
//#define IN_MESH "input/BRUMADINHO_space10_model_lucas_top_0p50.stl"
//#define IN_MESH "input/damErosion3D.stl"
//#define IN_MESH "input/hydro_3D_box.stl"

// Output folder
//#define OUT_FOLDER "BRUMADINHO_lo10p00_rho1500_v04p00e-01"
//#define OUT_FOLDER "BRUMADINHO_lo05p00_rho1500_v04p00e-01"
//#define OUT_FOLDER "dam1610_NONEW_02"
//#define OUT_FOLDER "dambreak01"
//#define OUT_FOLDER "hydro_3D_lo0p010_03"

void initMesh(MpsParticle* part, PolygonMesh* mesh);
void mainLoopOfSimulation(MpsParticle* part, PolygonMesh* mesh);

// Global variables
// Particles
MpsParticle *particles = nullptr;
// Polygons mesh
PolygonMesh *solidMesh = nullptr;
// Vector with ID of particles near to mesh
std::vector<int> partNearMesh;
// Forced rigid wall
double *nodeFRWX, *nodeFRWY, *nodeFRWZ; // Nodal positions

// Main
int main( int argc, char** argv) {

	printf("Start PolyMPS.\n");

	Eigen::initParallel();

	//////////////////////////////
	////// Initializations ///////
	//////////////////////////////
	// Create MpsParticle class
	particles = new MpsParticle();
	particles->init();

	if(particles->mpsType == calcPressType::IMPLICIT_PND || particles->mpsType == calcPressType::IMPLICIT_PND_DIVU) {
		if(particles->pndType == calcPNDType::DIFFUSIVE || particles->pndType == calcPNDType::MEAN_SUM_WIJ) {
			printf("\nError: Please, set 'pnd-type to 0' in the json file for incompressible MPS\n");
			throw std::exception();
		}
	}

	printf("Particle class initialized.\n");

	if(particles->wallType == boundaryWallType::POLYGON) {
		// Create an array of PolygonMesh class
		solidMesh = new PolygonMesh[particles->numOfMeshs];
		initMesh(particles, solidMesh);

		printf("Mesh class initialized.\n");
	}

	printf("Initial Step... ");
	// Update variables at 0th step
	particles->stepZero();
	printf("OK\n");

	printf("\nStart Main Loop of Simulation\n\n");
	///////////////////////////
	//////// Main loop ////////
	///////////////////////////
	mainLoopOfSimulation(particles, solidMesh);
	
	// Deallocate memory block
	free(nodeFRWX); free(nodeFRWY); free(nodeFRWZ);
	delete particles;
	delete[] solidMesh;

	printf("End PolyMPS.\n");
	return 0;
}

void initMesh(MpsParticle* part, PolygonMesh* mesh){
	
	printf("Reading STL File... \n");
	// It is necessary to read the msesh in the following order (maximum of 3 meshs)
	// 1st: Rigid mesh file
	// 2nd: Deformable mesh file
	// 3rd: Forced mesh file
	for(int me = 0; me < part->numOfMeshs; me++) {
		mesh[me].initPolygonMesh(part->numParticles);
		if(me == 0) {
			mesh[me].readPolygonMeshFile(part->meshRigidFilename);
		}
		if(part->femOn == true) {
			if(me == 1) {
				mesh[me].readPolygonMeshFile(part->meshDeformableFilename);
			}
		}
		if(part->forcedOn == true) {
			if(me == 2) {
				mesh[me].readPolygonMeshFile(part->meshForcedFilename);
			}
		}
	}

	// Libigl remove duplicates
	std::cout << " Number of meshes: " << part->numOfMeshs << std::endl;
	for(int me = 0; me < part->numOfMeshs; me++) {
		std::cout << " after removing duplicates, mesh containts " << mesh[me].NV.rows() << " vertices and " << mesh[me].NF.rows() << " faces" << std::endl;
	}
	printf("...OK\n");

	// Creation of the wall weight (Zij) and number of neighboors functions (numNeighWall)
	// std::cout << " reS: " << part->reS << " reL: " << part->reL << std::endl;
	for(int me = 0; me < part->numOfMeshs; me++) {
		mesh[me].initWijnNeigh((int)part->dim, part->weightType, part->partDist, part->reL, part->reS);
	}

	//////////////////////////////
	//////// Forced mesh /////////
	//////////////////////////////
	if(part->forcedOn == true) {
	
		nodeFRWX = (double*)malloc(sizeof(double)*mesh[meshType::FORCED].NV.rows());	// Node position
		nodeFRWY = (double*)malloc(sizeof(double)*mesh[meshType::FORCED].NV.rows());	// Node position
		nodeFRWZ = (double*)malloc(sizeof(double)*mesh[meshType::FORCED].NV.rows());	// Node position

		for(int nn=0;nn<mesh[meshType::FORCED].NV.rows();nn++) {
			nodeFRWX[nn] = mesh[meshType::FORCED].NV(nn,0);
			nodeFRWY[nn] = mesh[meshType::FORCED].NV(nn,1);
			nodeFRWZ[nn] = mesh[meshType::FORCED].NV(nn,2);
			// std::cout << "N: " << nn << " X: " << nodeX[nn] << " Y: " << nodeY[nn] << " Z: " << nodeZ[nn] << std::endl;
		}

		part->velVWall[0]=part->velVWall[1]=part->velVWall[2]=0.0;

	}

	//////////////////////////////
	////// Polygon wall mesh /////
	//////////////////////////////
	if(part->wallType == boundaryWallType::POLYGON) {
		// Positions of wall and mirror particles
		for(int me = 0; me < part->numOfMeshs; me++) {
			mesh[me].closestPointPNDBoundaryAABB(part->reS2, part->reL2, part->numParticles, part->weightType, 
				part->particleType, part->fluid, me, meshType::FIXED, meshType::DEFORMABLE, 
				meshType::FORCED, part->pos, part->particleAtWallPos, part->mirrorParticlePos, 
				part->distParticleWall2, part->elementID, part->nearMeshType, part->polygonNormal);
			// mesh[me].closestPointPNDBoundaryAABB(reS2, reL2, nP, WGT_TYP, 
			// Typ, FLD, GST, Pos, wallPos, mirrorPos, riw2, niw, numNeighw, elementID, partNearMesh);
		}
		partNearMesh.clear();
		//for(int me = 0; me < part->numOfMeshs; me++) {
			// mesh[me].updateParticlesNearPolygonMesh(part->reS2, part->reL2, part->numParticles, part->weightType, 
					// Typ, FLD, GST, riw2, niw, numNeighw, partNearMesh);
		// }
		// Only call once since the distance riw2 have the values from all the meshes (Rigid, Deformable and Forced)
		mesh[0].updateParticlesNearPolygonMesh(part->reS2, part->reL2, part->numParticles, part->weightType, 
			part->particleType, part->fluid, part->distParticleWall2, part->pndWallContribution,
			part->numNeighWallContribution, partNearMesh, part->particleNearWall);
		// NPCD PND due polygon wall
		part->calcWallNPCD();
	}	
}

// Main Loop
// Class MpsParticle: Particles (mps)
// Class PolygonMesh: Polygons from the point of view of particles (mps)
void mainLoopOfSimulation(MpsParticle* part, PolygonMesh* mesh) {

	// string -> char
	char *output_folder_char = new char[part->vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, part->vtuOutputFoldername.c_str());

	// Break if simulation reaches the final time
	while(true) {

		// Display simulation informations at each 100 iterations
		part->displayInfo(100);
		if(part->numOfIterations%part->iterOutput == 0) {
			// Write output particle files
			part->writeOutputFiles();
			if(part->wallType == boundaryWallType::POLYGON) {
				// Write output mesh files
				if(part->femOn == true) {
					// Write deformable mesh (STL files)
					mesh[meshType::DEFORMABLE].writePolygonMeshFile(meshType::DEFORMABLE, output_folder_char, part->fileNumber);
				}
				if(part->forcedOn == true) {
					// Write forced rigid mesh (STL files)
					mesh[meshType::FORCED].writePolygonMeshFile(meshType::FORCED, output_folder_char, part->fileNumber);
				}
			}
			part->fileNumber++; // Integer number
		}
		///////////////////////////
		/// Numerical simulation //
		///////////////////////////
		// Update particle ID's in buckets
		part->updateBuckets();
		// Non newtonian calculation
		if(part->fluidType == viscType::NON_NEWTONIAN) {
			part->calcVolumeFraction(); // Volume of fraction if phase II in the mixture
			part->calcViscosityInteractionVal();// Viscosity interaction values for "real" fluid particles
		}
		// Calculation of acceleration due laplacian of viscosity and gravity
		part->calcViscosityGravity();
		// Add acceleration due pressure gradient (Prediction)
		if(part->relaxPress < 1.0) {
			part->predictionPressGradient();
			if(part->wallType == boundaryWallType::POLYGON) {
				part->predictionWallPressGradient(); // Pressure gradient on polygon wall
			}
		}
		// Update velocity and positions. Set some variables to zero or inf
		part->updateVelocityPosition1st();
		// Verify if particle is out of the domain
		part->checkParticleOutDomain();
		// Check collision between particles
		if(part->collisionType == colType::PC) {
			part->checkParticleCollisions();
		}
		else {
			part->checkDynamicParticleCollisions();
		}
		// Contributions due polygon wall
		if(part->wallType == boundaryWallType::POLYGON) {
			// Fluid particles: Calculation of PND due wall and number of neighboors
			// Positions of wall and mirror particles
			for(int me = 0; me < part->numOfMeshs; me++) {
				mesh[me].closestPointPNDBoundaryAABB(part->reS2, part->reL2, part->numParticles, part->weightType, part->particleType, 
					part->fluid, me, meshType::FIXED, meshType::DEFORMABLE, meshType::FORCED, part->pos, part->particleAtWallPos, part->mirrorParticlePos, 
					part->distParticleWall2, part->elementID, part->nearMeshType, part->polygonNormal);
				// mesh[me].closestPointPNDBoundaryAABB(part->reS2, part->reL2, part->numParticles, part->weightType, part->particleType, 
				// 		part->fluid, part->pos, part->particleAtWallPos, part->mirrorParticlePos, part->distParticleWall2, 
				// 		part->pndWallContribution, part->numNeighWallContribution, elementID, partNearMesh);
			}
			partNearMesh.clear();
			//for(int me = 0; me < part->numOfMeshs; me++) {
			//	mesh[me].updateParticlesNearPolygonMesh(part->reS2, part->reL2, part->numParticles, part->weightType, part->particleType, 
			//	part->fluid, part->distParticleWall2, part->pndWallContribution, part->numNeighWallContribution, partNearMesh);
			//}
			// Only call once since the distance riw2 have the values from all the meshes
			mesh[0].updateParticlesNearPolygonMesh(part->reS2, part->reL2, part->numParticles, part->weightType, part->particleType, part->fluid, 
				part->distParticleWall2, part->pndWallContribution, part->numNeighWallContribution, partNearMesh, part->particleNearWall);

			// NPCD PND due polygon wall
			part->calcWallNPCD();
		}
		// Compute correction matrix
		if(part->gradientCorrection == true) {
			part->correctionMatrix();
		}
		// PND, number of neighbors and NPCD calculation
		part->calcPndnNeighNPCD();
		// Diffusion term
		if(part->pndType == calcPNDType::DIFFUSIVE) {
			part->calcPndDiffusiveTerm();
			if(part->wallType == boundaryWallType::POLYGON) {
				// Add contribution of PND due polygon wall
				if(part->slipCondition == slipBC::FREE_SLIP) {
					part->calcWallSlipPndDiffusiveTerm();
				}
				else if (part->slipCondition == slipBC::NO_SLIP) {
					part->calcWallNoSlipPndDiffusiveTerm();
				}
			}
			// Update PND
			part->updatePnd();
			// Mean PND at wall and dummy
			part->meanPndParticlesWallDummySurface();
		}
		// Mean PND
		else if(part->pndType == calcPNDType::MEAN_SUM_WIJ) {
			if(part->wallType == boundaryWallType::POLYGON) {
				part->meanWallPnd(); // Contribution due polygon wall
			}
			part->meanPnd();
		}
		// Mean fluid neighbor PND
		//part->meanNeighFluidPnd();
		// Update type of particle
		if(part->freeSurfType == calcBCType::PND_ARC) {
			// Compute fluid particles normal vector
			part->calcNormalParticles();
			if(part->wallType == boundaryWallType::POLYGON) {
				// Contribution to normal vector due polygon wall
				part->calcWallNormalParticles();
			}
		}
		part->updateParticleBC();
		// Pressure calculation
		if(part->mpsType == calcPressType::EXPLICIT) {
			part->calcPressEMPS();
		}
		else if(part->mpsType == calcPressType::WEAKLY) {
			part->calcPressWCMPS();
		}
		else if(part->mpsType == calcPressType::IMPLICIT_PND) {
			part->solvePressurePoissonPnd();
		}
		else if(part->mpsType == calcPressType::IMPLICIT_PND_DIVU) {
			part->calcVelDivergence();
			if(part->wallType == boundaryWallType::POLYGON) {
				if(part->slipCondition == slipBC::FREE_SLIP) {
					part->calcWallSlipVelDivergence(); // Free-Slip condition
				}
				else if(part->slipCondition == slipBC::NO_SLIP) {
					part->calcWallNoSlipVelDivergence(); // No-Slip condition
				}
			}
			part->solvePressurePoissonPndDivU();
		}
		// Extrapolate pressure to wall and dummy particles
		if(part->wallType == boundaryWallType::PARTICLE) {
			part->extrapolatePressParticlesWallDummy();
		}
		if(part->wallType == boundaryWallType::POLYGON) {
			// Pressure at inner particles near corners
			// Extrapolate pressure to inner particles near polygon walls (Not working !!)
			////part->extrapolatePressParticlesNearPolygonWall();
		}
		// Compute correction matrix
		if(part->gradientCorrection == true) {
			part->correctionMatrix();
		}
		// Calculation of acceleration due pressure gradient
		part->calcPressGradient();
		if(part->wallType == boundaryWallType::POLYGON) {
			// Add acceleration due pressure gradient on polygon wall
			part->calcWallPressGradient();
		}
		// Add acceleration due laplacian of viscosity on wall
		///////////////// Non-Newtonian flow /////////////////
		if(part->fluidType == viscType::NON_NEWTONIAN) {

			part->calcVolumeFraction(); // Calculation of the volume of fraction if phase II in the mixture
			part->calcViscosityInteractionVal(); // Viscosity interaction values for "real" fluid particles
			
			if(part->wallType == boundaryWallType::POLYGON) {
				if(part->slipCondition == slipBC::FREE_SLIP) {
					part->calcWallSlipViscosityInteractionVal(); // Free-slip condition. Viscosity interaction values
				}
				else if(part->slipCondition == slipBC::NO_SLIP) {
					part->calcWallNoSlipViscosityInteractionVal(); // No-slip condition. Viscosity interaction values
				}
			}
		}
		///////////////// Newtonian flow /////////////////////
		if(part->wallType == boundaryWallType::POLYGON) {
			if(part->slipCondition == slipBC::FREE_SLIP) {
				part->calcWallSlipViscosity(); // Free-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
			}
			else if(part->slipCondition == slipBC::NO_SLIP) {
				part->calcWallNoSlipViscosity(); // No-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
			}
		}
		if(part->forcedOn == true) {
			// Update forced mesh
			mesh[meshType::FORCED].updateForcedPolygonMesh(nodeFRWX, nodeFRWY, nodeFRWZ, part->velVWall, part->timeStep, part->timeCurrent);
		}
		// Update velocity and positions
		part->updateVelocityPosition2nd();
		// Verify if particle is out of the domain
		part->checkParticleOutDomain();
		// Shifting techniques
		// Adjust velocity
		if(part->shiftingType == 1) {
			//part->correctionMatrix();
			//part->calcNormalParticles();
			part->calcShifting();
			//part->calcWallShifting();
		}
		// Adjust position
		else if(part->shiftingType == 2) {
			//part->calcNormalConcentration();
			part->calcConcAndConcGradient();
			if(part->wallType == boundaryWallType::POLYGON) {
				part->calcWallConcAndConcGradient();
			}
		}
		// Wall and dummy velocity
		if(part->wallType == boundaryWallType::PARTICLE) {
			part->updateVelocityParticlesWallDummy();
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
		part->checkParticleOutDomain();

		// Update iteration and time
		part->numOfIterations++;
		part->timeCurrent += part->timeStep;

		// Break if simulation reachs the final time
		if(part->timeCurrent > part->timeSimulation ) {
			delete[] output_folder_char;
			output_folder_char = NULL;
			break;
		}
	}
}