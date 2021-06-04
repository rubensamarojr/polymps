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

// Global variables
// Vector with ID of particles near to mesh
std::vector<int> partNearMesh;
// Forced rigid wall
double *nodeFRWX, *nodeFRWY, *nodeFRWZ; // Nodal positions

// Main Loop
// Class MpsParticle: Particles (mps)
// Class PolygonMesh: Polygons from the point of view of particles (mps)
void mainLoopOfSimulation(MpsParticle &part, PolygonMesh* &mesh) {

	// string -> char
	char *output_folder_char = new char[part.vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, part.vtuOutputFoldername.c_str());

	// Break if simulation reaches the final time
	while(true) {

		// Display simulation informations at each 10 iterations
		part.displayInfo(10);
		// Write output particle files
		part.writeOutputFiles();
		// Write output mesh files
		if(part.numOfIterations%part.iterOutput == 0) {
			if(part.femOn == true) {
				// Write deformable mesh (STL files)
				mesh[meshType::DEFORMABLE].writePolygonMeshFile(meshType::DEFORMABLE, output_folder_char, part.fileNumber);
			}
			if(part.forcedOn == true) {
				// Write forced rigid mesh (STL files)
				mesh[meshType::FORCED].writePolygonMeshFile(meshType::FORCED, output_folder_char, part.fileNumber);
			}

			part.fileNumber++; // Integer number
			
			// Break if simulation reachs the final time
			if(part.timeCurrent >= part.timeSimulation ) {
				delete[] output_folder_char;
				output_folder_char = NULL;
				break;
			}
		}
		///////////////////////////
		/// Numerical simulation //
		///////////////////////////
		// Update particle ID's in buckets
		part.updateBuckets();

		// if(part.fluidType == viscType::NON_NEWTONIAN) 
		// {
		// 	part.calcVolumeFraction(); // Volume of fraction if phase II in the mixture
		// 	part.calcViscosityInteractionVal();// Viscosity interaction values for "real" fluid particles
		// }

		// Calculation of acceleration due laplacian of viscosity and gravity
		part.calcViscosityGravity();
		// Add acceleration due pressure gradient (Prediction)
		if(part.relaxPress < 1.0) {
			part.predictionPressGradient();
			if(part.wallType == boundaryWallType::POLYGON) {
				part.predictionWallPressGradient(); // Pressure gradient on polygon wall
			}
		}
		// Update velocity and positions. Set some variables to zero or inf
		part.updateVelocityPosition1st();
		// Check collision between particles
		part.checkCollisions();

		// Contributions due polygon wall
		if(part.wallType == boundaryWallType::POLYGON) {
			// Fluid particles: Calculation of PND due wall and number of neighboors
			// Positions of wall and mirror particles
			for(int me = 0; me < part.numOfMeshs; me++) {
				mesh[me].closestPointPNDBoundaryAABB(part.reS2, part.reL2, part.numParticles, part.weightType, part.particleType, 
					part.fluid, part.ghost, me, meshType::FIXED, meshType::DEFORMABLE, meshType::FORCED, part.pos, part.particleAtWallPos, part.mirrorParticlePos, 
					part.distParticleWall2, part.elementID, part.nearMeshType, part.polygonNormal);
				// mesh[me].closestPointPNDBoundaryAABB(part.reS2, part.reL2, part.numParticles, part.weightType, part.particleType, 
				// 		part.fluid, part.ghost, part.pos, part.particleAtWallPos, part.mirrorParticlePos, part.distParticleWall2, 
				// 		part.pndWallContribution, part.numNeighWallContribution, elementID, partNearMesh);
			}
			partNearMesh.clear();
			//for(int me = 0; me < part.numOfMeshs; me++) {
			//	mesh[me].updateParticlesNearPolygonMesh(part.reS2, part.reL2, part.numParticles, part.weightType, part.particleType, 
			//	part.fluid, part.ghost, part.distParticleWall2, part.pndWallContribution, part.numNeighWallContribution, partNearMesh);
			//}
			// Only call once since the distance riw2 have the values from all the meshes
			mesh[0].updateParticlesNearPolygonMesh(part.reS2, part.reL2, part.numParticles, part.weightType, part.particleType, part.fluid, 
				part.ghost, part.distParticleWall2, part.pndWallContribution, part.numNeighWallContribution, partNearMesh, part.particleNearWall);

			// NPCD PND due polygon wall
			part.calcWallNPCD();
		}

		// PND calculation
		part.calcPnd();
		// Diffusion term
		if(part.pndType == calcPNDType::DIFFUSIVE) {
			part.calcPndDiffusiveTerm();
			if(part.wallType == boundaryWallType::POLYGON) {
				// Add contribution of PND due polygon wall
				if(part.slipCondition == slipBC::FREE_SLIP) {
					part.calcWallSlipPndDiffusiveTerm();
				}
				else if (part.slipCondition == slipBC::NO_SLIP) {
					part.calcWallNoSlipPndDiffusiveTerm();
				}
			}
			// Update PND
			part.updatePnd();
			// Mean PND at wall and dummy
			part.meanPndParticlesWallDummySurface();
		}
		// Mean PND
		else if(part.pndType == calcPNDType::MEAN_SUM_WIJ) {
			if(part.wallType == boundaryWallType::POLYGON) {
				part.meanWallPnd(); // Contribution due polygon wall
			}
			part.meanPnd();
		}

		// Mean fluid neighbor PND
		//part.meanNeighFluidPnd();

		// Update type of particle
		part.updateParticleBC();
		// Pressure calculation
		if(part.mpsType == calcPressType::EXPLICIT) {
			part.calcPressEMPS();
		}
		else if(part.mpsType == calcPressType::WEAKLY) {
			part.calcPressWCMPS();
		}
		else if(part.mpsType == calcPressType::IMPLICIT_PND)
		{
			part.solvePressurePoissonPnd();
		}
		else if(part.mpsType == calcPressType::IMPLICIT_PND_DIVU)
		{
			part.calcVelDivergence();
			if(part.wallType == boundaryWallType::POLYGON) {
				if(part.slipCondition == slipBC::FREE_SLIP) {
					part.calcWallSlipVelDivergence(); // Free-Slip condition
				}
				else if(part.slipCondition == slipBC::NO_SLIP) {
					part.calcWallNoSlipVelDivergence(); // No-Slip condition
				}
			}
			part.solvePressurePoissonPndDivU();
		}
		
		if(part.wallType == boundaryWallType::PARTICLE) {
			part.extrapolatePressParticlesWallDummy(); // Extrapolate pressure to wall and dummy particles
		}
		if(part.wallType == boundaryWallType::POLYGON) {
			// Pressure at inner particles near corners
			// Extrapolate pressure to inner particles near polygon walls (Not working !!)
			////part.extrapolatePressParticlesNearPolygonWall();
		}
		// Compute correction matrix
		if(part.gradientCorrection == true) {
			part.correctionMatrix();
		}
		// Calculation of acceleration due pressure gradient
		part.calcPressGradient();
		if(part.wallType == boundaryWallType::POLYGON) {
			// Add acceleration due pressure gradient on polygon wall
			part.calcWallPressGradient();
		}

		// Add acceleration due laplacian of viscosity on wall
		///////////////// Non-Newtonian flow /////////////////
		if(part.fluidType == viscType::NON_NEWTONIAN) {

			part.calcVolumeFraction(); // Calculation of the volume of fraction if phase II in the mixture
			part.calcViscosityInteractionVal(); // Viscosity interaction values for "real" fluid particles
			
			if(part.wallType == boundaryWallType::POLYGON) {
				if(part.slipCondition == slipBC::FREE_SLIP) {
					part.calcWallSlipViscosityInteractionVal(); // Free-slip condition. Viscosity interaction values
				}
				else if(part.slipCondition == slipBC::NO_SLIP) {
					part.calcWallNoSlipViscosityInteractionVal(); // No-slip condition. Viscosity interaction values
				}
			}
		}
		///////////////// Newtonian flow /////////////////////
		if(part.wallType == boundaryWallType::POLYGON) {
			if(part.slipCondition == slipBC::FREE_SLIP) {
				part.calcWallSlipViscosity(); // Free-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
			}
			else if(part.slipCondition == slipBC::NO_SLIP) {
				part.calcWallNoSlipViscosity(); // No-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
			}
		}
		if(part.forcedOn == true) {
			// Update forced mesh
			mesh[meshType::FORCED].updateForcedPolygonMesh(nodeFRWX, nodeFRWY, nodeFRWZ, part.velVWall, part.timeStep, part.timeCurrent);
		}
		
		// Update velocity and positions
		part.updateVelocityPosition2nd();
		// Shifting techniques
		// Adjust velocity
		if(part.shiftingType == 1) {
			//part.correctionMatrix();
			//part.calcNormalParticles();
			part.calcShifting();
			//part.calcWallShifting();
		}
		// Adjust position
		else if(part.shiftingType == 2) {
			//part.calcNormalConcentration();
			part.calcConcAndConcGradient();
			if(part.wallType == boundaryWallType::POLYGON) {
				part.calcWallConcAndConcGradient();
			}
		}
		// Wall and dummy velocity
		// if(part.wallType == boundaryWallType::PARTICLE) {
		// 	part.updateVelocityParticlesWallDummy();
		// }
		// Pressure calculation
		// if(part.mpsType == calcPressType::EXPLICIT) {
		// 	part.calcPressEMPSandParticleBC();
		// }
		// else if(part.mpsType == calcPressType::WEAKLY) {
		// 	part.calcPressWCMPSandParticleBC();
		// }

		// for(int i=0; i<part.numParticles; i++) {
		// 	part.pressAverage[i] += part.press[i];
		// }

		// Update iteration and time
		part.numOfIterations++;
		part.timeCurrent += part.timeStep;

	}
}

// Main
int main( int argc, char** argv) {

	Eigen::initParallel();

	printf("Start MPS.\n");
	//////////////////////////////
	////// Initializations////////
	//////////////////////////////

	// Create MpsParticle class
	//particle* particles = NULL;
	//particles = new particle[1];
	MpsParticle particles;

	// Read and allocate memory for data
	particles.readInputFile();
	// Write header of output txt files (force and pressure)
	// particles.writeHeaderTxtFiles(); (NOT WORKING !!!)

	// Allocation of buckets
	particles.allocateBuckets();
	// Setting parameters
	particles.setParameters();
	// Update particle ID's in buckets
	particles.updateBuckets();

	// Create an array of PolygonMesh class
	PolygonMesh* solidMesh = NULL;
	solidMesh = new PolygonMesh[particles.numOfMeshs];

	// It is necessary to read the msesh in the following order (maximum of 3 meshs)
	// 1st: Rigid mesh file
	// 2nd: Deformable mesh file
	// 3rd: Forced mesh file
	for(int me = 0; me < particles.numOfMeshs; me++) {
		solidMesh[me].initPolygonMesh(particles.numParticles);
		if(me == 0) {
			solidMesh[me].readPolygonMeshFile(particles.meshRigidFilename);
		}
		if(particles.femOn == true) {
			if(me == 1) {
				solidMesh[me].readPolygonMeshFile(particles.meshDeformableFilename);
			}
		}
		if(particles.forcedOn == true) {
			if(me == 2) {
				solidMesh[me].readPolygonMeshFile(particles.meshForcedFilename);
			}
		}
	}

	// Libigl remove duplicates
	std::cout << " Number of meshes: " << particles.numOfMeshs << std::endl;
	for(int me = 0; me < particles.numOfMeshs; me++) {
		std::cout << " after removing duplicates, mesh containts " << solidMesh[me].NV.rows() << " vertices and " << solidMesh[me].NF.rows() << " faces" << std::endl;
	}

	// Creation of the wall weight (Zij) and number of neighboors functions (numNeighWall)
	// std::cout << " reS: " << particles.reS << " reL: " << particles.reL << std::endl;
	for(int me = 0; me < particles.numOfMeshs; me++) {
		solidMesh[me].initWijnNeigh(particles.dim, particles.weightType, particles.partDist, particles.reL, particles.reS);
	}

	//////////////////////////////
	//////// Forced mesh /////////
	//////////////////////////////
	if(particles.forcedOn == true) {
	
		nodeFRWX = (double*)malloc(sizeof(double)*solidMesh[meshType::FORCED].NV.rows());	// Node position
		nodeFRWY = (double*)malloc(sizeof(double)*solidMesh[meshType::FORCED].NV.rows());	// Node position
		nodeFRWZ = (double*)malloc(sizeof(double)*solidMesh[meshType::FORCED].NV.rows());	// Node position

		for(int nn=0;nn<solidMesh[meshType::FORCED].NV.rows();nn++) {
			nodeFRWX[nn] = solidMesh[meshType::FORCED].NV(nn,0);
			nodeFRWY[nn] = solidMesh[meshType::FORCED].NV(nn,1);
			nodeFRWZ[nn] = solidMesh[meshType::FORCED].NV(nn,2);
			// std::cout << "N: " << nn << " X: " << nodeX[nn] << " Y: " << nodeY[nn] << " Z: " << nodeZ[nn] << std::endl;
		}

		particles.velVWall[0]=particles.velVWall[1]=particles.velVWall[2]=0.0;

	}

	//////////////////////////////
	////// Polygon wall mesh /////
	//////////////////////////////
	if(particles.wallType == boundaryWallType::POLYGON) {
		// Positions of wall and mirror particles
		for(int me = 0; me < particles.numOfMeshs; me++) {
			solidMesh[me].closestPointPNDBoundaryAABB(particles.reS2, particles.reL2, particles.numParticles, particles.weightType, 
				particles.particleType, particles.fluid, particles.ghost, me, meshType::FIXED, meshType::DEFORMABLE, 
				meshType::FORCED, particles.pos, particles.particleAtWallPos, particles.mirrorParticlePos, 
				particles.distParticleWall2, particles.elementID, particles.nearMeshType, particles.polygonNormal);
			// solidMesh[me].closestPointPNDBoundaryAABB(reS2, reL2, nP, WGT_TYP, 
			// Typ, FLD, GST, Pos, wallPos, mirrorPos, riw2, niw, numNeighw, elementID, partNearMesh);
		}
		partNearMesh.clear();
		//for(int me = 0; me < particles.numOfMeshs; me++) {
			// solidMesh[me].updateParticlesNearPolygonMesh(particles.reS2, particles.reL2, particles.numParticles, particles.weightType, 
					// Typ, FLD, GST, riw2, niw, numNeighw, partNearMesh);
		// }
		// Only call once since the distance riw2 have the values from all the meshes (Rigid, Deformable and Forced)
		solidMesh[0].updateParticlesNearPolygonMesh(particles.reS2, particles.reL2, particles.numParticles, particles.weightType, 
			particles.particleType, particles.fluid, particles.ghost, particles.distParticleWall2, particles.pndWallContribution,
			particles.numNeighWallContribution, partNearMesh, particles.particleNearWall);
		// NPCD PND due polygon wall
		particles.calcWallNPCD();
	}

	// Initial PND
	particles.setInitialPndNumberOfNeigh();
	if(particles.wallType == boundaryWallType::POLYGON) {
		// Contribution to mean PND due polygon wall
		particles.meanWallPnd();
	}
	// Mean of PND
	particles.meanPnd();

	// Mean fluid neighbor PND
	particles.meanNeighFluidPnd();

	// Update type of particle
	particles.updateParticleBC();
	// Compute pressure
	if(particles.mpsType == calcPressType::EXPLICIT) {
		particles.calcPressEMPS();
	}
	else if(particles.mpsType == calcPressType::WEAKLY) {
		particles.calcPressWCMPS();
	}
	else if(particles.mpsType == calcPressType::IMPLICIT_PND)
	{
		particles.solvePressurePoissonPnd();
	}
	else if(particles.mpsType == calcPressType::IMPLICIT_PND_DIVU)
	{
		particles.calcVelDivergence();
		if(particles.wallType == boundaryWallType::POLYGON) {
			if(particles.slipCondition == slipBC::FREE_SLIP) {
				particles.calcWallSlipVelDivergence(); // Free-Slip condition
			}
			else if(particles.slipCondition == slipBC::NO_SLIP) {
				particles.calcWallNoSlipVelDivergence(); // No-Slip condition
			}
		}
		particles.solvePressurePoissonPndDivU();
	}
	// Write header for vtu files
	particles.writePvd();

	///////////////////////////
	//////// Main loop ////////
	///////////////////////////
	mainLoopOfSimulation(particles, solidMesh);
	
	// Deallocate memory block
	free(nodeFRWX); free(nodeFRWY); free(nodeFRWZ);

	delete[] solidMesh;

	printf("End MPS.\n");
	return 0;
}
