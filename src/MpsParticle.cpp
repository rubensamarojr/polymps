// Copyright (c) 2021 Rubens AMARO
// Distributed under the MIT License.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
// strings and c-strings
#include <iostream>
#include <cstring>
#include <string>
// input file
#include "json.hpp"
// mkdir for Linux
#include <sys/stat.h>
// mkdir for Windows
#if defined(_WIN32) || defined(WIN32) || defined(__MINGW32__) || defined(__BORLANDC__)
#include <direct.h>
#endif
#include <sys/time.h>
#include "MpsParticle.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/LU>

#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;

// Constructor declaration
MpsParticle::MpsParticle()
{
}
// Destructor declaration
MpsParticle::~MpsParticle()
{
}

// Return time
double MpsParticle::getTime() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 1.0e-6);
}

void MpsParticle::displayInfo(const int intervalIter) {
	if(numOfIterations%intervalIter == 0) {
		timerEnd = getTime();
		int seconds, hours, minutes;
		seconds = int(timerEnd - timerStart);
		//seconds = int(timer_end - timer_sta);
		minutes = seconds / 60;
		hours = minutes / 60;
		printf("Iteration: %5dth Time: %lfsec Num. Particles: %d Max Velocity: %lfm/s Courant: %lf", 
			numOfIterations, timeCurrent, numParticles, velMax, CFLcurrent);
		if(mpsType == calcPressType::IMPLICIT_PND || mpsType == calcPressType::IMPLICIT_PND_DIVU){
			printf(" Solver iterations: %3d Estimated error: %.2e", solverIter, solverError);
		}
		printf(" RunTime: %02dh%02dm%02dsec\n", int(hours), int(minutes%60), int(seconds%60));
	}
}

// Initialize elements of the class
void MpsParticle::init() {
	// Read and allocate memory for data
	readInputFile();
	// Write header of output txt files (force and pressure)
	// part->writeHeaderTxtFiles(); (NOT WORKING !!!)
	// Allocation of buckets
	allocateBuckets();
	// Setting parameters
	setParameters();
	// Update particle ID's in buckets
	if((int)dim == 2) {
		updateBuckets2D();
	}
	else {
		updateBuckets3D();
	}
}

// Update variables at 0th step
void MpsParticle::stepZero() {

	// Initial PND
	setInitialPndNumberOfNeigh();
	if(wallType == boundaryWallType::POLYGON) {
		// Contribution to mean PND due polygon wall
		meanWallPnd();
	}
	// Mean of PND
	meanPnd();
	// Mean fluid neighbor PND
	meanNeighFluidPnd();
	// Update type of particle
	if(freeSurfType == calcBCType::PND_ARC) {
		// Compute fluid particles normal vector
		calcNormalParticles();
		if(wallType == boundaryWallType::POLYGON) {
			// Contribution to normal vector due polygon wall
			calcWallNormalParticles();
		}
	}
	updateParticleBC();
	// Compute pressure
	if(mpsType == calcPressType::EXPLICIT) {
		calcPressEMPS();
	}
	else if(mpsType == calcPressType::WEAKLY) {
		calcPressWCMPS();
	}
	else if(mpsType == calcPressType::IMPLICIT_PND)
	{
		solvePressurePoissonPnd();
	}
	else if(mpsType == calcPressType::IMPLICIT_PND_DIVU)
	{
		calcVelDivergence();
		if(wallType == boundaryWallType::POLYGON) {
			if(slipCondition == slipBC::FREE_SLIP) {
				calcWallSlipVelDivergence(); // Free-Slip condition
			}
			else if(slipCondition == slipBC::NO_SLIP) {
				calcWallNoSlipVelDivergence(); // No-Slip condition
			}
		}
		solvePressurePoissonPndDivU();
	}
	// Write header for vtu files
	writePvd();
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

// Return the bucket coordinates for particle "i"
void MpsParticle::bucketCoordinates(int &bx, int &by, int &bz,
	const double rxi, const double ryi, const double rzi) {
	bx = (int)((rxi - domainMinX)*invBucketSide) + 1;
	by = (int)((ryi - domainMinY)*invBucketSide) + 1;
	bz = (int)((rzi - domainMinZ)*invBucketSide) + 1;
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
		default:
			return -re/(dst*dst);
	}
}

////////////////////////////////////////////////////////////
// Functions called only at the initial instant (t=0)
////////////////////////////////////////////////////////////

// Read input data from file .json to class MpsParticle
void MpsParticle::readInputFile() {

	// Runtime start
	timerStart = getTime();

	char json_folder[] = "input/";
	char json_file_char [1000];
	char json_path_char [1000];
	bool readOK = false;
	printf(" ____________________________________________________________ \n");
	printf("|                                                            |\n");
	printf("|                      E-MPS/WC-MPS/MPS                      |\n");
	printf("|         Explicit/Weakly-compressible/Semi-implicit         |\n");
	printf("|                Moving Particle Simulation                  |\n");
	printf("|                                                            |\n");
	printf("|               University of Sao Paulo - Brazil             |\n");
	printf("|                                                            |\n");
	printf("|  by Rubens Augusto Amaro Junior                            |\n");
	printf("|____________________________________________________________|\n\n");
	while (readOK == false)
	{
		printf("Enter the name of the MPS input file:\n");
		scanf("%s", json_file_char);
		printf("\n");

		// char *json_file_char = new char[json_file.length()+1];
		// strcpy (json_file_char, json_file.c_str());
		// json_file_char now contains a c-string copy of json_file
		strcat(json_file_char, ".json");
		snprintf(json_path_char, 1000, "%s%s", json_folder, json_file_char);

		//tries to read the input json file
		js = fopen(json_path_char, "r");
		if (js == NULL) {
			printf("Error reading the input file %s. Try again.\n", json_path_char);
		}
		else {
			readOK = true;
		}
	}
	
	//printf("Input file: %s\n", json_file_char);
	printf("Reading JSON File... ");

	using json = nlohmann::json;
	// read a JSON file
	ifstream ifs(json_path_char);
	//json je;
	//ifs >> je;
	
	//json je = json::parse(ifs);
	// Set parameter ignore_comments to true in the parse function to ignore // or /* */ comments. 
	// Comments will then be treated as whitespace.
	// If skip_comments is set to true, the comments are skipped during parsing
	json je = json::parse(ifs,
					/* callback */ nullptr,
					/* allow exceptions */ true,
					/* skip_comments */ true);

	// Access the values
	// Types of simulations and output
	wallType = je.at("flags").value("wall_type", 1);
	femOn = je.at("flags").value("fem_MESH", false);
	forcedOn = je.at("flags").value("forced_MESH", false);
	vtuType = je.at("flags").at("output_VTU").value("type", 0);
	freeSurfWall = je.at("flags").at("output_VTU").value("only_freeSurface", false);
	outputPnd = je.at("flags").at("output_VTU").value("pnd", false);
	outputNeigh = je.at("flags").at("output_VTU").value("neigh", false);
	outputDeviation = je.at("flags").at("output_VTU").value("deviation", false);
	outputConcentration = je.at("flags").at("output_VTU").value("concentration", false);
	outputAuxiliar = je.at("flags").at("output_VTU").value("auxiliar", false);
	outputNonNewtonian = je.at("flags").at("output_VTU").value("non_newtonian", false);
	txtPress = je.at("flags").at("output_TXT").value("press", false);
	txtForce = je.at("flags").at("output_TXT").value("force", false);
	// Paths of input and output files/folders
	gridFilename = je.at("pathNames").value("particle_grid_file", "oops");
	meshRigidFilename = je.at("pathNames").value("mesh_rigid_file", "oops");
	meshDeformableFilename = je.at("pathNames").value("mesh_deformable_file", "oops");
	meshForcedFilename = je.at("pathNames").value("mesh_forced_file", "oops");
	vtuOutputFoldername = je.at("pathNames").value("vtu_output_folder", "oops");
	forceTxtFilename = je.at("pathNames").value("forceTxt_file", "oops");
	pressTxtFilename = je.at("pathNames").value("pressTxt_file", "oops");
	// Geometry dimension limits
	domainMinX = je.at("domain").at("min").value("x", 0.0);
	domainMinY = je.at("domain").at("min").value("y", 0.0);
	domainMinZ = je.at("domain").at("min").value("z", 0.0);
	domainMaxX = je.at("domain").at("max").value("x", 0.0);
	domainMaxY = je.at("domain").at("max").value("y", 0.0);
	domainMaxZ = je.at("domain").at("max").value("z", 0.0);
	// Physical parameters
	densityFluid = je.at("physical").value("fluid_density", 1000.0);
	densityWall = je.at("physical").value("wall_density", 1000.0);
	KNM_VS1 = je.at("physical").value("kinematic_visc", 0.000001);
	gravityX = je.at("physical").at("gravity").value("x", 0.0);
	gravityY = je.at("physical").at("gravity").value("y", 0.0);
	gravityZ = je.at("physical").at("gravity").value("z", -9.8);
	// Rheological parameters
	KNM_VS2 = je.at("physical").at("rheological").value("kinematic_visc_phase_2", 0.000001);
	DNS_FL1 = je.at("physical").at("rheological").value("fluid_density_phase_1", 1000.0);
	DNS_FL2 = je.at("physical").at("rheological").value("fluid_density_phase_2", 1540.0);
	DNS_SDT = je.at("physical").at("rheological").value("sediment_density", 1540.0);
	fluidType = je.at("physical").at("rheological").value("fluid_type", 0);
	N = je.at("physical").at("rheological").value("power_law_index", 1.2);
	MEU0 = je.at("physical").at("rheological").value("consistency_index", 0.03);
	PHI = je.at("physical").at("rheological").at("phi").value("lower", 0.541);
	PHI_WAL = je.at("physical").at("rheological").at("phi").value("wall", 0.541);
	PHI_BED = je.at("physical").at("rheological").at("phi").value("bed", 0.541);
	PHI_2 = je.at("physical").at("rheological").at("phi").value("second", 0.6);
	cohes = je.at("physical").at("rheological").value("cohes_coeff", 0.0);
	Fraction_method = je.at("physical").at("rheological").value("fraction_method", 2);
	//visc_max = je.at("physical").at("rheological").value("visc_max", 20);
	DG = je.at("physical").at("rheological").value("grain_size", 0.0035);
	I0 = je.at("physical").at("rheological").value("I0", 0.75);
	mm = je.at("physical").at("rheological").value("mm", 100.0);
	stress_calc_method = je.at("physical").at("rheological").value("stress_calc_method", 1);
	visc_itr_num = je.at("physical").at("rheological").at("viscosity").value("iter_num", 1);
	visc_error = je.at("physical").at("rheological").at("viscosity").value("error", 0.0);
	visc_ave = je.at("physical").at("rheological").at("viscosity").value("average", 0.0);
	Cd = je.at("physical").at("rheological").value("drag_coeff", 0.47);
	VF_min = je.at("physical").at("rheological").at("volume_fraction").value("min", 0.25);
	VF_max = je.at("physical").at("rheological").at("volume_fraction").value("max", 0.65);
	// Numerical parameters
	dim = je.at("numerical").value("dimension", 3.0);
	partDist = je.at("numerical").value("particle_dist", 0.01);
	timeStep = je.at("numerical").value("time_step", 0.0005);
	timeSimulation = je.at("numerical").value("final_time", 1.0);
	iterOutput = je.at("numerical").value("iter_output", 80);
	cflNumber = je.at("numerical").value("CFL_number", 0.2);
	weightType = je.at("numerical").value("weight_type", 0);
	slipCondition = je.at("numerical").value("slip_condition", 0);
	reS = je.at("numerical").at("effective_radius").value("small", 2.1);
	reL = je.at("numerical").at("effective_radius").value("large", 2.1);
	gradientType = je.at("numerical").at("gradient").value("type", 3);
	gradientCorrection = je.at("numerical").at("gradient").value("correction", false);
	relaxPress = je.at("numerical").at("gradient").value("relax_fact", 1.0);
	mpsType = je.at("numerical").value("mps_type", 1);
	soundSpeed = je.at("numerical").at("explicit_mps").at("equation_state").value("speed_sound", 15.0);
	gamma = je.at("numerical").at("explicit_mps").at("equation_state").value("gamma", 7.0);
	solverType = je.at("numerical").at("semi_implicit_mps").value("solver_type", 0);
	alphaCompressibility = je.at("numerical").at("semi_implicit_mps").at("weak_compressibility").value("alpha", 0.000001);
	relaxPND = je.at("numerical").at("semi_implicit_mps").at("source_term").value("relax_pnd", 0.001);
	shiftingType = je.at("numerical").at("particle_shifting").value("type", 2);
	dri = je.at("numerical").at("particle_shifting").value("DRI", 0.01);
	coefA = je.at("numerical").at("particle_shifting").value("coef_A", 2.0);
	machNumber = je.at("numerical").at("particle_shifting").value("mach_number", 0.1);
	VEL_A = je.at("numerical").at("particle_shifting").value("adj_vel_A", 0.1);
	pndType = je.at("numerical").at("pnd").value("type", 0);
	diffusiveCoef = je.at("numerical").at("pnd").value("diffusive_coeff", 0.35);
	repulsiveForceType = je.at("numerical").at("wall_repulsive_force").value("type", 2);
	reRepulsiveForce = je.at("numerical").at("wall_repulsive_force").value("re", 0.5);
	expectMaxVelocity = je.at("numerical").at("wall_repulsive_force").value("maxVel", 6.0);
	repForceCoefMitsume = je.at("numerical").at("wall_repulsive_force").at("coefficient").value("Mitsume", 40000000.0);
	repForceCoefLennardJones = je.at("numerical").at("wall_repulsive_force").at("coefficient").value("Lennard-Jones", 2.0);
	repForceCoefMonaghanKajtar = je.at("numerical").at("wall_repulsive_force").at("coefficient").value("Monaghan-Kajtar", 1.0);
	EPS_RE = je.at("numerical").at("wall_repulsive_force").value("eps_re", 0.01);
	freeSurfType = je.at("numerical").at("free_surface_threshold").value("type", 0);
	pndThreshold = je.at("numerical").at("free_surface_threshold").value("pnd", 0.98);
	neighThreshold = je.at("numerical").at("free_surface_threshold").value("neigh", 0.85);
	npcdThreshold = je.at("numerical").at("free_surface_threshold").value("NPCD", 0.20);
	thetaThreshold = je.at("numerical").at("free_surface_threshold").value("ARC", 45.0);
	normThreshold = je.at("numerical").at("free_surface_threshold").value("normal", 0.1);
	collisionType = je.at("numerical").at("particle_collision").value("type", 0);
	collisionRatio = je.at("numerical").at("particle_collision").value("ratio", 0.20);
	distLimitRatio = je.at("numerical").at("particle_collision").value("dist_limit_ratio", 0.85);
	lambdaCollision = je.at("numerical").at("particle_collision").value("lambda", 0.20);
	ghost = je.at("numerical").at("particle_type").value("ghost", -1);
	fluid = je.at("numerical").at("particle_type").value("fluid",  0);
	wall = je.at("numerical").at("particle_type").value("wall",   1);
	dummyWall = je.at("numerical").at("particle_type").value("dummyWall",   2);
	surface = je.at("numerical").at("boundary_type").value("free_surface", 1);
	inner = je.at("numerical").at("boundary_type").value("inner", 0);
	other = je.at("numerical").at("boundary_type").value("other", -1);
	numPartTypes = 3;

	printf("OK\n");

	printf("Reading GRID File... ");
	readMpsParticleFile(gridFilename);
	printf("OK\n");

	// Extend domain
	domainMinX = domainMinX - partDist*3;
	domainMinY = domainMinY - partDist*3;
	domainMinZ = domainMinZ - partDist*3;
	domainMaxX = domainMaxX + partDist*3;
	domainMaxY = domainMaxY + partDist*3;
	domainMaxZ = domainMaxZ + partDist*3;
	if((int)dim == 2) {	
		domainMinZ = 0.0;
		domainMaxZ = 0.0;
	}

	// Number of meshs
	numOfRigidMesh = 0;	numOfDeformableMesh = 0;	numOfForcedMesh = 0;
	if(wallType == 1) numOfRigidMesh = 1;
	if(femOn == true) numOfDeformableMesh = 1;
	if(forcedOn == true) numOfForcedMesh = 1;
	numOfMeshs = numOfRigidMesh + numOfDeformableMesh + numOfForcedMesh;

	// // Print the values
	// cout << "INPUT FILE .JSON" << endl;
	// cout << "Number of Meshs: " << numOfMeshs << " | ";
	// cout << "WallType:" << wallType << " | ";
	// cout << "GiraffeOn: " << femOn << " | ";
	// cout << "forcedOn: " << forcedOn << " | ";
	// cout << "vtuType: " << vtuType << " | ";
	// cout << "freeSurfWall: " << freeSurfWall << " | ";
	// cout << "outputPnd: " << outputPnd << " | ";
	// cout << "outputNeigh: " << outputNeigh << " | ";
	// cout << "outputDeviation: " << outputDeviation << " | ";
	// cout << "outputConcentration: " << outputConcentration << " | ";
	// cout << "outputAuxiliar: " << outputAuxiliar << " | ";
	// cout << "outputNonNewtonian: " << outputNonNewtonian << " | ";
	// cout << "txtPress: " << txtPress << " | ";
	// cout << "txtForce: " << txtForce << endl;
	// cout << "gridFilename: " << gridFilename << endl;
	// cout << "meshRigidFilename: " << meshRigidFilename << endl;
	// cout << "meshDeformableFilename: " << meshDeformableFilename << endl;
	// cout << "meshForcedFilename: " << meshForcedFilename << endl;
	// cout << "vtuOutputFoldername: " << vtuOutputFoldername << endl;
	// cout << "forceTxtFilename: " << forceTxtFilename << endl;
	// cout << "pressTxtFilename: " << pressTxtFilename << endl;
	// cout << "domainMin: " << domainMinX << ": " << domainMinY << ": " << domainMinZ << endl;
	// cout << "domainMax: " << domainMaxX << ": " << domainMaxY << ": " << domainMaxZ << endl;
	// cout << "gravity: " << gravityX << ": " << gravityY << ": " << gravityZ << endl;
	// cout << "densityFluid: " << densityFluid << " | ";
	// cout << "densityWall: " << densityWall << " | ";
	// cout << "KNM_VS1-2: " << KNM_VS1 << ": " << KNM_VS2 << " | ";
	// cout << "DNS_FL1-2-DST: " << DNS_FL1 << ": " << DNS_FL2 << ": " << DNS_SDT << endl;
	// cout << "fluidType: " << fluidType << " | ";
	// cout << "N: " << N << " | ";
	// cout << "MEU0: " << MEU0 << " | ";
	// cout << "PHI-FL-WAL-BED-2: " << PHI << ": " << PHI_WAL << ": " << PHI_BED << ": " << PHI_2 << endl;
	// cout << "cohes: " << cohes << " | ";
	// cout << "Fraction_method: " << Fraction_method << " | ";
	// cout << "DG: " << DG << " | ";
	// cout << "I0: " << I0 << " | ";
	// cout << "mm: " << mm << " | ";
	// cout << "stress_calc_method: " << stress_calc_method << " | ";
	// cout << "visc_itr_num-error-ave: " << visc_itr_num << ": " << visc_error << ": " << visc_ave << endl;
	// cout << "Cd: " << Cd << " | ";
	// cout << "VFminmax: " << VF_min << ": " << VF_max << endl;
	// cout << "dim: " << dim << " | ";
	// cout << "lo: " << partDist << " | ";
	// cout << "dt: " << timeStep << " | ";
	// cout << "tf: " << timeSimulation << " | ";
	// cout << "itO: " << iterOutput << " | ";
	// cout << "CFL: " << cflNumber << endl;
	// cout << "mpsType: " << mpsType << " | ";
	// cout << "weightType: " << weightType << " | ";
	// cout << "slip: " << slipCondition << " | ";
	// cout << "reSL: " << reS << ": " << reL << endl;
	// cout << "gradientType: " << gradientType << " | ";
	// cout << "gradientCorrection: " << gradientCorrection << " | ";
	// cout << "relaxPress: " << relaxPress << " | ";
	// cout << "soundSpeed: " << soundSpeed << " | ";
	// cout << "solverType: " << solverType << endl;
	// cout << "alphaCompressibility: " << alphaCompressibility << " | ";
	// cout << "relaxPND: " << relaxPND << endl;
	// cout << "shiftingType: " << shiftingType << " | ";
	// cout << "dri: " << dri << " | ";
	// cout << "coefA: " << coefA << " | ";
	// cout << "machNumber: " << machNumber << " | ";
	// cout << "VEL_A: " << VEL_A << " | ";
	// cout << "pndType: " << pndType << endl;
	// cout << "diffusiveCoef: " << diffusiveCoef << " | ";
	// cout << "repulsiveForceType: " << repulsiveForceType << " | ";
	// cout << "reRepulsiveForce: " << reRepulsiveForce << " | ";
	// cout << "expectMaxVelocity: " << expectMaxVelocity << " | ";
	// cout << "repForceCoefMitsume" << repForceCoefMitsume << " | ";
	// cout << "repForceCoefLennardJones: " << repForceCoefLennardJones << " | ";
	// cout << "repForceCoefMonaghanKajtar: " << repForceCoefMonaghanKajtar << endl;
	// cout << "EPS_RE: " << EPS_RE << " | ";
	// cout << "freeSurfType: " << freeSurfType << " | ";
	// cout << "pndThreshold: " << pndThreshold << " | ";
	// cout << "neighThreshold: " << neighThreshold << " | ";
	// cout << "npcdThreshold: " << npcdThreshold << " | ";
	// cout << "thetaThreshold: " << thetaThreshold << " | ";
	// cout << "normThreshold: " << normThreshold << " | ";
	// cout << "collisionType: " << collisionType << " | ";
	// cout << "collisionRatio: " << collisionRatio << " | ";
	// cout << "distLimitRatio: " << distLimitRatio << " | ";
	// cout << "lambdaCollision: " << lambdaCollision << endl;
	// cout << "ghost: " << ghost << " | ";
	// cout << "fluid: " << fluid << " | ";
	// cout << "wall: " << wall << " | ";
	// cout << "surface: " << surface << " | ";
	// cout << "inner: " << inner << " | ";
	// cout << "other: " << other << " | ";
	// cout << "numParticles: " << numPartTypes << endl;
		
	// // cout << endl;
	// Close .json file
	fclose(js);

}

// Read data from file .grid to class MpsParticle
void MpsParticle::readMpsParticleFile(const std::string& grid_file) {
	char *grid_file_char = new char[grid_file.length()+1];
	strcpy (grid_file_char, grid_file.c_str());
	// grid_file_char now contains a c-string copy of grid_file

	fp = fopen(grid_file_char, "r");
	if(fp == NULL) perror ("Error opening grid file");

	int zeroZero;
	fscanf(fp,"%d",&zeroZero);
	fscanf(fp,"%d",&numParticles);									// Read number of particles
	// printf("Number of particles: %d\n",numParticles);

	// Memory allocation
	// Scalars
	particleType = (int*)malloc(sizeof(int)*numParticles);			// Particle type
	particleBC = (int*)malloc(sizeof(int)*numParticles);			// BC particle type
	numNeigh = (int*)malloc(sizeof(int)*numParticles);				// Number of neighbors

	press = (double*)malloc(sizeof(double)*numParticles);			// Particle pressure
	pressAverage = (double*)malloc(sizeof(double)*numParticles);	// Time averaged particle pressure
	pndi = (double*)malloc(sizeof(double)*numParticles);			// PND
	pndki = (double*)malloc(sizeof(double)*numParticles);			// PND step k
	pndski = (double*)malloc(sizeof(double)*numParticles);			// Mean fluid neighbor PND step k
	pndSmall = (double*)malloc(sizeof(double)*numParticles);		// PND small = sum(wij)
	npcdDeviation2 = (double*)malloc(sizeof(double)*numParticles);	// NPCD deviation modulus
	concentration = (double*)malloc(sizeof(double)*numParticles);	// Concentration
	velDivergence = (double*)malloc(sizeof(double)*numParticles);	// Divergence of velocity
	diffusiveTerm = (double*)malloc(sizeof(double)*numParticles);	// Diffusive term
	
	Dns = (double*)malloc(sizeof(double)*numPartTypes);				// Density
	invDns = (double*)malloc(sizeof(double)*numPartTypes);			// Inverse of Density

	// Vectors
	acc = (double*)malloc(sizeof(double)*numParticles*3);			// Particle acceleration
	accStar = (double*)malloc(sizeof(double)*numParticles*3);		// Particle acceleration due gravity and viscosity
	pos = (double*)malloc(sizeof(double)*numParticles*3);			// Particle position
	vel = (double*)malloc(sizeof(double)*numParticles*3);			// Particle velocity
	npcdDeviation = (double*)malloc(sizeof(double)*numParticles*3);			// NPCD deviation
	gradConcentration = (double*)malloc(sizeof(double)*numParticles*3);		// Gradient of concentration
	correcMatrixRow1 = (double*)malloc(sizeof(double)*numParticles*3);		// Correction matrix - Row 1
	correcMatrixRow2 = (double*)malloc(sizeof(double)*numParticles*3);		// Correction matrix - Row 2
	correcMatrixRow3 = (double*)malloc(sizeof(double)*numParticles*3);		// Correction matrix - Row 3
	normal = (double*)malloc(sizeof(double)*numParticles*3);		// Particle normal
	dvelCollision = (double*)malloc(sizeof(double)*numParticles*3);			// Variation of velocity due collision

	// Polygons
	// Scalars
	nearMeshType = (int*)malloc(sizeof(int)*numParticles);				// Type of mesh near particle
	particleNearWall = (bool*)malloc(sizeof(bool)*numParticles);		// Particle near polygon wall
	numNeighWallContribution = (int*)malloc(sizeof(int)*numParticles);	// Number of neighbors due wall

	pndWallContribution = (double*)malloc(sizeof(double)*numParticles);			// PND wall
	deviationDotPolygonNormal = (double*)malloc(sizeof(double)*numParticles);	// Deviation vector X polygonal wall
	numNeighborsSurfaceParticles = (double*)malloc(sizeof(double)*numParticles);// Number of free-surface particle neighbors
	distParticleWall2 = (double*)malloc(sizeof(double)*numParticles);			// Squared distance of particle to triangle mesh
	// Vectors
	particleAtWallPos = (double*)malloc(sizeof(double)*numParticles*3);		// Particle at wall coordinate
	mirrorParticlePos = (double*)malloc(sizeof(double)*numParticles*3);		// Mirrored particle coordinate
	wallParticleForce1 = (double*)malloc(sizeof(double)*numParticles*3);	// Wall-Particle force
	wallParticleForce2 = (double*)malloc(sizeof(double)*numParticles*3);	// Wall-Particle force
	polygonNormal = (double*)malloc(sizeof(double)*numParticles*3);		// Polygon normal
	
//	Posk = (double*)malloc(sizeof(double)*numParticles*3);		// Particle coordinates
//	Velk = (double*)malloc(sizeof(double)*numParticles*3);		// Particle velocity
//	Acv = (double*)malloc(sizeof(double)*numParticles*3);		// Part

	// Non-Newtonian
	// Scalars
	PTYPE = (int*)malloc(sizeof(int)*numParticles);				// Type of fluid

	Cv = (double*)malloc(sizeof(double)*numParticles);			// Concentration
	II = (double*)malloc(sizeof(double)*numParticles);			// Invariant
	MEU = (double*)malloc(sizeof(double)*numParticles);			// Dynamic viscosity
	MEU_Y = (double*)malloc(sizeof(double)*numParticles);		// Dynamic viscosity ??
	Inertia = (double*)malloc(sizeof(double)*numParticles);		//
	pnew = (double*)malloc(sizeof(double)*numParticles);		// New pressure
	p_rheo_new = (double*)malloc(sizeof(double)*numParticles);	//
	RHO = (double*)malloc(sizeof(double)*numParticles);			// Fluid density
	p_smooth = (double*)malloc(sizeof(double)*numParticles);	//
	VF = (double*)malloc(sizeof(double)*numParticles);			//
	S12 = (double*)malloc(sizeof(double)*numParticles);			//
	S13 = (double*)malloc(sizeof(double)*numParticles);			//
	S23 = (double*)malloc(sizeof(double)*numParticles);			//
	S11 = (double*)malloc(sizeof(double)*numParticles);			//
	S22 = (double*)malloc(sizeof(double)*numParticles);			//
	S33 = (double*)malloc(sizeof(double)*numParticles);			//

	// FSI
	// Scalars
	elementID = (int*)malloc(sizeof(int)*numParticles);			// Element ID
	// Vectors
	forceWall = (double*)malloc(sizeof(double)*numParticles*3);	// Force on wall

	// Solver PPE
	pressurePPE = Eigen::VectorXd::Zero(numParticles);
	sourceTerm = Eigen::VectorXd::Zero(numParticles);

	// Set values from .grid file
	for(int i=0; i<numParticles; i++) {
		int a[2];
		double b[8];

		// Uncomment here to read .prof file
		//fscanf(fp," %d %d %lf %lf %lf %lf %lf %lf %lf %lf",&a[0],&a[1],&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7]);
		// Uncomment here to read .grid file
		a[0] = 0;
		fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf",&a[1],&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7]);
		particleType[i]=a[1];
		pos[i*3]=b[0];	pos[i*3+1]=b[1];	pos[i*3+2]=b[2];
		vel[i*3]=b[3];	vel[i*3+1]=b[4];	vel[i*3+2]=b[5];
		//	|				|					|
		//	i*3 = x			i*3+1 = y			i*3+2 = z
		press[i]=b[6];	pressAverage[i]=b[7];

//		printf("X: %d %lf %lf %lf %lf %lf %lf %lf %lf\n",particleType[i],pos[3*i],pos[3*i+1],pos[3*i+2],vel[3*i],vel[3*i+1],vel[3*i+2],press[i],pressAverage[i]);
//		Posk[i*3]=b[0];	Posk[i*3+1]=b[1];	Posk[i*3+2]=b[2];
//		Velk[i*3]=b[3];	Velk[i*3+1]=b[4];	Velk[i*3+2]=b[5];
	}
	// Close .grid file
	fclose(fp);
	// Verify if particle is out of domain
	checkParticleOutDomain();
	// Set vectors to zero
	for(int i=0;i<numParticles*3;i++) {
		acc[i]=0.0;accStar[i]=0.0;npcdDeviation[i]=0.0;gradConcentration[i]=0.0;
		correcMatrixRow1[i]=0.0;correcMatrixRow2[i]=0.0;correcMatrixRow3[i]=0.0;normal[i]=0.0;dvelCollision[i]=0.0;//Acv[i]=0.0;
		particleAtWallPos[i]=0.0;mirrorParticlePos[i]=0.0;wallParticleForce1[i]=0.0;wallParticleForce2[i]=0.0;polygonNormal[i]=0.0;
		forceWall[i]=0.0;
	}

	// Set scalars to zero or infinity(10e8)
	for(int i=0; i<numParticles; i++) {
		particleBC[i]=0;numNeigh[i]=0;numNeighWallContribution[i]=0;elementID[i]=0;
		particleNearWall[i]=false;
		nearMeshType[i]=meshType::FIXED;

		pndi[i]=0.0;pndki[i]=0.0;pndski[i]=0.0;pndSmall[i]=0.0;npcdDeviation2[i]=0.0;concentration[i]=0.0;
		velDivergence[i]=0.0;diffusiveTerm[i]=0.0;pndWallContribution[i]=0.0;deviationDotPolygonNormal[i]=0.0;
		numNeighborsSurfaceParticles[i]=0.0;Cv[i]=0.0;II[i]=0.0;MEU_Y[i]=0.0;Inertia[i]=0.0;pnew[i]=0.0;
		p_rheo_new[i]=0.0;p_smooth[i]=0.0;VF[i]=0.0;S12[i]=0.0;S13[i]=0.0;S23[i]=0.0;S11[i]=0.0;S22[i]=0.0;S33[i]=0.0;

		distParticleWall2[i]=10e8*partDist;
	}
	// Assign type and density
	for(int i=0; i<numParticles; i++) {
		/*
		// Assign type and density
		if(pos[i*3+2] <= 0.3) {
			PTYPE[i]=2;
			RHO[i] = DNS_FL2;
			// CHANGED Only at the first time step
			MEU[i] = KNM_VS2 * DNS_FL2;
		}
		else {
			PTYPE[i]=1;
			RHO[i] = DNS_FL1;
			// CHANGED Only at the first time step
			MEU[i] = KNM_VS1 * DNS_FL1;
		}
		*/

		if(fluidType == viscType::NEWTONIAN) {
			RHO[i] = DNS_FL1;
			PTYPE[i] = 1;
			MEU[i] = KNM_VS1 * DNS_FL1;
		}
		// Multiphase simulations - Granular Fluid
		if(fluidType == viscType::NON_NEWTONIAN) {
			// Assign type and density
			if(particleType[i] == 1) {
				particleType[i] = 0;
				PTYPE[i] = 2;
				RHO[i] = DNS_FL2;
				// CHANGED Only at the first time step
				MEU[i] = KNM_VS2 * DNS_FL2;
			}
			else {
				particleType[i] = 0;
				PTYPE[i] = 1;
				RHO[i] = DNS_FL1;
				// CHANGED Only at the first time step
				MEU[i] = KNM_VS1 * DNS_FL1;
			}
		}
	}
}

// Allocation of buckets
// Murotani et al., 2015. Performance improvements of differential operators code for MPS method on GPU.
void MpsParticle::allocateBuckets() {
	reS = partDist*reS;								// Influence radius small
	reL = partDist*reL;								// Influence radius large
	reS2 = reS*reS;									// Influence radius small to square
	reL2 = reL*reL;									// Influence radius large to square
	EPS_RE = EPS_RE*reS2/4.0;
	reRepulsiveForce = partDist*reRepulsiveForce;	// Influence radius for repulsive force
	bucketSide = reL*(1.0+cflNumber);				// Length of one bucket side
	bucketSide2 = bucketSide*bucketSide;
	invBucketSide = 1.0/bucketSide;

	numBucketsX = (int)((domainMaxX - domainMinX)*invBucketSide) + 3;		// Number of buckets in the x direction in the analysis domain
	numBucketsY = (int)((domainMaxY - domainMinY)*invBucketSide) + 3;		// Number of buckets in the y direction in the analysis domain
	numBucketsZ = (int)((domainMaxZ - domainMinZ)*invBucketSide) + 3;		// Number of buckets in the z direction in the analysis domain
	if((int)dim == 2) {	numBucketsZ = 1; }
	numBucketsXY = numBucketsX*numBucketsY;
	numBucketsXYZ = numBucketsX*numBucketsY*numBucketsZ;					// Number of buckets in analysis area
	
	firstParticleInBucket = 	(int*)malloc(sizeof(int) * numBucketsXYZ);	// First particle number stored in the bucket
	lastParticleInBucket = 		(int*)malloc(sizeof(int) * numBucketsXYZ);	// Last particle number stored in the bucket
	nextParticleInSameBucket  = (int*)malloc(sizeof(int) * numParticles);	// Next particle number in the same bucket
}

// Set parameters
void MpsParticle::setParameters() {
	pndSmallZero = pndLargeZero = pndGradientZero = lambdaZero = numNeighZero = 0.0;
	int lmin = ceil(reL/partDist) + 1;
	int lmax = ceil(reL/partDist) + 2;
	int flag2D = 0;
	int flag3D = 1;
	if((int)dim == 2) {
		flag2D = 1;
		flag3D = 0;
	}
	for(int ix= -lmin; ix<lmax; ix++) {
	for(int iy= -lmin; iy<lmax; iy++) {
	for(int iz= -lmin*flag3D; iz<lmax*flag3D+flag2D; iz++) {
		double x = partDist* (double)ix;
		double y = partDist* (double)iy;
		double z = partDist* (double)iz;
		double dst2 = x*x+y*y+z*z;
		if(dst2 <= reL2) {
			if(dst2 == 0.0) continue;
			double dst = sqrt(dst2);
			pndLargeZero += weight(dst, reL, weightType);			// Initial particle number density (large)
			lambdaZero += dst2 * weight(dst, reL, weightType);
			numNeighZero += 1;										// Initial number of neighbors
			if(dst2 <= reS2) {
				pndSmallZero += weight(dst, reS, weightType);		// Initial particle number density (small)
				pndGradientZero += weightGradient(dst, reS, weightType);	// Initial particle number density (gradient operator)
			}
		}
	}}}
	lambdaZero = lambdaZero/pndLargeZero;							// Coefficient λ of Laplacian model
	coeffViscosity = 2.0*KNM_VS1*dim/(pndLargeZero*lambdaZero);		// Coefficient used to calculate viscosity term
	coeffViscMultiphase = 2.0*dim/(pndLargeZero*lambdaZero);		// Coefficient used to calculate viscosity term Multiphase
	coeffPressEMPS = soundSpeed*soundSpeed/pndSmallZero;			// Coefficient used to calculate pressure E-MPS
	coeffPressGrad = -dim/pndGradientZero;							// Coefficient used to calculate pressure gradient term
	coeffPressWCMPS = soundSpeed*soundSpeed;						// Coefficient used to calculate pressure WC-MPS
	coeffShifting1 = dri*partDist/pndSmallZero;						// Coefficient used to adjust velocity type 1
	coeffShifting2 = coefA*partDist*partDist*cflNumber*machNumber;	// Coefficient used to adjust velocity type 2
	coeffPPE = 2.0*dim/(pndLargeZero*lambdaZero);					// Coefficient used to PPE
	coeffPPESource = relaxPND/(timeStep*timeStep*pndSmallZero);		// Coefficient used to PPE source term
	Dns[partType::FLUID]=densityFluid;			Dns[partType::WALL]=densityWall;
	invDns[partType::FLUID]=1.0/densityFluid;	invDns[partType::WALL]=1.0/densityWall;
	invPartDist = 1.0/partDist;
	distCollisionLimit = partDist*distLimitRatio;					// A distance that does not allow further access between particles
	distCollisionLimit2 = distCollisionLimit*distCollisionLimit;
	restitutionCollision = 1.0 + collisionRatio;
	numOfIterations = 0;											// Number of iterations
	fileNumber = 0;													// File number
	timeCurrent = 0.0;												// Simulation time
	velMax = 0.0;													// Maximum flow velocity
	CFLcurrent = cflNumber;											// Current Courant number
	betaPnd = pndThreshold*pndSmallZero;							// Surface cte PND
	betaNeigh = neighThreshold*numNeighZero;						// Surface cte Neighbors
	delta2 = npcdThreshold*npcdThreshold*partDist*partDist;			// Surface cte NPCD 
	thetaArc = thetaThreshold/180.0*3.14159265;						// Surface cte theta ARC
	hThreshold2 = 1.33*1.33*partDist*partDist;						// Surface cte radius ARC
	dstThreshold2 = 2.0*hThreshold2;								// Surface cte radius ARC
	normThreshold2 = normThreshold*normThreshold;					// Surface cte Normal
	
	//cout << "lo: " << partDist << " m, dt: " << timeStep << " s, PND0Small: " << pndSmallZero << " PND0Large: " << pndLargeZero << " PND0Grad: " << pndGradientZero << " lambda: " << lambdaZero << std::endl;
	//cout << "bPnd: " << betaPnd << "betaNeigh: " << betaNeigh << endl;
}

// Set initial PND and number of neighbors
void MpsParticle::setInitialPndNumberOfNeigh() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double wSum = 0.0;
		numNeigh[i] = 0;
		npcdDeviation[i*3] = npcdDeviation[i*3+1] = npcdDeviation[i*3+2] = 0;
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
					if(j != i) {
						numNeigh[i] += 1;
						if(dstij2 < reS2) {
							double dst = sqrt(dstij2);
							double wS = weight(dst, reS, weightType);
							pndi[i] += wS;

							npcdDeviation[i*3  ] += v0ij*wS*invPartDist;
							npcdDeviation[i*3+1] += v1ij*wS*invPartDist;
							npcdDeviation[i*3+2] += v2ij*wS*invPartDist;
							wSum += wS;
						}
					}
				}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// Add PND due wall polygon
		pndi[i] += pndWallContribution[i];
		if(particleType[i] == wall)
			pndi[i] = pndSmallZero;
		pndSmall[i] = pndi[i];
		pndki[i] = pndi[i];
		// Add Number of neighbors due wall polygon
		numNeigh[i] += numNeighWallContribution[i];

		if(wSum > 1.0e-8) {
			//npcdDeviation[i*3  ] /= pndSmall[i];
			//npcdDeviation[i*3+1] /= pndSmall[i];
			//npcdDeviation[i*3+2] /= pndSmall[i];
			npcdDeviation[i*3  ] /= wSum;
			npcdDeviation[i*3+1] /= wSum;
			npcdDeviation[i*3+2] /= wSum;
		}

		npcdDeviation2[i] = npcdDeviation[i*3]*npcdDeviation[i*3] + npcdDeviation[i*3+1]*npcdDeviation[i*3+1] +
			npcdDeviation[i*3+2]*npcdDeviation[i*3+2];

		//deviationDotPolygonNormal[i] = npcdDeviation[i*3]*polygonNormal[i*3]+npcdDeviation[i*3+1]*polygonNormal[i*3+1]+npcdDeviation[i*3+2]*polygonNormal[i*3+2];
		if(npcdDeviation[i*3]*polygonNormal[i*3]+npcdDeviation[i*3+1]*polygonNormal[i*3+1]+npcdDeviation[i*3+2]*polygonNormal[i*3+2] < 0.0)
			deviationDotPolygonNormal[i] = 1;
		else
			deviationDotPolygonNormal[i] = -1;
	}
}


////////////////////////////////////////////////////////////
// Functions called during the simulation (main loop)
////////////////////////////////////////////////////////////

// Verify if particle is out of domain
void MpsParticle::checkParticleOutDomain() {
	for(int i=0; i<numParticles; i++) {
		if(	pos[i*3  ]>domainMaxX || pos[i*3  ]<domainMinX ||
			pos[i*3+1]>domainMaxY || pos[i*3+1]<domainMinY ||
			pos[i*3+2]>domainMaxZ || pos[i*3+2]<domainMinZ) {

			// ID of last particle
			int iLastParticle = numParticles - 1;
			// Move the data from "Last Particle" to i-th ghost particle

			// Scalars
			particleType[i]=particleType[iLastParticle];
			particleBC[i]=particleBC[iLastParticle];
			numNeigh[i]=numNeigh[iLastParticle];
			press[i]=press[iLastParticle];
			pressAverage[i]=pressAverage[iLastParticle];
			pndi[i]=pndi[iLastParticle];
			pndki[i]=pndki[iLastParticle];
			pndski[i]=pndski[iLastParticle];
			pndSmall[i]=pndSmall[iLastParticle];
			npcdDeviation2[i]=npcdDeviation2[iLastParticle];
			concentration[i]=concentration[iLastParticle];
			velDivergence[i]=velDivergence[iLastParticle];
			diffusiveTerm[i]=diffusiveTerm[iLastParticle];
			nearMeshType[i]=nearMeshType[iLastParticle];
			particleNearWall[i]=particleNearWall[iLastParticle];
			numNeighWallContribution[i]=numNeighWallContribution[iLastParticle];
			pndWallContribution[i]=pndWallContribution[iLastParticle];
			deviationDotPolygonNormal[i]=deviationDotPolygonNormal[iLastParticle];
			numNeighborsSurfaceParticles[i]=numNeighborsSurfaceParticles[iLastParticle];
			distParticleWall2[i]=distParticleWall2[iLastParticle];
			PTYPE[i]=PTYPE[iLastParticle];
			Cv[i]=Cv[iLastParticle];
			II[i]=II[iLastParticle];
			MEU[i]=MEU[iLastParticle];
			MEU_Y[i]=MEU_Y[iLastParticle];
			Inertia[i]=Inertia[iLastParticle];
			pnew[i]=pnew[iLastParticle];
			p_rheo_new[i]=p_rheo_new[iLastParticle];
			RHO[i]=RHO[iLastParticle];
			p_smooth[i]=p_smooth[iLastParticle];
			VF[i]=VF[iLastParticle];
			S12[i]=S12[iLastParticle];
			S13[i]=S13[iLastParticle];
			S23[i]=S23[iLastParticle];
			S11[i]=S11[iLastParticle];
			S22[i]=S22[iLastParticle];
			S33[i]=S33[iLastParticle];
			elementID[i]=elementID[iLastParticle];

			pressurePPE(i)=pressurePPE(iLastParticle);
			sourceTerm(i)=sourceTerm(iLastParticle);
			
			// Vectors
			for (int j = 0; j < 3; j++)
			{
				acc[i*3+j]=acc[iLastParticle*3+j];
				accStar[i*3+j]=accStar[iLastParticle*3+j];
				pos[i*3+j]=pos[iLastParticle*3+j];
				vel[i*3+j]=vel[iLastParticle*3+j];
				npcdDeviation[i*3+j]=npcdDeviation[iLastParticle*3+j];
				gradConcentration[i*3+j]=gradConcentration[iLastParticle*3+j];
				correcMatrixRow1[i*3+j]=correcMatrixRow1[iLastParticle*3+j];
				correcMatrixRow2[i*3+j]=correcMatrixRow2[iLastParticle*3+j];
				correcMatrixRow3[i*3+j]=correcMatrixRow3[iLastParticle*3+j];
				normal[i*3+j]=normal[iLastParticle*3+j];
				dvelCollision[i*3+j]=dvelCollision[iLastParticle*3+j];
				particleAtWallPos[i*3+j]=particleAtWallPos[iLastParticle*3+j];
				mirrorParticlePos[i*3+j]=mirrorParticlePos[iLastParticle*3+j];
				wallParticleForce1[i*3+j]=wallParticleForce1[iLastParticle*3+j];
				wallParticleForce2[i*3+j]=wallParticleForce2[iLastParticle*3+j];
				polygonNormal[i*3+j]=polygonNormal[iLastParticle*3+j];
				forceWall[i*3+j]=forceWall[iLastParticle*3+j];
				//Posk[i*3+j]=Posk[iLastParticle*3+j];
				//Velk[i*3+j]=Velk[iLastParticle*3+j];
				//Acv[i*3+j]=Acv[iLastParticle*3+j];
			}

			// Update some data of "Last Particle"
			particleType[iLastParticle]=ghost;
			particleBC[iLastParticle]=other;
			particleNearWall[iLastParticle]=false;
			nearMeshType[iLastParticle]=meshType::FIXED;
			distParticleWall2[iLastParticle]=10e8*partDist;
			// Set zero to position of lastParticle
			for (int j = 0; j < 3; j++){
				pos[iLastParticle*3+j]=0.0;
				//Posk[iLastParticle*3+j]=0.0;
			}
			
			// Decrease number of particles
			numParticles--;
		}
	}
}

// Update particle ID's in buckets
// Update particle ID's in buckets
void MpsParticle::updateBuckets2D() {
	for(int i=0; i<numBucketsXY ;i++) 	{	firstParticleInBucket[i] = -1;	}
	for(int i=0; i<numBucketsXY ;i++) 	{	lastParticleInBucket[i] = -1;	}
	for(int i=0; i<numParticles ;i++) 	{	nextParticleInSameBucket[i] = -1;	}
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == ghost) continue;
		int ix = (int)((pos[i*3  ] - domainMinX)*invBucketSide) + 1;
		int iy = (int)((pos[i*3+1] - domainMinY)*invBucketSide) + 1;
		int ib = iy*numBucketsX + ix;
		int j = lastParticleInBucket[ib];
		lastParticleInBucket[ib] = i;
		if(j == -1) {	firstParticleInBucket[ib] = i;	}
		else 		{	nextParticleInSameBucket[j] = i;}
	}
}

void MpsParticle::updateBuckets3D() {
	for(int i=0; i<numBucketsXYZ ;i++) 	{	firstParticleInBucket[i] = -1;	}
	for(int i=0; i<numBucketsXYZ ;i++) 	{	lastParticleInBucket[i] = -1;	}
	for(int i=0; i<numParticles ;i++) 	{	nextParticleInSameBucket[i] = -1;	}
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == ghost) continue;
		int ix = (int)((pos[i*3  ] - domainMinX)*invBucketSide) + 1;
		int iy = (int)((pos[i*3+1] - domainMinY)*invBucketSide) + 1;
		int iz = (int)((pos[i*3+2] - domainMinZ)*invBucketSide) + 1;
		int ib = iz*numBucketsXY + iy*numBucketsX + ix;
		int j = lastParticleInBucket[ib];
		lastParticleInBucket[ib] = i;
		if(j == -1) {	firstParticleInBucket[ib] = i;	}
		else 		{	nextParticleInSameBucket[j] = i;}
	}
}

// Acceleration due Laplacian of velocity and gravity
void MpsParticle::calcViscosityGravity() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
//		if(particleType[i] == fluid) {
		double meu_i = MEU[i];
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
					if(j != i) {
						double dst = sqrt(dstij2);
						double wL = weight(dst, reL, weightType);
						if ((meu_i + MEU[j]) > 1.0e-8)
							NEU = 2 * meu_i * MEU[j] / (meu_i + MEU[j]);
						else
							NEU = 0.0;
						//NEU = KNM_VS2 * DNS_FL2;
						if(PTYPE[i] == 1) NEU = NEU/DNS_FL1;
						else NEU = NEU/DNS_FL2;
	//					NEU = NEU/RHO[i];
						//if((NEUt[i] + NEUt[j]) > 0) NEU = NEU + (2 * NEUt[i] * RHO[j] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j])) / RHO[i];
						// Original
	//					accX +=(vel[j*3  ]-velXi)*w;
	//					accY +=(vel[j*3+1]-velYi)*w;
	//					accZ +=(vel[j*3+2]-velZi)*w;
						// Modified
						accX +=(vel[j*3  ]-velXi)*wL*NEU;
						accY +=(vel[j*3+1]-velYi)*wL*NEU;
						accZ +=(vel[j*3+2]-velZi)*wL*NEU;
					}
				}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// Original
//		acc[i*3  ]=accX*coeffViscosity + gravityX;
//		acc[i*3+1]=accY*coeffViscosity + gravityY;
//		acc[i*3+2]=accZ*coeffViscosity + gravityZ;
		// Modified
		//if(timeCurrent > 0.3) {
		// coeffViscMultiphase = 2.0*dim/(pndLargeZero*lambdaZero);
		acc[i*3  ] = coeffViscMultiphase*accX + gravityX;
		acc[i*3+1] = coeffViscMultiphase*accY + gravityY;
		acc[i*3+2] = coeffViscMultiphase*accZ + gravityZ;
		//}		
		accStar[i*3  ] = acc[i*3  ];
		accStar[i*3+1] = acc[i*3+1];
		accStar[i*3+2] = acc[i*3+2];
	}
}

// Prediction of pressure gradient
void MpsParticle::predictionPressGradient() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++){
//		if(particleType[i] == fluid) {
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double Pi = press[i];			double ni = pndi[i];		double pressMin = Pi;
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		if(gradientType == 0 || gradientType == 2) {
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
						if(j != i) {
							if(pressMin > press[j]) pressMin = press[j];
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
		}
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
					if(j != i) {
						double dst = sqrt(dstij2);
						double wS = weightGradient(dst, reS, weightType);
						if(gradientType == 0)
							wS *= (press[j] - pressMin)/dstij2;
						else if(gradientType == 1)
							wS *= (press[j] + Pi)/dstij2;
						else if(gradientType == 2)
							wS *= (press[j] + Pi - 2.0*pressMin)/dstij2;
						else if(gradientType == 3) {
							double nj = pndi[j];
							if(ni > 1.0e-8 && nj > 1.0e-8)
								wS *= (ni*press[j]/nj + nj*Pi/ni)/dstij2;
						}
						accX += v0ij*wS;	accY += v1ij*wS;	accZ += v2ij*wS;
					}
				}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// coeffPressGrad is a negative cte (-dim/noGrad)
		// Original
//		acc[i*3  ]+=(1.0-relaxPress)*accX*invDns[partType::FLUID]*coeffPressGrad;
//		acc[i*3+1]+=(1.0-relaxPress)*accY*invDns[partType::FLUID]*coeffPressGrad;
//		acc[i*3+2]+=(1.0-relaxPress)*accZ*invDns[partType::FLUID]*coeffPressGrad;
		// Modified
		acc[i*3  ]+=(1.0-relaxPress)*accX*coeffPressGrad/RHO[i];
		acc[i*3+1]+=(1.0-relaxPress)*accY*coeffPressGrad/RHO[i];
		acc[i*3+2]+=(1.0-relaxPress)*accZ*coeffPressGrad/RHO[i];
	}
}

// Prediction of pressure gradient (Polygon wall)
void MpsParticle::predictionWallPressGradient() {
	// Maximum velocity is the minimum of the computed and expected maximum velocities
	double maxVelocity = min(velMax, expectMaxVelocity);
	double velMax2 = maxVelocity*maxVelocity;
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
		//particleNearWall[i]=true; // Only to show particles near polygon
		if(particleType[i] == fluid && particleNearWall[i] == true) {
			double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
			double Pi = press[i];			double ni = pndi[i];		double pressMin = Pi;
			
			// Wall gradient Mitsume`s model
		    double Rref_i[9], normaliw[3], normaliwSqrt;
		    // normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		    normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
		    normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);
		    if(normaliwSqrt > 1.0e-8) {
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
			Rai[0] = Rref_i[0]*accStar[i*3] + Rref_i[1]*accStar[i*3+1] + Rref_i[2]*accStar[i*3+2];
			Rai[1] = Rref_i[3]*accStar[i*3] + Rref_i[4]*accStar[i*3+1] + Rref_i[5]*accStar[i*3+2];
			Rai[2] = Rref_i[6]*accStar[i*3] + Rref_i[7]*accStar[i*3+1] + Rref_i[8]*accStar[i*3+2];

			int ix, iy, iz;
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			if(gradientType == 0 || gradientType == 2) {
				for(int jz=minZ;jz<=maxZ;jz++) {
				for(int jy=iy-1;jy<=iy+1;jy++) {
				for(int jx=ix-1;jx<=ix+1;jx++) {
					int jb = jz*numBucketsXY + jy*numBucketsX + jx;
					int j = firstParticleInBucket[jb];
					if(j == -1) continue;
					while(true) {
						double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
						// Particle distance r_ij = Xj - Xi_temporary_position
						sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
						// Mirror particle distance r_imj = Xj - Xim_temporary_position
						sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

						// If j is inside the neighborhood of i and 
						// is not at the same side of im (avoid real j in the virtual neihborhood)
						if(dstij2 < reS2 && dstij2 < dstimj2) {
						if(j != i) {
							if(pressMin > press[j]) pressMin = press[j];
						}}
						j = nextParticleInSameBucket[j];
						if(j == -1) break;
					}
				}}}
			}
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
						if(j != i) {

							double dst = sqrt(dstimj2);
							double wS = weightGradient(dst, reS, weightType);
							
							// Taylor pressure Pj
							double Pj;
		//					Pj = Pi + RHO[i]*(Rai[0]*v0 + Rai[1]*v1 + Rai[2]*v2);
							Pj = press[j];

							if(gradientType == 0)
								wS *= (Pj - pressMin)/dstimj2;//(press[j] - pressMin)/dstimj2;
							else if(gradientType == 1)
								wS *= (Pj + Pi)/dstimj2;//(press[j] + Pi)/dstimj2;
							else if(gradientType == 2)
								wS *= (Pj + Pi - 2.0*pressMin)/dstimj2;//(press[j] + Pi - 2.0*pressMin)/dstimj2;
							else if(gradientType == 3) {
								double nj = pndi[j];
								if(ni > 1.0e-8 && nj > 1.0e-8)
									wS *= (ni*Pj/nj + nj*Pi/ni)/dstimj2;
							}
						accX += v0imj*wS;	accY += v1imj*wS;	accZ += v2imj*wS;
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
		  	
			if(dstimi2 < reS2) {
				double dst = sqrt(dstimi2);
				double wS = weightGradient(dst, reS, weightType);

				// Taylor pressure Pj
				double Pj;
	//			Pj = Pi + RHO[i]*(Rai[0]*v0 + Rai[1]*v1 + Rai[2]*v2);
				Pj = Pi;

				if(gradientType == 0)
					wS *= (Pj - pressMin)/dstimi2;//(Pi - pressMin)/dstimi2
				else if(gradientType == 1)
					wS *= (Pj + Pi)/dstimi2;//(Pi + Pi)/dstimi2
				else if(gradientType == 2)
					wS *= (Pj + Pi - 2.0*pressMin)/dstimi2;//(Pi + Pi - 2.0*pressMin)/dstimi2
				else if(gradientType == 3) {
					double nj = pndi[i];
					if(ni > 1.0e-8 && nj > 1.0e-8)
						wS *= (ni*Pj/nj + nj*Pi/ni)/dstimi2;
				}
				accX += v0imi*wS;	accY += v1imi*wS;	accZ += v2imi*wS;
		  	}

			// Repulsive force
			double rpsForce[3];
			rpsForce[0]=rpsForce[1]=rpsForce[2] = 0.0;

			if(repulsiveForceType == repForceType::HARADA) {
				// Parallel analysis system for free-surface flow using MPS method with explicitly represented polygon wall boundary model
				// https://doi.org/10.1007/s40571-019-00269-6
				if(normaliwSqrt < reRepulsiveForce && normaliwSqrt > 1.0e-8) {
					double wijRep = RHO[i]/(timeStep*timeStep)*(reRepulsiveForce-normaliwSqrt);
					rpsForce[0] = - wijRep*normaliw[0];
					rpsForce[1] = - wijRep*normaliw[1];
					rpsForce[2] = - wijRep*normaliw[2];
					//if(i == 6) {
						//printf("x:%lf y:%lf z:%lf x:%lf y:%lf z:%lf\n", pos[i*3],pos[i*3+1],pos[i*3+2],mirrorParticlePos[i*3],mirrorParticlePos[i*3+1],mirrorParticlePos[i*3+2]);
						//printf("nx:%lf ny: %lf nz: %lf nN:%lf\n", normaliw[0],normaliw[1],normaliw[2],normaliwSqrt);
						//printf("Fx:%lf Fy:%lf Fz:%lf\n", rpsForce[0],rpsForce[1],rpsForce[2]);
					//}
				}
			}
			else if(repulsiveForceType == repForceType::MITSUME) {
				// Explicitly represented polygon wall boundary model for the explicit MPS method
				// https://doi.org/10.1007/s40571-015-0037-8
				if(normaliwSqrt < reRepulsiveForce && normaliwSqrt > 1.0e-8) {
					double wijRep = repForceCoefMitsume*weightGradient(normaliwSqrt, reRepulsiveForce, weightType);
					rpsForce[0] = - wijRep*normaliw[0];
					rpsForce[1] = - wijRep*normaliw[1];
					rpsForce[2] = - wijRep*normaliw[2];
					//if(i == 6) {
						//printf("x:%lf y:%lf z:%lf x:%lf y:%lf z:%lf\n", pos[i*3],pos[i*3+1],pos[i*3+2],mirrorParticlePos[i*3],mirrorParticlePos[i*3+1],mirrorParticlePos[i*3+2]);
						//printf("nx:%lf ny: %lf nz: %lf nN:%lf\n", normaliw[0],normaliw[1],normaliw[2],normaliwSqrt);
						//printf("Fx:%lf Fy:%lf Fz:%lf\n", rpsForce[0],rpsForce[1],rpsForce[2]);
					//}
				}
			}
			else if(repulsiveForceType == repForceType::LENNARD_JONES) {
				// Simulating Free Surface Flows with SPH
				// https://doi.org/10.1006/jcph.1994.1034
				if(normaliwSqrt < reRepulsiveForce && normaliwSqrt > 1.0e-8) {
					double R1 = (reRepulsiveForce/normaliwSqrt)*(reRepulsiveForce/normaliwSqrt);
					double R2 = R1*R1;
					double wijRep = (repForceCoefLennardJones*velMax2/normaliwSqrt)*(R2-R1)*RHO[i];
					rpsForce[0] = - wijRep*normaliw[0];
					rpsForce[1] = - wijRep*normaliw[1];
					rpsForce[2] = - wijRep*normaliw[2];
				}
			}
		  	else {
				// SPH particle boundary forces for arbitrary boundaries 
				// https://doi.org/10.1016/j.cpc.2009.05.008
				if(normaliwSqrt < reRepulsiveForce && normaliwSqrt > 1.0e-8) {
					double W1 = (1.0+3.0/2.0*normaliwSqrt/(reRepulsiveForce));
					double W2 = (1.0-normaliwSqrt/(reRepulsiveForce))*(1.0-normaliwSqrt/(reRepulsiveForce))*(1.0-normaliwSqrt/(reRepulsiveForce));
					double wijRep = (repForceCoefMonaghanKajtar*velMax2/(normaliwSqrt - 0.0*partDist))*(1.0/8.0)*(W1)*(W2)*RHO[i];
					rpsForce[0] = - wijRep*normaliw[0];
					rpsForce[1] = - wijRep*normaliw[1];
					rpsForce[2] = - wijRep*normaliw[2];
				}
		  	}
			// coeffPressGrad is a negative cte (-dim/noGrad)
			// Original
	//		acc[i*3  ] += ((1.0-relaxPress)*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*coeffPressGrad - rpsForce[0])*invDns[partType::FLUID];
	//		acc[i*3+1] += ((1.0-relaxPress)*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*coeffPressGrad - rpsForce[1])*invDns[partType::FLUID];
	//		acc[i*3+2] += ((1.0-relaxPress)*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*coeffPressGrad - rpsForce[2])*invDns[partType::FLUID];
			// Modified
			acc[i*3  ] += ((1.0-relaxPress)*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*coeffPressGrad - rpsForce[0])/RHO[i];
			acc[i*3+1] += ((1.0-relaxPress)*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*coeffPressGrad - rpsForce[1])/RHO[i];
			acc[i*3+2] += ((1.0-relaxPress)*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*coeffPressGrad - rpsForce[2])/RHO[i];

			//Fwall[i*3  ] =  (Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*invDns[partType::FLUID]*coeffPressGrad - rpsForce[0]*invDns[partType::FLUID];
			//Fwall[i*3+1] =  (Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*invDns[partType::FLUID]*coeffPressGrad - rpsForce[1]*invDns[partType::FLUID];
			//Fwall[i*3+2] =  (Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*invDns[partType::FLUID]*coeffPressGrad - rpsForce[2]*invDns[partType::FLUID];
		}
	}
}

// Update velocity and position
void MpsParticle::updateVelocityPosition1st() {
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
			vel[i*3  ] += acc[i*3  ]*timeStep;	vel[i*3+1] += acc[i*3+1]*timeStep;	vel[i*3+2] += acc[i*3+2]*timeStep;
			//if(particleType[i] == fluid) {
			pos[i*3  ] += vel[i*3  ]*timeStep;	pos[i*3+1] += vel[i*3+1]*timeStep;	pos[i*3+2] += vel[i*3+2]*timeStep;
			//}
		}
		acc[i*3]=acc[i*3+1]=acc[i*3+2]=0.0;
		dvelCollision[i*3]=dvelCollision[i*3+1]=dvelCollision[i*3+2]=0.0;
		wallParticleForce1[i*3]=wallParticleForce1[i*3+1]=wallParticleForce1[i*3+2]=0.0;
		wallParticleForce2[i*3]=wallParticleForce2[i*3+1]=wallParticleForce2[i*3+2]=0.0;
		npcdDeviation[i*3]=npcdDeviation[i*3+1]=npcdDeviation[i*3+2]=0.0;
		numNeighWallContribution[i]=0;
		particleNearWall[i]=false;
		// Set squared distance of particle to triangle mesh to ~infinite
		distParticleWall2[i] = 10e8*partDist;
		numNeighborsSurfaceParticles[i]=0.0;

		if(wallType == boundaryWallType::POLYGON) {
			// Set mirrored particle to ~infinite if wall particles are used
			mirrorParticlePos[i*3  ] = 10e8*partDist; mirrorParticlePos[i*3+1] = 10e8*partDist; mirrorParticlePos[i*3+2] = 10e8*partDist;
		}
	}
}

// Check collisions between particles
// Step-by-step improvement of MPS method in simulating violent free-surface motions and impact-loads
// https://doi.org/10.1016/j.cma.2010.12.001
void MpsParticle::checkParticleCollisions() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
	//		double mi = Dns[partType::FLUID];
			double mi;
			if(PTYPE[i] == 1) mi = DNS_FL1;
			else mi = DNS_FL2;
			double posXi = pos[i*3  ];double posYi = pos[i*3+1];double posZi = pos[i*3+2];
			double velXi = vel[i*3  ];double velYi = vel[i*3+1];double velZi = vel[i*3+2];
			double dVelXi = 0.0;double dVelYi = 0.0;double dVelZi = 0.0;
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
			
			int ix, iy, iz;
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < distCollisionLimit2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
						if(j != i) {
							double fDT = (velXi-vel[j*3  ])*v0ij+(velYi-vel[j*3+1])*v1ij+(velZi-vel[j*3+2])*v2ij;
							if(fDT > 0.0) {
								double mj;
								if(particleType[j]==fluid)
								{
									if(PTYPE[j] == 1) mj = DNS_FL1;
									else mj = DNS_FL2;
								}
								else
								{
									mj = Dns[partType::WALL];
								}

								fDT *= restitutionCollision*mj/(mi+mj)/dstij2;
								if(particleType[j]==fluid)
								{
									dVelXi -= v0ij*fDT;		dVelYi -= v1ij*fDT;		dVelZi -= v2ij*fDT;
								}
								else
								{
									dVelXi -= 2.0*v0ij*fDT;	dVelYi -= 2.0*v1ij*fDT;	dVelZi -= 2.0*v2ij*fDT;
								}
							}
							/*
							double fDT = (vel[j*3  ]-velXi)*v0+(vel[j*3+1]-velYi)*v1+(vel[j*3+2]-velZi)*v2;
							double mj;
							if(particleType[j]==fluid)
								mj = Dns[partType::FLUID];
							else
								mj = Dns[partType::WALL];
							fDT *= restitutionCollision*mj/(mi+mj)/dst2;
							velXi2 += v0*fDT;		vecYi2 += v1*fDT;		vecZi2 += v2*fDT;
							*/
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}

			dvelCollision[i*3  ]=dVelXi;	dvelCollision[i*3+1]=dVelYi;	dvelCollision[i*3+2]=dVelZi;
			//accStar[i*3  ]=vel[i*3  ]+dVelXi;	accStar[i*3+1]=vel[i*3+1]+dVelYi;	accStar[i*3+2]=vel[i*3+2]+dVelZi;
		}
	}
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
			// CHANGED !!!
			//pos[i*3  ]+=(acc[i*3  ]-vel[i*3  ])*timeStep; pos[i*3+1]+=(acc[i*3+1]-vel[i*3+1])*timeStep; pos[i*3+2]+=(acc[i*3+2]-vel[i*3+2])*timeStep;
			vel[i*3  ]+=dvelCollision[i*3  ];	vel[i*3+1]+=dvelCollision[i*3+1];	vel[i*3+2]+=dvelCollision[i*3+2];
			
			//Velk[i*3  ]=vel[i*3  ];	Velk[i*3+1]=vel[i*3+1];	Velk[i*3+2]=vel[i*3+2];
			//pos[i*3  ]=Posk[i*3  ]+vel[i*3  ]*timeStep; pos[i*3+1]=Posk[i*3+1]+vel[i*3+1]*timeStep; pos[i*3+2]=Posk[i*3+2]+vel[i*3+2]*timeStep;
		}
		dvelCollision[i*3  ]=0.0;	dvelCollision[i*3+1]=0.0;	dvelCollision[i*3+2]=0.0;
	}
}

// Check collisions between particles (Dynamic Particle Collision)
// Enhanced weakly-compressible MPS method for violent free-surface flows: Role of particle regularization techniques
// https://doi.org/10.1016/j.jcp.2021.110202
void MpsParticle::checkDynamicParticleCollisions() {
	
	double Wij5 = 0.5*0.5*0.5*0.5*3.0;
	//double Wij5 = 0.5*0.5;
	double pmax = 0.0;
	
	// Compute maximum pressure on the walls
#pragma omp parallel
{
	double local_pmax = 0.0;
#pragma omp for
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == wall) {
			local_pmax = max(local_pmax, press[i]);
		}
	}
#pragma omp critical
	{
		if (local_pmax > pmax)
			pmax = local_pmax;
	}
}
	// Compute collision and repulsive terms and the dynamic coefficients
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
	//		double mi = Dns[partType::FLUID];
			double mi;
			if(PTYPE[i] == 1) mi = DNS_FL1;
			else mi = DNS_FL2;

			double posXi = pos[i*3  ];double posYi = pos[i*3+1];double posZi = pos[i*3+2];
			double velXi = vel[i*3  ];double velYi = vel[i*3+1];double velZi = vel[i*3+2];
			double dVelXi = 0.0;double dVelYi = 0.0;double dVelZi = 0.0;
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];

			int ix, iy, iz;
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < partDist*partDist && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
						if(j != i) {
							double mj;
							if(particleType[j]==fluid)
							{
								if(PTYPE[j] == 1) mj = DNS_FL1;
								else mj = DNS_FL2;
							}
							else
							{
								mj = Dns[partType::WALL];
							}
							double dst = sqrt(dstij2);
							
							// inter-particle distance
							double Wij = 0.0;
							if(dst > 1.0e-8 && dst < partDist)
							{
								double w1 = 1.0 - dst*invPartDist;
								double w2 = 4.0*dst*invPartDist + 1;
								Wij = w1*w1*w1*w1*w2;
								//Wij = w1*w1;
							}
							double chi = sqrt(Wij/Wij5);
							double kappa = 0.0;
							if(dst < 0.5*partDist)
							{
								kappa = 1.0;
							}
							else if(dst >= 0.5*partDist && dst < partDist)
							{
								kappa = chi;
							}

							double fDT = (velXi-vel[j*3  ])*v0ij+(velYi-vel[j*3+1])*v1ij+(velZi-vel[j*3+2])*v2ij;
							if(fDT > 0.0) {
								fDT *= kappa*2.0*mj/(mi+mj)/dstij2;
								//fDT *= restitutionCollision*mj/(mi+mj)/dstij2;
								if(particleType[j]==fluid)
								{
									dVelXi -= v0ij*fDT;		dVelYi -= v1ij*fDT;		dVelZi -= v2ij*fDT;
								}
								else //if(particleBC[i] == surface)
								{
									dVelXi -= 2.0*v0ij*fDT;		dVelYi -= 2.0*v1ij*fDT;		dVelZi -= 2.0*v2ij*fDT;
								}
							}
							else
							{
								// Dynamic background pressure
								double pmax = 2.0/3.0*DNS_FL1*gravityY*0.3;
								double pmin = DNS_FL1*gravityY*partDist;
								double ptil = max(min(lambdaCollision*fabs(press[i]+press[j]), lambdaCollision*pmax), pmin);
								double pb = ptil*chi;
								double rep = timeStep/mi*chi*pb/dstij2;
								if(particleType[j]==fluid)
								{
									dVelXi -= v0ij*rep;		dVelYi -= v1ij*rep;		dVelZi -= v2ij*rep;
								}
								else //if(particleBC[i] == surface)
								{
									dVelXi -= 2.0*v0ij*rep;		dVelYi -= 2.0*v1ij*rep;		dVelZi -= 2.0*v2ij*rep;
								}
							}
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}

			dvelCollision[i*3  ]=dVelXi;	dvelCollision[i*3+1]=dVelYi;	dvelCollision[i*3+2]=dVelZi;
		}
	}
	// Update velocity and position
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
			
			vel[i*3  ]+=dvelCollision[i*3  ];	vel[i*3+1]+=dvelCollision[i*3+1];	vel[i*3+2]+=dvelCollision[i*3+2];
			pos[i*3  ]+=dvelCollision[i*3  ]*timeStep; pos[i*3+1]+=dvelCollision[i*3+1]*timeStep; pos[i*3+2]+=dvelCollision[i*3+2]*timeStep;
			/*
			double drNew[3], drMod, drMin, duNew[3];
			duNew[0] = dvelCollision[i*3  ];
			duNew[1] = dvelCollision[i*3+1];
			duNew[2] = dvelCollision[i*3+2];
			drNew[0] = dvelCollision[i*3  ]*timeStep;
			drNew[1] = dvelCollision[i*3+1]*timeStep;
			drNew[2] = dvelCollision[i*3+2]*timeStep;
			drMod = (drNew[0]*drNew[0] + drNew[1]*drNew[1] + drNew[2]*drNew[2]);
			drMin = min(0.1*partDist, drMod);
			if(drMin > 1.0e-8)
			{
				pos[i*3  ]+=drMin*drNew[0]/drMod;
				pos[i*3+1]+=drMin*drNew[1]/drMod;
				pos[i*3+2]+=drMin*drNew[2]/drMod;
				vel[i*3  ]+=drMin*duNew[0]/drMod;
				vel[i*3+1]+=drMin*duNew[1]/drMod;
				vel[i*3+2]+=drMin*duNew[2]/drMod;
			}
			*/
		}
		dvelCollision[i*3  ]=0.0;	dvelCollision[i*3+1]=0.0;	dvelCollision[i*3+2]=0.0;
	}
}

// Set force on wall to zero
//void MpsParticle::WallZeroForce_omp(int nNodes, int nSolids, solid_fem * &solid) {
void MpsParticle::setWallForceZero(const int nNodes, double *nodeforceX, double *nodeforceY, double *nodeforceZ) {
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
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
}

// Free-surface particles. NPCD (Polygon wall)
void MpsParticle::calcWallNPCD() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid && particleNearWall[i] == true) {
			double npcdDeviationXi = 0.0;	double npcdDeviationYi = 0.0;	double npcdDeviationZi = 0.0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
			// Wall gradient Mitsume`s model
			double Rref_i[9], normaliw[3], normaliwSqrt;
			// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
			normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
			normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);

			if(normaliwSqrt > 1.0e-8) {
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
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
						if(j != i) {
							double dst = sqrt(dstimj2);
							double wS = weightGradient(dst, reS, weightType);
							npcdDeviationXi += v0imj*wS;
							npcdDeviationYi += v1imj*wS;
							npcdDeviationZi += v2imj*wS;
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}

			}}}
			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
			
			if(dstimi2 < reS2) {
				double dst = sqrt(dstimi2);
				double wS = weightGradient(dst, reS, weightType);
				npcdDeviationXi += v0imi*wS;
				npcdDeviationYi += v1imi*wS;
				npcdDeviationZi += v2imi*wS;
			}
			npcdDeviation[i*3  ] += Rref_i[0]*npcdDeviationXi + Rref_i[1]*npcdDeviationYi + Rref_i[2]*npcdDeviationZi;
			npcdDeviation[i*3+1] += Rref_i[3]*npcdDeviationXi + Rref_i[4]*npcdDeviationYi + Rref_i[5]*npcdDeviationZi;
			npcdDeviation[i*3+2] += Rref_i[6]*npcdDeviationXi + Rref_i[7]*npcdDeviationYi + Rref_i[8]*npcdDeviationZi;
		}
	}
}

// Compute PND, number of neighbors and NPCD
void MpsParticle::calcPndnNeighNPCD() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double ni = 0.0; double wSum = 0.0;
		numNeigh[i] = 0;
		// Add Number of neighbors due Wall polygon
		numNeigh[i] += numNeighWallContribution[i];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
					if(j != i) {
						numNeigh[i] += 1;
						//double dst = sqrt(dst2);
						//double wL = weight(dst, reL*invPartDist, weightType);
						//npcdDeviation[i*3  ] += v0*wL*invPartDist;
						//npcdDeviation[i*3+1] += v1*wL*invPartDist;
						//npcdDeviation[i*3+2] += v2*wL*invPartDist;
						//wSum += wL;
						if(dstij2 < reS2) {
							double dst = sqrt(dstij2);
							double wS = weight(dst, reS, weightType);
							ni += wS;
							dst = dst*invPartDist;
							wS = weight(dst, reS*invPartDist, weightType);
							npcdDeviation[i*3  ] += v0ij*wS*invPartDist;
							npcdDeviation[i*3+1] += v1ij*wS*invPartDist;
							npcdDeviation[i*3+2] += v2ij*wS*invPartDist;
							wSum += wS;
						}
					}
				}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

//		double mi;
//		if(PTYPE[i] == 1) mi = DNS_FL1;
//		else mi = DNS_FL2;
		//if(particleType[i]==fluid)
		//	mi = Dns[partType::FLUID];
		//else
		//	mi = Dns[partType::WALL];

		if(pndType == calcPNDType::SUM_WIJ || pndType == calcPNDType::MEAN_SUM_WIJ) {
//		if(pndType == calcPNDType::SUM_WIJ) {
			// PND at initial of step k
			pndki[i] = pndi[i];
			// New PND due particles and Wall polygon
			pndi[i] = ni + pndWallContribution[i];
		}

//		if(particleType[i] == wall) {
			// PND due particles and Wall polygon
//			pndi[i] = ni + pndWallContribution[i];
//			if(pndi[i] < pndSmallZero)
//				pndi[i] = pndSmallZero;

//				pndi[i] = pndSmallZero*pow((press[i]*gamma/(mi*coeffPressWCMPS)+1),gamma);
//		}

		// Add PND due Wall polygon
		pndSmall[i] = ni + pndWallContribution[i];
		// Prevent pndSmall[i] = 0
//		if(numNeigh[i]>1) {
		if(wSum > 1.0e-8) {
			//npcdDeviation[i*3  ] /= pndSmall[i];
			//npcdDeviation[i*3+1] /= pndSmall[i];
			//npcdDeviation[i*3+2] /= pndSmall[i];
			npcdDeviation[i*3  ] /= wSum;
			npcdDeviation[i*3+1] /= wSum;
			npcdDeviation[i*3+2] /= wSum;
		}

		npcdDeviation2[i] = npcdDeviation[i*3]*npcdDeviation[i*3]+npcdDeviation[i*3+1]*npcdDeviation[i*3+1]+npcdDeviation[i*3+2]*npcdDeviation[i*3+2];

		//deviationDotPolygonNormal[i] = npcdDeviation[i*3]*polygonNormal[i*3]+npcdDeviation[i*3+1]*polygonNormal[i*3+1]+npcdDeviation[i*3+2]*polygonNormal[i*3+2];
		if(npcdDeviation[i*3]*polygonNormal[i*3]+npcdDeviation[i*3+1]*polygonNormal[i*3+1]+npcdDeviation[i*3+2]*polygonNormal[i*3+2]< 0.0)
			deviationDotPolygonNormal[i] = 1;
		else
			deviationDotPolygonNormal[i] = -1;
		
		// First check based on particle number density
//		if(pndSmall[i] < pndThreshold*pndSmallZero)
//			particleBC[i] = surface;
//		else
//			particleBC[i] = inner;

		// Boundary particle verification based on relative distance and weight (NPCD)
//		if(particleBC[i] == surface) {
//			if(numNeigh[i] > 4 && npcdDeviation2[i] < delta2)
//			{
//				particleBC[i] = inner;
				//printf(" inner %d \n", i);
//			}
//		}

//		if(pndSmall[i] < pndThreshold*pndSmallZero && numNeigh[i] < neighThreshold*numNeighZero)
//			particleBC[i] = surface;
//		else
//			particleBC[i] = inner;
	}
}

// Diffusive term of density/PND (pndType = calcPNDType::DIFFUSIVE)
// An enhanced weakly-compressible MPS method for free-surface flows
// https://doi.org/10.1016/j.cma.2019.112771
void MpsParticle::calcPndDiffusiveTerm() {
	// coeffViscMultiphase = 2.0*dim/(pndLargeZero*lambdaZero);
	double C1 = diffusiveCoef*timeStep*soundSpeed*soundSpeed*coeffViscMultiphase/(pndLargeZero);
	double C2 = diffusiveCoef*partDist*soundSpeed*coeffViscMultiphase/(pndLargeZero);
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
//	if(particleType[i] == fluid) {
		double Di = 0.0; double DivV = 0.0; double flagDi = 1.0; double pndAux = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		//double M1[3][3];
		//for(int im=0; im<3; im++)
		//{
		//	for(int jm=0; jm<3; jm++)
		//	{
		//		M1[im][jm] = 0.0;
		//	}
		//}

		double ni = pndi[i];
		if(ni < 1.0e-8) continue;
		double Pi = press[i];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
					if(j != i) {
	//				if(j != i && particleType[j] == fluid) {
						double dst = sqrt(dstij2);
						double wL = weight(dst, reL, weightType);
						double nj = pndi[j];
						if(particleType[i] == fluid && particleType[j] == fluid) {
						//if(particleType[i] == inner) {
	//					if(particleType[i] == fluid && particleType[j] == fluid && particleBC[i] == inner) {
							// coeffViscMultiphase = 2.0*dim/(pndLargeZero*lambdaZero);
	//						double pgh = RHO[i]*(gravityX*v0+gravityY*v1+gravityZ*v2);
	//						double CB = RHO[i]*soundSpeed*soundSpeed/gamma;
	//						double nijH = pndSmallZero*(pow((pgh+1.0)/CB,1.0/gamma)-1.0);
	//						double nijH = pndSmallZero*(pow((pgh)/CB+1.0,1.0/gamma)-1.0);
	//
						//	pow( ( pgh + 1.0 ) / CB - 1.0 , 1.0  )

							//double CB = soundSpeed*soundSpeed*RHO[i];
							//double nijH = pndSmallZero*((PijH+1.0)/CB-1);
	//						if(isnan(nijH) == 0)
	//							Di += C1*(nj - ni - nijH)*wL;
	//						else
								////////Di += C1*(nj - ni)*wL;
							//if(isnan(nijH) == 1)
							//if(i == 200)
							//	printf(" pgh %e CB %e ni %e nj %e nijH %e res %e \n", pgh, CB, ni, nj, nijH, pndSmallZero*(pow((pgh+1)/CB,1/gamma)-1));
							// PND
	//						Di += C1*(nj-ni)*wL;
							//Di += C2*(nj-ni)*wL;
							// Pressure
							// Delta Voronoi smoothed particle hydrodynamics, δ-VSPH
							// https://doi.org/10.1016/j.jcp.2019.109000
							double pgh = -RHO[i]*(gravityX*v0ij+gravityY*v1ij+gravityZ*v2ij);
							double Pj = press[j];
							Di += timeStep/RHO[i]*coeffViscMultiphase*(Pj-Pi-pgh)*wL;
							//Di += (partDist/soundSpeed)/RHO[i]*coeffViscMultiphase*(Pj-Pi+pgh)*wL;
						}
						//else
						//	flagDi = 0.0;
						if(dstij2 < reS2) {
							double vijx = vel[j*3  ]-velXi;
							double vijy = vel[j*3+1]-velYi;
							double vijz = vel[j*3+2]-velZi;
							double wS = weight(dst, reS, weightType);
							if(ni > 1.0e-8)
							{
								DivV += (dim/pndSmallZero)*(nj/ni)*(vijx*v0ij+vijy*v1ij+vijz*v2ij)*wS/dstij2;
							}

	//						M1[0][0] += (dim/pndSmallZero)*(nj/ni)*(v0*vijx)*wS/dst2; M1[0][1] += (dim/pndSmallZero)*(nj/ni)*(v0*vijy)*wS/dst2; M1[0][2] += (dim/pndSmallZero)*(nj/ni)*(v0*vijz)*wS/dst2;
	//						M1[1][0] += (dim/pndSmallZero)*(nj/ni)*(v1*vijx)*wS/dst2; M1[1][1] += (dim/pndSmallZero)*(nj/ni)*(v1*vijy)*wS/dst2; M1[1][2] += (dim/pndSmallZero)*(nj/ni)*(v1*vijz)*wS/dst2;
	//						M1[2][0] += (dim/pndSmallZero)*(nj/ni)*(v2*vijx)*wS/dst2; M1[2][1] += (dim/pndSmallZero)*(nj/ni)*(v2*vijy)*wS/dst2; M1[2][2] += (dim/pndSmallZero)*(nj/ni)*(v2*vijz)*wS/dst2;

	//						if(particleType[i] == wall)
	//						if(particleType[i] == wall || particleBC[i] == surface)
	//							pndAux += wS;
						}
					}
				}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}


//		DivV = 	correcMatrixRow1[i*3  ]*M1[0][0] + correcMatrixRow1[i*3+1]*M1[1][0] + correcMatrixRow1[i*3+2]*M1[2][0] +
//				correcMatrixRow2[i*3  ]*M1[0][1] + correcMatrixRow2[i*3+1]*M1[1][1] + correcMatrixRow2[i*3+2]*M1[2][1] +
//				correcMatrixRow3[i*3  ]*M1[0][2] + correcMatrixRow3[i*3+1]*M1[1][2] + correcMatrixRow3[i*3+2]*M1[2][2];

		acc[i*3] = pndi[i]*(1.0+timeStep*(-DivV+Di*flagDi));


//		if(isnan(DivV) || isnan(Di))
//			printf(" i %d \n", i);
		velDivergence[i] = DivV;
		diffusiveTerm[i] = Di;

	//acc[i*3] = pndi[i]*(1.0+timeStep*(-(1.0-diffusiveCoef)*DivV+Di*flagDi));
	//acc[i*3] = pndSmall[i]*(1.0+timeStep*(-DivV+Di*flagDi));
	//acc[i*3] =     pndSmallZero*(1.0+timeStep*(-DivV+Di*flagDi)); // Ruim
//		if(particleType[i] == wall)
//		{
//			if(pndAux < pndSmallZero)
//				pndAux = pndSmallZero;
//			acc[i*3] = pndAux;
//		}
//		if(particleBC[i] == surface)
//		{
//			acc[i*3] = pndAux;
//		}
	}
//#pragma omp parallel for
//	for(int i=0; i<numParticles; i++) {
/////	if(particleType[i] == fluid) {
//		pndi[i] = acc[i*3];
//		acc[i*3]=0.0;
//	}
}

// Diffusive term of density/PND (Polygon wall) - Free-slip (pndType = calcPNDType::DIFFUSIVE)
void MpsParticle::calcWallSlipPndDiffusiveTerm() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
		double ni = pndi[i];
		//if(particleType[i] == fluid && ni > 1.0e-8) {
		if(particleType[i] == fluid && particleNearWall[i] == true && ni > 1.0e-8) {
	//	if(particleType[i] == fluid) {
			double DivV = 0.0;
			//double Pi = press[i];
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
	//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2

			// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
			double Rref_i[9], normaliw[3], normalMod2;
			// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
			normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
			normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];
			if(normalMod2 > 1.0e-8) {
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
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
						if(j != i) {
							double dst = sqrt(dstimj2);
							double wS = weight(dst, reS, weightType);
							double nj = pndi[j];
							double vijx = vel[j*3  ]-velMirrorXi;
							double vijy = vel[j*3+1]-velMirrorYi;
							double vijz = vel[j*3+2]-velMirrorZi;
							DivV += (dim/pndSmallZero)*(nj/ni)*(vijx*v0imj+vijy*v1imj+vijz*v2imj)*wS/dstimj2;
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}

			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
			
			if(dstimi2 < reS2) {
				double dst = sqrt(dstimi2);
				double wS = weight(dst, reS, weightType);
				double nj = pndi[i];
				double vijx = velXi-velMirrorXi;
				double vijy = velYi-velMirrorYi;
				double vijz = velZi-velMirrorZi;
				DivV += (dim/pndSmallZero)*(nj/ni)*(vijx*v0imi+vijy*v1imi+vijz*v2imi)*wS/dstimi2;
		  	}

			acc[i*3] += -pndi[i]*timeStep*DivV;
			//acc[i*3] = pndSmallZero*(1.0+timeStep*(-DivV+Di*flagDi));
		}
	}
//#pragma omp parallel for
//	for(int i=0; i<numParticles; i++) {
//	if(particleType[i] == fluid) {
//		pndi[i] += acc[i*3];
//		acc[i*3]=0.0;
//	}}
}

// Diffusive term of density/PND (Polygon wall) - No-slip (pndType = calcPNDType::DIFFUSIVE)
void MpsParticle::calcWallNoSlipPndDiffusiveTerm() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
		double ni = pndi[i];
		//if(particleType[i] == fluid && ni > 1.0e-8) {
		if(particleType[i] == fluid && particleNearWall[i] == true && ni > 1.0e-8) {
	//	if(particleType[i] == fluid) {
			double DivV = 0.0;
			//double Pi = press[i];
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
	//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2

			// Inverse matrix Rinv_i = - I
			double Rinv_i[9], Rref_i[9], normaliw[3], normalMod2;
		    // normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		    normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
		    normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];
		    if(normalMod2 > 1.0e-8) {
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

			if(nearMeshType[i] == meshType::FORCED) {
				viwall[0] = velVWall[0];
				viwall[1] = velVWall[1];
				viwall[2] = velVWall[2];
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
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
						if(j != i) {
							double dst = sqrt(dstimj2);
							double wS = weight(dst, reS, weightType);
							double nj = pndi[j];
							double vijx = -(vel[j*3  ]-velMirrorXi);
							double vijy = -(vel[j*3+1]-velMirrorYi);
							double vijz = -(vel[j*3+2]-velMirrorZi);
							// Refelected rij' = Rref_i * ri'j
	      					double v0m = (Rref_i[0]*v0imj + Rref_i[1]*v1imj + Rref_i[2]*v2imj);
							double v1m = (Rref_i[3]*v0imj + Rref_i[4]*v1imj + Rref_i[5]*v2imj);
							double v2m = (Rref_i[6]*v0imj + Rref_i[7]*v1imj + Rref_i[8]*v2imj);
							DivV += (dim/pndSmallZero)*(nj/ni)*(vijx*v0m+vijy*v1m+vijz*v2m)*wS/dstimj2;
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
					}
			}}}

			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
			
			if(dstimi2 < reS2) {
				double dst = sqrt(dstimi2);
				double wS = weight(dst, reS, weightType);
				double nj = pndi[i];
				double vijx = -(velXi-velMirrorXi);
				double vijy = -(velYi-velMirrorYi);
				double vijz = -(velZi-velMirrorZi);
				// Refelected rij' = Rref_i * ri'j
				double v0m = (Rref_i[0]*v0imi + Rref_i[1]*v1imi + Rref_i[2]*v2imi);
				double v1m = (Rref_i[3]*v0imi + Rref_i[4]*v1imi + Rref_i[5]*v2imi);
				double v2m = (Rref_i[6]*v0imi + Rref_i[7]*v1imi + Rref_i[8]*v2imi);
				DivV += (dim/pndSmallZero)*(nj/ni)*(vijx*v0m+vijy*v1m+vijz*v2m)*wS/dstimi2;
		  	}

			acc[i*3] += -pndi[i]*timeStep*DivV;
			//acc[i*3] = pndSmallZero*(1.0+timeStep*(-DivV+Di*flagDi));
		}
	}
//#pragma omp parallel for
//	for(int i=0; i<numParticles; i++) {
//	if(particleType[i] == fluid) {
//		pndi[i] += acc[i*3];
//		acc[i*3]=0.0;
//	}}
}

// Update PND (pndType = calcPNDType::DIFFUSIVE)
void MpsParticle::updatePnd() {
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
	//		if(particleType[i] == fluid)
			pndi[i] = acc[i*3];
//		else
//		{
//			pndi[i] = acc[i*3];
			/*
			double mi;
			if(PTYPE[i] == 1) 
				mi = DNS_FL1;
			else 
				mi = DNS_FL2;
			if(mpsType == calcPressType::EXPLICIT)
				pndi[i] = pndSmallZero*(press[i]/(mi*coeffPressWCMPS)+1);
			else if(mpsType == calcPressType::WEAKLY)
				pndi[i] = pndSmallZero*pow(press[i]*gamma/(mi*coeffPressWCMPS)+1,gamma);
				*/
//		}
		acc[i*3]=0.0;
	}
}

// Mean PND at wall and dummy particles (pndType = calcPNDType::DIFFUSIVE)
void MpsParticle::meanPndParticlesWallDummySurface() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
//	if(particleType[i] == wall) {
		if(particleType[i] == wall || particleBC[i] == surface) {
//	if(particleBC[i] == surface) {
			double PNDup = 0.0;
			double PNDdo = 0.0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
			
			int ix, iy, iz;
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
						if(j != i) {
							double dst = sqrt(dstij2);
							double wS = weight(dst, reS, weightType);
							PNDup += pndi[j]*wS;
							PNDdo += wS;
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			acc[i*3  ] = PNDup;
			acc[i*3+1] = PNDdo;
			//acc[i*3] = pndSmallZero*(1.0+timeStep*(-DivV+Di*flagDi));
	//	}}}
		}
	}
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
//	if(particleType[i] == wall) {
		if(particleType[i] == wall || particleBC[i] == surface) {
//	if(particleBC[i] == surface) {
		// Prevent PNDdo = 0
//		if(numNeigh[i] < 1)
			if(acc[i*3+1] < 1.0e-8)
				pndi[i] = acc[i*3];
			else
				pndi[i] = acc[i*3]/(acc[i*3+1]);
//			pndi[i] = acc[i*3]/(acc[i*3+1] + 0.01*reS2/4.0);
		}
//	}
		acc[i*3]=0.0;acc[i*3+1]=0.0;
	}
}

// Mean PND (Polygon wall) (pndType = calcPNDType::MEAN_SUM_WIJ)
void MpsParticle::meanWallPnd() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
//	if(particleType[i] == fluid) {
		if(particleNearWall[i] == true) {
			double PNDup = 0.0;
			double PNDdo = 0.0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];

			int ix, iy, iz;
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
						if(j != i) {
							double dst = sqrt(dstimj2);
							double wS = weight(dst, reS, weightType);
							PNDup += pndi[j]*wS;
							PNDdo += wS;
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}

			}}}
			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
			
			if(dstimi2 < reS2) {
				double dst = sqrt(dstimi2);
				double wS = weight(dst, reS, weightType);
				PNDup += pndi[i]*wS;
				PNDdo += wS;
			}
			acc[i*3  ] = PNDup;
			acc[i*3+1] = PNDdo;
		}
	}
}

// Mean PND (pndType = calcPNDType::MEAN_SUM_WIJ)
void MpsParticle::meanPnd() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		double PNDup = 0.0;
		double PNDdo = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
					if(j != i) {
						double dst = sqrt(dstij2);
						double wS = weight(dst, reS, weightType);
						PNDup += pndi[j]*wS;
						PNDdo += wS;
					}
				}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		acc[i*3  ] += PNDup;
		acc[i*3+1] += PNDdo;
		//acc[i*3] = pndSmallZero*(1.0+timeStep*(-DivV+Di*flagDi));
	}
#pragma omp parallel for
for(int i=0; i<numParticles; i++) {
//	if(particleType[i] == fluid) {
	// Prevent PNDdo = 0
		if(numNeigh[i] < 1) {
			pndi[i] = acc[i*3];
		}
		else {
			pndi[i] = acc[i*3]/acc[i*3+1];
		}
		acc[i*3]=0.0;acc[i*3+1]=0.0;
	}
}

// Mean Fluid PND (pndType = calcPNDType::MEAN_SUM_WIJ) // CHANGED
void MpsParticle::meanNeighFluidPnd() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		double PNDup = 0.0;
		double PNDdo = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
					if(j != i && particleType[j] == fluid) {
						double dst = sqrt(dstij2);
						double wS = weight(dst, reS, weightType);
						PNDup += pndki[j]*wS;
						PNDdo += wS;
					}
				}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		if(PNDdo > 1.0e-8) {
			pndski[i] = PNDup/PNDdo;
		}
		else {
			pndski[i] = PNDup;
		}
	}
}

// Update type of particle
void MpsParticle::updateParticleBC() {
// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
/*
	double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
	double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
	double ni = 0.0; double wSum = 0.0;
	numNeigh[i] = 0;
	// Add Number of neighbors due Wall polygon
	numNeigh[i] += numNeighWallContribution[i];
	
	int ix, iy, iz;
	bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
	int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
	for(int jz=minZ;jz<=maxZ;jz++) {
	for(int jy=iy-1;jy<=iy+1;jy++) {
	for(int jx=ix-1;jx<=ix+1;jx++) {
		int jb = jz*numBucketsXY + jy*numBucketsX + jx;
		int j = firstParticleInBucket[jb];
		if(j == -1) continue;
		while(true) {
			double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
			
			// Particle distance r_ij = Xj - Xi_temporary_position
			sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
			// Mirror particle distance r_imj = Xj - Xim_temporary_position
			sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

			// If j is inside the neighborhood of i and 
			// is not at the same side of im (avoid real j in the virtual neihborhood)
			if(dstij2 < reL2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
			if(j != i) {
				numNeigh[i] += 1;
				//double dst = sqrt(dstij2);
				//double wL = weight(dst, reL*invPartDist, weightType);
				//npcdDeviation[i*3  ] += v0ij*wL*invPartDist;
				//npcdDeviation[i*3+1] += v1ij*wL*invPartDist;
				//npcdDeviation[i*3+2] += v2ij*wL*invPartDist;
				//wSum += wL;
				if(dstij2 < reS2) {
					double dst = sqrt(dstij2);
					double wS = weight(dst, reS, weightType);
					ni += wS;
					dst = dst*invPartDist;
					wS = weight(dst, reS*invPartDist, weightType);;
					npcdDeviation[i*3  ] += v0ij*wS*invPartDist;
					npcdDeviation[i*3+1] += v1ij*wS*invPartDist;
					npcdDeviation[i*3+2] += v2ij*wS*invPartDist;
					wSum += wS;
			}}}
			j = nextParticleInSameBucket[j];
			if(j == -1) break;
		}
	}}}

	double mi;
	if(PTYPE[i] == 1) mi = DNS_FL1;
	else mi = DNS_FL2;
	//if(particleType[i]==fluid)
	//	mi = Dns[partType::FLUID];
	//else
	//	mi = Dns[partType::WALL];

	if(pndType == calcPNDType::SUM_WIJ || pndType == calcPNDType::MEAN_SUM_WIJ)
//		if(pndType == calcPNDType::SUM_WIJ)
	{
		// PND due particles and Wall polygon
		pndi[i] = ni + pndWallContribution[i];
	}

//		if(particleType[i] == wall) {
		// PND due particles and Wall polygon
//			pndi[i] = ni + pndWallContribution[i];
//			if(pndi[i] < pndSmallZero)
//				pndi[i] = pndSmallZero;

//				pndi[i] = pndSmallZero*pow((press[i]*gamma/(mi*coeffPressWCMPS)+1),gamma);
//		}

	// Add PND due Wall polygon
	pndSmall[i] = ni + pndWallContribution[i];
	// Prevent pndSmall[i] = 0
//		if(numNeigh[i] > 1) {
	if(wSum > 1.0e-8) {
		//npcdDeviation[i*3  ] /= pndSmall[i];
		//npcdDeviation[i*3+1] /= pndSmall[i];
		//npcdDeviation[i*3+2] /= pndSmall[i];
		npcdDeviation[i*3  ] /= wSum;
		npcdDeviation[i*3+1] /= wSum;
		npcdDeviation[i*3+2] /= wSum;
	}

	npcdDeviation2[i] = npcdDeviation[i*3]*npcdDeviation[i*3]+npcdDeviation[i*3+1]*npcdDeviation[i*3+1]+npcdDeviation[i*3+2]*npcdDeviation[i*3+2];

	//deviationDotPolygonNormal[i] = npcdDeviation[i*3]*polygonNormal[i*3]+npcdDeviation[i*3+1]*polygonNormal[i*3+1]+npcdDeviation[i*3+2]*polygonNormal[i*3+2];
	if(npcdDeviation[i*3]*polygonNormal[i*3]+npcdDeviation[i*3+1]*polygonNormal[i*3+1]+npcdDeviation[i*3+2]*polygonNormal[i*3+2] < 0.0)
		deviationDotPolygonNormal[i] = 1;
	else
		deviationDotPolygonNormal[i] = -1;

	// coeffPressWCMPS = soundSpeed*soundSpeed
	double pressure = 0.0;
*/		
		if(particleType[i] == dummyWall)
		{
			particleBC[i] = other;
			continue;
		}
		// First check based on particle number density
		if(pndSmall[i] < betaPnd) {
		//if(pndi[i] < betaPnd) {
			particleBC[i] = surface;
		}
		else {
			particleBC[i] = inner;
		}

		if(freeSurfType == calcBCType::PND_NEIGH)
		{
			if(pndSmall[i] < betaPnd && numNeigh[i] < betaNeigh) {
			//if(pndi[i] < betaPnd && numNeigh[i] < betaNeigh) {
				particleBC[i] = surface;
			}
			else {
				particleBC[i] = inner;
			}
		}
		else if(freeSurfType == calcBCType::PND_NPCD)
		{
			// Boundary particle verification based on relative distance and weight (NPCD)
			// 2016 - Fluid interface detection technique based on neighborhood particles 
			// centroid deviation (NPCD) for particle methods
			if(particleBC[i] == surface) {
				if(numNeigh[i] > 4 && npcdDeviation2[i] < delta2) {
					particleBC[i] = inner;
				}
			}
		}
		else if(freeSurfType == calcBCType::PND_ARC)
		{
			double normalXi = normal[i*3  ];	double normalYi = normal[i*3+1];	double normalZi = normal[i*3+2];
			double norm2 = normalXi*normalXi + normalYi*normalYi + normalZi*normalZi;
			// 2017 - A multiphase MPS solver for modeling multi-fluid interaction with 
			// free surface and its application in oil spill
			//if((pndSmall[i] >= betaPnd || numNeigh[i] >= betaNeigh) && norm2 <= normThreshold2) {
			/*if(pndSmall[i] >= betaPnd && numNeigh[i] >= betaNeigh && norm2 <= normThreshold2) {
			//if(pndSmall[i] >= betaPnd || numNeigh[i] >= betaNeigh) {
			*/
			if(pndi[i] >= betaPnd && numNeigh[i] >= betaNeigh && norm2 <= normThreshold2) {
				particleBC[i] = inner;
			}
			else {
				double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
				double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
				//double normalXi = normal[i*3  ];	double normalYi = normal[i*3+1];	double normalZi = normal[i*3+2];
				double norm = sqrt(norm2);
				normalXi /= norm;	normalYi /= norm;	normalZi /= norm;
				// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
				double Rref_i[9], normaliw[3], normalMod2;
				// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
				normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
				normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];
				if(normalMod2 > 1.0e-8) {
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

				// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
				Rref_i[0] = 1.0 - 2.0*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2.0*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2.0*normaliw[0]*normaliw[2];
				Rref_i[3] = 0.0 - 2.0*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2.0*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2.0*normaliw[1]*normaliw[2];
				Rref_i[6] = 0.0 - 2.0*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2.0*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2.0*normaliw[2]*normaliw[2];

				// Normal mirror i
				double normalMirrorXi = Rref_i[0]*normalXi + Rref_i[1]*normalYi + Rref_i[2]*normalZi;
				double normalMirrorYi = Rref_i[3]*normalXi + Rref_i[4]*normalYi + Rref_i[5]*normalZi;
				double normalMirrorZi = Rref_i[6]*normalXi + Rref_i[7]*normalYi + Rref_i[8]*normalZi;

				int ix, iy, iz;
				bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
				int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
				for(int jz=minZ;jz<=maxZ;jz++) {
				for(int jy=iy-1;jy<=iy+1;jy++) {
				for(int jx=ix-1;jx<=ix+1;jx++) {
					int jb = jz*numBucketsXY + jy*numBucketsX + jx;
					int j = firstParticleInBucket[jb];
					if(j == -1) continue;
					while(true) {
						double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
						
						// Particle distance r_ij = Xj - Xi_temporary_position
						sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
						// Mirror particle distance r_imj = Xj - Xim_temporary_position
						sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

						// Real particle i and neighbor j
						// If j is inside the neighborhood of i and 
						// is not at the same side of im (avoid real j in the virtual neihborhood)
						if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
						if(j != i) {
							double v0inj, v1inj, v2inj, dstinj2;
							sqrDistBetweenParticles(j, posXi + partDist*normalXi, posYi + partDist*normalYi, 
								posZi + partDist*normalZi, v0inj, v1inj, v2inj, dstinj2);
							
							double rijn = v0ij*normalXi + v1ij*normalYi + v2ij*normalZi;

							double ang = acos(rijn/sqrt(dstij2));

							/*
							if (dstij2 >= dstThreshold2 && dstinj2 < hThreshold2) {
								particleBC[i] = inner;
								//goto endloop;
							}
							else if (dstij2 < dstThreshold2 && ang < 3.14159265/4.0) {
							//else if (dstij2 < dstThreshold2 && rijn*rijn < dstij2*0.5) {
								particleBC[i] = inner;
								//goto endloop;
							}
							else {
								particleBC[i] = surface;
								goto endloop;
							}
							*/
							if (dstij2 >= dstThreshold2 && dstinj2 < hThreshold2){
								particleBC[i] == inner;
								goto endloop;
							}
							else if (dstij2 < dstThreshold2 && ang < thetaArc) {
								particleBC[i] == inner;
								goto endloop;
							}
							else {
								particleBC[i] == surface;
							}
/*
							//if ((dstij2 < dstThreshold2 && ang < thetaArc) || (numNeigh[i] >= betaNeigh)) {
							if (ang < thetaArc && numNeigh[i] >= 4) {
							//else if (dstij2 < dstThreshold2 && rijn*rijn < dstij2*0.5) {
								particleBC[i] = inner;
								goto endloop;
							}
							else {
								particleBC[i] = surface;
								//goto endloop;
							}*/
						}}

						if(wallType == boundaryWallType::POLYGON){
							// Virtual particle i and real neighbor j
							// If j is inside the neighborhood of i and im (intersection) and 
							// is not at the same side of im (avoid real j in the virtual neihborhood)
							if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
								if(j != i) {
									double v0imnj, v1imnj, v2imnj, dstimnj2;
									sqrDistBetweenParticles(j, posMirrorXi + partDist*normalMirrorXi, posMirrorYi + partDist*normalMirrorYi, 
										posMirrorZi + partDist*normalMirrorZi, v0imnj, v1imnj, v2imnj, dstimnj2);
									
									double rimjn = v0imj*normalMirrorXi + v1imj*normalMirrorYi + v2imj*normalMirrorZi;

									double angm = acos(rimjn/sqrt(dstimj2));

									/*
									if (dstimj2 >= dstThreshold2 && dstimnj2 < hThreshold2) {
										particleBC[i] = inner;
										//goto endloop;
									}
									//else if (dstimj2 < dstThreshold2 && rimjn*rimjn < dstimj2*0.6675*0.6675) {
									else if (dstimj2 < dstThreshold2 && angm < thetaArc) {
									//else if (dstimj2 < dstThreshold2 && rimjn*rimjn < dstimj2*0.5) {
										particleBC[i] = inner;
										//goto endloop;
									}
									else {
										particleBC[i] = surface;
										goto endloop;
									}

									*/

									if (dstimj2 >= dstThreshold2 && dstimnj2 < hThreshold2){
										particleBC[i] == inner;
										goto endloop;
									}
									else if (dstimj2 < dstThreshold2 && angm < thetaArc) {
										particleBC[i] == inner;
										goto endloop;
									}
									else {
										particleBC[i] == surface;
									}
									/*
									//if ((dstimj2 < dstThreshold2 && angm < 3.14159265/4.0) || (numNeigh[i] >= betaNeigh)) {
									if (angm < thetaArc && numNeigh[i] >= 4) {
										particleBC[i] = inner;
										goto endloop;
									}
									else {
										particleBC[i] = surface;
										//goto endloop;
									}*/
								}
							}
						}
						j = nextParticleInSameBucket[j];
						if(j == -1) break;
					}
				}}}

				endloop: ;
			}
		}
	}
}

// Compute pressure EMPS (mpsType = calcPressType::EXPLICIT)
void MpsParticle::calcPressEMPS() {
// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		double mi;
		if(PTYPE[i] == 1) mi = DNS_FL1;
		else mi = DNS_FL2;
		//if(particleType[i]==fluid)
		//	mi = Dns[partType::FLUID];
		//else
		//	mi = Dns[partType::WALL];

		double pressure = 0.0;
		if(particleBC[i] == inner) {
			pressure = (pndi[i] - pndSmallZero) * coeffPressEMPS * mi;
		}

//		if(pndSmall[i] < pndThreshold*pndSmallZero && numNeigh[i] < neighThreshold*numNeighZero)
//			particleBC[i] = surface;
//		else
//		{
//			particleBC[i] = inner;
//			pressure = (pndi[i] - pndSmallZero) * coeffPressEMPS * mi;
//		}
//		if(wallType == boundaryWallType::POLYGON) {
//			if(pndi[i] > pndThreshold*pndSmallZero) {
//				pressure = -mi*gravityZ*(0.3-posZi);
//			}
//		}
//		else if(wallType == boundaryWallType::PARTICLE) {
//			if(pndi[i] > pndThreshold*pndSmallZero || numNeigh[i] > neighThreshold*numNeighZero) {
//				pressure = -mi*gravityZ*(0.3-posZi);
//			}
//		}

		if(pressure < 0.0) {
			pressure = 0.0;
		}
		press[i] = pressure;
	}
}

// Compute pressure WCMPS (mpsType = calcPressType::WEAKLY)
void MpsParticle::calcPressWCMPS() {
// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		double mi;
		if(PTYPE[i] == 1) mi = DNS_FL1;
		else mi = DNS_FL2;
		//if(particleType[i]==fluid)
		//	mi = Dns[partType::FLUID];
		//else
		//	mi = Dns[partType::WALL];

		double pressure = 0.0;
		if(particleBC[i] == inner) {
		// if(particleBC[i] == inner && particleType[i] == fluid)
			pressure = (mi*coeffPressWCMPS/gamma)*(pow(pndi[i]/pndSmallZero,gamma)-1);
		}
	
//		if(pndSmall[i] < pndThreshold*pndSmallZero && numNeigh[i] < neighThreshold*numNeighZero)
//			particleBC[i] = surface;
//		else
//		{
//			particleBC[i] = inner;
//			pressure = (mi*coeffPressWCMPS/gamma)*(pow(pndi[i]/pndSmallZero,gamma)-1);
		//pressure = (ni - n0) * coeffPressEMPS * mi;
		//pressure = (pndi[i] - pndSmallZero) * coeffPressEMPS * mi;
//		}
	
//			if(wallType == boundaryWallType::POLYGON) {
//				if(pndi[i] > pndThreshold*pndSmallZero){
//					pressure = -mi*gravityZ*(0.20 - 0.5*partDist -posZi);
//					pressure = -mi*gravityZ*(0.18 - 0.5*partDist -posZi); // lat
//				}
//			}
//		else if(wallType == boundaryWallType::PARTICLE) 
//				if(pndi[i] > pndThreshold*pndSmallZero || numNeigh[i] > neighThreshold*numNeighZero) {
//					pressure = -mi*gravityZ*(0.20 - 0.5*partDist -posZi);
//					pressure = -mi*gravityZ*(0.18 - 0.5*partDist -posZi); // lat
//				}
//			}
		
		if(pressure < 0.0) {
			pressure = 0.0;
		}
		press[i] = pressure;
	}
}

// Solve linear system solver PPE (mpsType = calcPressType::IMPLICIT_PND)
// Perform conjugate gradient method on symmetry matrix A to solve Ax=b
// matA			symmetric (sparse) matrix
// sourceTerm	vector
void MpsParticle::solvePressurePoissonPnd() {

	using T = Eigen::Triplet<double>;
	double lap_r = reL*invPartDist;
	int n_size = (int)(pow(lap_r * 2, dim)); // maximum number of neighbors
	Eigen::SparseMatrix<double> matA(numParticles, numParticles); // declares a column-major sparse matrix type of double
	sourceTerm.resize(numParticles); // Resizing a dynamic-size matrix
	sourceTerm.setZero(); // Right hand side-vector set to zero
	vector<T> coeffs(numParticles * n_size); // list of non-zeros coefficients

// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		if(particleBC[i] == other || particleBC[i] == surface) {
			coeffs.push_back(T(i, i, 1.0));
			continue;
		}

		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double ni = pndSmall[i];
		double sum = 0.0;

		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					double dst = sqrt(dstij2);
					double wL = weight(dst, reL, weightType);
					// coeffPPE = 2.0*dim/(pndLargeZero*lambdaZero)
					double mat_ij = wL*coeffPPE;
					if (particleType[j] == dummyWall) {
						double pgh = RHO[j]*(v0ij*gravityX + v1ij*gravityY + v2ij*gravityZ)*wL;
						sourceTerm(i) -= coeffPPE*pgh;
					}
					else
					{
						sum -= mat_ij;
						if (particleBC[j] == inner) {
							coeffs.push_back(T(i, j, mat_ij));
						}
					}
				}}
				
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		double density = DNS_FL1;
		if(PTYPE[i] != 1) density = DNS_FL2;

		// Increase diagonal
		sum -= alphaCompressibility*density/(timeStep*timeStep);

		//double Cdiag = 1.0;
		//double beta = 0.9;
		//double DI2 = Cdiag*beta*coeffPPE*relaxPND*pndWallContribution[i]*density;
		//sum -= DI2;

		coeffs.push_back(T(i, i, sum));

		//coeffPPESource = relaxPND/(timeStep*timeStep*pndSmallZero)
		sourceTerm(i) = - coeffPPESource*density*(ni - pndSmallZero);

		// 2019 - Enhancement of stabilization of MPS to arbitrary geometries with a generic wall boundary condition
		//double pndc = 0.0;
		//if(ni > 0)
		//	pndc = (ni - pndWallContribution[i])/ni;
		//sourceTerm(i) = pndc*((1.0-relaxPND)*(pndki[i] - ni) + relaxPND*(pndSmallZero - pndski[i]))*ddt/pndSmallZero;
		//sourceTerm(i) = - relaxPND*ddt*(pndski[i] - pndSmallZero)/pndSmallZero;

		//double riw[3], riwSqrt;
		// Normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		//riw[0] = 0.5*(posXi - posMirrorXi); riw[1] = 0.5*(posYi - posMirrorYi); riw[2] = 0.5*(posZi - posMirrorZi);
		//riwSqrt = sqrt(riw[0]*riw[0] + riw[1]*riw[1] + riw[2]*riw[2]);
		//double ST1 = - coeffPPESource*density*(ni - pndSmallZero);
		//double ST2 = 0.0;
		//if(riwSqrt < 0.5*partDist)
		//	ST2 = - Cdiag*(1.0-beta)*2.0*partDist/lambdaZero*(0.5*partDist - riwSqrt)/(timeStep*timeStep);
		//sourceTerm(i) = ST1 + ST2;
	}

	// Finished setup matrix
	matA.setFromTriplets(coeffs.begin(), coeffs.end());
	// Solve PPE
	if(solverType == solvPressType::CG)
		solveConjugateGradient(matA);
	else if (solverType == solvPressType::BICGSTAB)
		solveBiConjugateGradientStabilized(matA);
	// Set zero to negative pressures
	setZeroOnNegativePressure();
}

// Solve linear system solver PPE (mpsType = calcPressType::IMPLICIT_PND_DIVU)
void MpsParticle::solvePressurePoissonPndDivU() {

	using T = Eigen::Triplet<double>;
	double lap_r = reL*invPartDist;
	int n_size = (int)(pow(lap_r * 2, dim)); // maximum number of neighbors
	Eigen::SparseMatrix<double> matA(numParticles, numParticles); // declares a column-major sparse matrix type of double
	sourceTerm.resize(numParticles); // Resizing a dynamic-size matrix
	sourceTerm.setZero(); // Right hand side-vector set to zero
	vector<T> coeffs(numParticles * n_size); // list of non-zeros coefficients

// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		if(particleBC[i] == other || particleBC[i] == surface) {
			coeffs.push_back(T(i, i, 1.0));
			continue;
		}

		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double ni = pndSmall[i];
		double sum = 0.0;

		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					double dst = sqrt(dstij2);
					double wL = weight(dst, reL, weightType);
					// coeffPPE = 2.0*dim/(pndLargeZero*lambdaZero)
					double mat_ij = wL*coeffPPE;
					if (particleType[j] == dummyWall) {
						double pgh = RHO[j]*(v0ij*gravityX + v1ij*gravityY + v2ij*gravityZ)*wL;
						sourceTerm(i) -= coeffPPE*pgh;
					}
					else
					{
						sum -= mat_ij;
						if (particleBC[j] == inner) {
							coeffs.push_back(T(i, j, mat_ij));
						}
					}
				}}
				
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		double density = DNS_FL1;
		if(PTYPE[i] != 1) density = DNS_FL2;

		// Increase diagonal
		sum -= alphaCompressibility*density/(timeStep*timeStep);

		//double Cdiag = 1.0;
		//double beta = 0.9;
		//double DI2 = Cdiag*beta*coeffPPE*relaxPND*pndWallContribution[i]*density;
		//sum -= DI2;

		coeffs.push_back(T(i, i, sum));

		//coeffPPESource = relaxPND/(timeStep*timeStep*pndSmallZero)
		//sourceTerm(i) = - coeffPPESource*density*(ni - pndSmallZero) + (1.0-relaxPND)*density*velDivergence[i]/timeStep;
		sourceTerm(i) = - coeffPPESource*density*(pndki[i] - pndSmallZero) + (1.0-relaxPND)*density*velDivergence[i]/timeStep;
		//sourceTerm(i) = - 4.0*density*(ni - pndSmallZero)/(partDist*partDist*pndSmallZero) 
		//				+ 2.0*density*velDivergence[i]*invPartDist;

		// Sun et al., 2015. Modified MPS method for the 2D fluid structure interaction problem with free surface
		////double dtPhysical = partDist/20.0;
		//double dtPhysical = timeStep;
		//double a1 = fabs(ni - pndSmallZero)/pndSmallZero;
		//if ((pndSmallZero-ni)*velDivergence[i] > 1e-6)
		//{
		//	a1 += dtPhysical*fabs(velDivergence[i]);
		//}
		////double a2 = fabs((ni - pndSmallZero)/pndSmallZero);
		//sourceTerm(i) = - a1*density/(dtPhysical*dtPhysical)*(ni - pndSmallZero)/pndSmallZero 
		//	+ density*velDivergence[i]/dtPhysical;

		// 2019 - Enhancement of stabilization of MPS to arbitrary geometries with a generic wall boundary condition
		//double pndc = 0.0;
		//if(ni > 0)
		//	pndc = (ni - pndWallContribution[i])/ni;
		//sourceTerm(i) = pndc*((1.0-relaxPND)*(pndki[i] - ni) + relaxPND*(pndSmallZero - pndski[i]))*ddt/pndSmallZero;
		//sourceTerm(i) = - relaxPND*ddt*(pndski[i] - pndSmallZero)/pndSmallZero;

		//double riw[3], riwSqrt;
		// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		//riw[0] = 0.5*(posXi - posMirrorXi); riw[1] = 0.5*(posYi - posMirrorYi); riw[2] = 0.5*(posZi - posMirrorZi);
		//riwSqrt = sqrt(riw[0]*riw[0] + riw[1]*riw[1] + riw[2]*riw[2]);
		//double ST1 = - coeffPPESource*density*(ni - pndSmallZero);
		//double ST2 = 0.0;
		//if(riwSqrt < 0.5*partDist)
		//	ST2 = - Cdiag*(1.0-beta)*2.0*partDist/lambdaZero*(0.5*partDist - riwSqrt)/(timeStep*timeStep);
		//sourceTerm(i) = ST1 + ST2;
	}

	// Finished setup matrix
	matA.setFromTriplets(coeffs.begin(), coeffs.end());
	// Solve PPE
	if(solverType == solvPressType::CG)
		solveConjugateGradient(matA);
	else if (solverType == solvPressType::BICGSTAB)
		solveBiConjugateGradientStabilized(matA);
	// Set zero to negative pressures
	setZeroOnNegativePressure();
}

// Solve linear system using Conjugate Gradient (solverType = solvPressType::CG)
void MpsParticle::solveConjugateGradient(Eigen::SparseMatrix<double> p_mat) {
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
	//cg.setTolerance(1.0e-9);
	cg.compute(p_mat);
	if (cg.info() != Eigen::ComputationInfo::Success) {
		cerr << "Error: Failed decompostion." << endl;
	}
	//pressurePPE = cg.solve(sourceTerm);
	pressurePPE = cg.solveWithGuess(sourceTerm,pressurePPE);

#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		press[i] = pressurePPE(i);
	}

	if (cg.info() != Eigen::ComputationInfo::Success) {
		cerr << "Error: Failed solving." << endl;
	}
	solverIter = cg.iterations();
	solverError = cg.error();
}

// Solve linear system using Bi Conjugate Gradient Stabiçized (solverType = solvPressType::BICGSTAB)
void MpsParticle::solveBiConjugateGradientStabilized(Eigen::SparseMatrix<double> p_mat) {
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> bicg;
	//bicg.setTolerance(1.0e-9);
	bicg.compute(p_mat);
	if (bicg.info() != Eigen::ComputationInfo::Success) {
		cerr << "Error: Failed decompostion." << endl;
	}

	pressurePPE = bicg.solveWithGuess(sourceTerm,pressurePPE);

#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		press[i] = pressurePPE(i);
	}

	if (bicg.info() != Eigen::ComputationInfo::Success) {
		cerr << "Error: Failed solving." << endl;
	}
	solverIter = bicg.iterations();
	solverError = bicg.error();
}

// Set negative pressure to zero
void MpsParticle::setZeroOnNegativePressure(){
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		if (press[i] < 0) press[i] = 0.0;
	}
}

// Divergence of velocity
void MpsParticle::calcVelDivergence() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		velDivergence[i] = 0.0;
//		if(particleType[i] == fluid) {
		double DivV = 0.0;
		double ni = pndi[i];
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
					if(j != i) {
						if(particleType[i] == fluid && particleType[j] == fluid) {
							double nj = pndi[j];
							double vijx = vel[j*3  ]-velXi;
							double vijy = vel[j*3+1]-velYi;
							double vijz = vel[j*3+2]-velZi;
							double dst = sqrt(dstij2);
							double wS = weight(dst, reS, weightType);
							if(ni > 1.0e-8)
								DivV += (dim/pndSmallZero)*(nj/ni)*(vijx*v0ij+vijy*v1ij+vijz*v2ij)*wS/dstij2;
						}
					}
				}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		
		velDivergence[i] = DivV;
	}
}

// Divergence of velocity (Polygon wall) - Free-slip
void MpsParticle::calcWallSlipVelDivergence() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
//	if(particleType[i] == fluid) {
		double ni = pndi[i];
		if(particleType[i] == fluid && particleNearWall[i] == true && ni > 1.0e-8) {
			double DivV = 0.0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
			
			// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
			double Rref_i[9], normaliw[3], normalMod2;
			// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
			normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
			normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];
			if(normalMod2 > 1.0e-8) {
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
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
						if(j != i) {
							double dst = sqrt(dstimj2);
							double wS = weight(dst, reS, weightType);
							double nj = pndi[j];
							double vijx = vel[j*3  ]-velMirrorXi;
							double vijy = vel[j*3+1]-velMirrorYi;
							double vijz = vel[j*3+2]-velMirrorZi;
							DivV += (dim/pndSmallZero)*(nj/ni)*(vijx*v0imj+vijy*v1imj+vijz*v2imj)*wS/dstimj2;
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}

			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
			
			if(dstimi2 < reS2) {
				double dst = sqrt(dstimi2);
				double wS = weight(dst, reS, weightType);
				double nj = pndi[i];
				double vijx = velXi-velMirrorXi;
				double vijy = velYi-velMirrorYi;
				double vijz = velZi-velMirrorZi;
				DivV += (dim/pndSmallZero)*(nj/ni)*(vijx*v0imi+vijy*v1imi+vijz*v2imi)*wS/dstimi2;
			}

			velDivergence[i] += DivV;
		}
	}
}


// Divergence of velocity (Polygon wall) - No-slip
void MpsParticle::calcWallNoSlipVelDivergence() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
		double ni = pndi[i];
		//if(particleType[i] == fluid && ni > 1.0e-8) {
		if(particleType[i] == fluid && particleNearWall[i] == true && ni > 1.0e-8) {
	//	if(particleType[i] == fluid) {
			double DivV = 0.0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];

			// Inverse matrix Rinv_i = - I
			double Rinv_i[9], Rref_i[9], normaliw[3], normalMod2;
			// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
			normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
			normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];
			if(normalMod2 > 1.0e-8) {
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

			if(nearMeshType[i] == meshType::FORCED) {
				viwall[0] = velVWall[0];
				viwall[1] = velVWall[1];
				viwall[2] = velVWall[2];
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
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and im (intersection) and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
						if(j != i) {
							double dst = sqrt(dstimj2);
							double wS = weight(dst, reS, weightType);
							double nj = pndi[j];
							double vijx = -(vel[j*3  ]-velMirrorXi);
							double vijy = -(vel[j*3+1]-velMirrorYi);
							double vijz = -(vel[j*3+2]-velMirrorZi);
							// Refelected rij' = Rref_i * ri'j
							double v0m = (Rref_i[0]*v0imj + Rref_i[1]*v1imj + Rref_i[2]*v2imj);
							double v1m = (Rref_i[3]*v0imj + Rref_i[4]*v1imj + Rref_i[5]*v2imj);
							double v2m = (Rref_i[6]*v0imj + Rref_i[7]*v1imj + Rref_i[8]*v2imj);
							DivV += (dim/pndSmallZero)*(nj/ni)*(vijx*v0m+vijy*v1m+vijz*v2m)*wS/dstimj2;
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
					}
			}}}

			// Add "i" contribution ("i" is a neighbor of "mirror i")
			double v0imi, v1imi, v2imi, dstimi2;
			sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
			
			if(dstimi2 < reS2) {
				double dst = sqrt(dstimi2);
				double wS = weight(dst, reS, weightType);
				double nj = pndi[i];
				double vijx = -(velXi-velMirrorXi);
				double vijy = -(velYi-velMirrorYi);
				double vijz = -(velZi-velMirrorZi);
				// Refelected rij' = Rref_i * ri'j
				double v0m = (Rref_i[0]*v0imi + Rref_i[1]*v1imi + Rref_i[2]*v2imi);
				double v1m = (Rref_i[3]*v0imi + Rref_i[4]*v1imi + Rref_i[5]*v2imi);
				double v2m = (Rref_i[6]*v0imi + Rref_i[7]*v1imi + Rref_i[8]*v2imi);
				DivV += (dim/pndSmallZero)*(nj/ni)*(vijx*v0m+vijy*v1m+vijz*v2m)*wS/dstimi2;
			}

			velDivergence[i] += DivV;
		}
	}
}


// Extrapolate pressure to wall and dummy particles
void MpsParticle::extrapolatePressParticlesWallDummy() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == dummyWall || (mpsType == calcPressType::WEAKLY && particleType[i] == wall)) {
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double ni = 0.0;
			double pressure = 0.0;
			
			int ix, iy, iz;
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, dstij2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);

					if(dstij2 < reS2) {
	//				if(j != i) {
						if(j != i && particleType[j] == fluid) {
							double dst = sqrt(dstij2);
							double wS = weight(dst, reS, weightType);
							ni += wS;
							pressure += (press[j] - RHO[j]*(gravityX*v0ij+gravityY*v1ij+gravityZ*v2ij))*wS;
							//pressure += (press[j] + RHO[j]*(gravityX*v0ij+gravityY*v1ij+gravityZ*v2ij))*wS;
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			if(pressure < 0.0)
				pressure = 0.0;
			if(ni > 0)
				press[i] = pressure/ni;
			else
				press[i] = pressure;
		}
	}
}

// Extrapolate pressure to inner particles near polygon walls
void MpsParticle::extrapolatePressParticlesNearPolygonWall() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
	//if(particleType[i] == fluid && press[i] == 0 /*&& particleNearWall[i] == true*/) {
		if(particleType[i] == fluid && press[i] == 0 && particleNearWall[i] == true) {
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
			double pressure = 0.0;
			double sumWij = 0.0;
			int nTotal = 0;
			int nFree = 0;
			
			int ix, iy, iz;
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					// Mirror particle distance r_imj = Xj - Xim_temporary_position
					sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

					// If j is inside the neighborhood of i and 
					// is not at the same side of im (avoid real j in the virtual neihborhood)
					if(dstij2 < reS2 && dstij2 < dstimj2) {
	//				if(dstij2 < 1.2*partDist) {
						if(j != i) {
							double dst = sqrt(dstij2);
							double wS = weight(dst, reS, weightType);
		//					double wS = WEI_WEND(dst, 1.2*partDist);
							sumWij += wS;
		//					pressure += press[j]*wS;
							pressure += (press[j] - RHO[j]*(gravityX*v0ij+gravityY*v1ij+gravityZ*v2ij))*wS;
							nTotal += 1;
							if(particleBC[j] == surface)
								nFree += 1;
						}
					}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			
			if(nTotal > 0)
				numNeighborsSurfaceParticles[i] = double(nFree)/nTotal;
			else
				numNeighborsSurfaceParticles[i] = 1.0;

	//		if(numNeighborsSurfaceParticles[i]<=0.5){
			if(pressure > 0) {
		  		if(sumWij > 1.0e-8)
					press[i] = pressure/sumWij;
				else
					press[i] = pressure;
			}
	//		}
		}
	}
}

// Determinant of matrix
double MpsParticle::detMatrix(double M11, double M12, double M13, double M21, double M22, double M23, double M31, double M32, double M33) {
	return (M11*M22*M33 + M12*M23*M31 + M13*M21*M32)
			- (M13*M22*M31 + M12*M21*M33 + M11*M23*M32);
}

// Inverse of matrix
int MpsParticle::inverseMatrix(int dim, double &M11, double &M12, double &M13, double &M21, double &M22, double &M23, double &M31, double &M32, double &M33) {
	double M[3][3], Maux[3][3];

	Maux[0][0] = M11;	Maux[0][1] = M12;	Maux[0][2] = M13;
	Maux[1][0] = M21;	Maux[1][1] = M22;	Maux[1][2] = M23;
	Maux[2][0] = M31;	Maux[2][1] = M32;	Maux[2][2] = M33;

	if(dim == 2)
		Maux[2][2] = 1.0;

	// Convert matrix to identity
	for(int i = 0; i < 3; i++)
	for(int j = 0; j < 3; j++) {
		if(i == j) M[i][j] = 1.0;
		else M[i][j] = 0.0;
	}

	for(int k = 0; k < dim; k++) {

		if(fabs(Maux[k][k]) <= 1.0e-8) {
			M11 = 1.0;	M12 = 0.0;	M13 = 0.0;
			M21 = 0.0;	M22 = 1.0;	M23 = 0.0;
			M31 = 0.0;	M32 = 0.0;	M33 = 1.0;
			return 0;
		}

		for(int i = 0; i < dim; i++) {
			if(i == k)
			continue;

			double m = Maux[i][k]/Maux[k][k];

			for(int j = 0; j < dim; j++) {
				Maux[i][j] -= m*Maux[k][j];
				M[i][j] -= m*M[k][j];
			}
		}
	}

	for(int i = 0; i < dim; i++)
		for(int j = 0; j < dim; j++)
			M[i][j] /= Maux[i][i];

	M11 = M[0][0];	M12 = M[0][1];	M13 = M[0][2];
	M21 = M[1][0];	M22 = M[1][1];	M23 = M[1][2];
	M31 = M[2][0];	M32 = M[2][1];	M33 = M[2][2];

	return 1;
}

// Correction matrix
void MpsParticle::correctionMatrix() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
//	if(particleType[i] == fluid) {
		correcMatrixRow1[i*3  ] = 0.0; correcMatrixRow1[i*3+1] = 0.0; correcMatrixRow1[i*3+2] = 0.0;
		correcMatrixRow2[i*3  ] = 0.0; correcMatrixRow2[i*3+1] = 0.0; correcMatrixRow2[i*3+2] = 0.0;
		correcMatrixRow3[i*3  ] = 0.0; correcMatrixRow3[i*3+1] = 0.0; correcMatrixRow3[i*3+2] = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					double dst = sqrt(dstij2);
					double wS = weightGradient(dst, reS, weightType);
					double invDstij2 = 1.0/dstij2;
					correcMatrixRow1[i*3  ] += wS*v0ij*v0ij*invDstij2;	correcMatrixRow1[i*3+1] += wS*v0ij*v1ij*invDstij2;	correcMatrixRow1[i*3+2] += wS*v0ij*v2ij*invDstij2;
					correcMatrixRow2[i*3  ] += wS*v1ij*v0ij*invDstij2;	correcMatrixRow2[i*3+1] += wS*v1ij*v1ij*invDstij2;	correcMatrixRow2[i*3+2] += wS*v1ij*v2ij*invDstij2;
					correcMatrixRow3[i*3  ] += wS*v2ij*v0ij*invDstij2;	correcMatrixRow3[i*3+1] += wS*v2ij*v1ij*invDstij2;	correcMatrixRow3[i*3+2] += wS*v2ij*v2ij*invDstij2;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// coeffPressGrad is a negative cte (-dim/noGrad)
		correcMatrixRow1[i*3  ] *= -coeffPressGrad;	correcMatrixRow1[i*3+1] *= -coeffPressGrad;	correcMatrixRow1[i*3+2] *= -coeffPressGrad;
		correcMatrixRow2[i*3  ] *= -coeffPressGrad;	correcMatrixRow2[i*3+1] *= -coeffPressGrad;	correcMatrixRow2[i*3+2] *= -coeffPressGrad;
		correcMatrixRow3[i*3  ] *= -coeffPressGrad;	correcMatrixRow3[i*3+1] *= -coeffPressGrad;	correcMatrixRow3[i*3+2] *= -coeffPressGrad;

		// Inverse of the matrix
		int rcv = inverseMatrix((int)(dim),correcMatrixRow1[i*3],correcMatrixRow1[i*3+1],correcMatrixRow1[i*3+2],correcMatrixRow2[i*3],correcMatrixRow2[i*3+1],correcMatrixRow2[i*3+2],correcMatrixRow3[i*3],correcMatrixRow3[i*3+1],correcMatrixRow3[i*3+2]);

//		if(i == 200) {
//			printf("\n X %e %e %e ", correcMatrixRow1[i*3  ], correcMatrixRow1[i*3+1], correcMatrixRow1[i*3+2]);
//			printf("\n Y %e %e %e ", correcMatrixRow2[i*3  ], correcMatrixRow2[i*3+1], correcMatrixRow2[i*3+2]);
//			printf("\n Z %e %e %e \n", correcMatrixRow3[i*3  ], correcMatrixRow3[i*3+1], correcMatrixRow3[i*3+2]);
//		}
	}
}

// Acceleration due to pressure gradient
void MpsParticle::calcPressGradient() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
//	if(particleType[i] == fluid) {
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double pressMin = press[i];
		double Pi = press[i];
		double ni = pndi[i];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		if(gradientType == 0 || gradientType == 2) {
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					if(pressMin > press[j]) pressMin = press[j];
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}}
		
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					double dst = sqrt(dstij2);
					double wS = weightGradient(dst, reS, weightType);
					if(gradientType == 0)
						wS *= (press[j] - pressMin)/dstij2;
					else if(gradientType == 1)
						wS *= (press[j] + Pi)/dstij2;
					else if(gradientType == 2)
						wS *= (press[j] + Pi - 2.0*pressMin)/dstij2;
					else if(gradientType == 3) {
						double nj = pndi[j];
						if(ni > 1.0e-8 && nj > 1.0e-8)
							wS *= (ni*press[j]/nj + nj*Pi/ni)/dstij2;
					}
					accX += v0ij*wS;	accY += v1ij*wS;	accZ += v2ij*wS;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// coeffPressGrad is a negative cte (-dim/noGrad)
		// Original
//		acc[i*3  ]=relaxPress*accX*invDns[partType::FLUID]*coeffPressGrad;
//		acc[i*3+1]=relaxPress*accY*invDns[partType::FLUID]*coeffPressGrad;
//		acc[i*3+2]=relaxPress*accZ*invDns[partType::FLUID]*coeffPressGrad;
		// Modified
		if(gradientCorrection == false) {
			acc[i*3  ]=relaxPress*accX*coeffPressGrad/RHO[i];
			acc[i*3+1]=relaxPress*accY*coeffPressGrad/RHO[i];
			acc[i*3+2]=relaxPress*accZ*coeffPressGrad/RHO[i];
		}
		else {
		//	if(correcMatrixRow1[1*3] > 1.0) {
		//		printf("\n X %e %e %e ", correcMatrixRow1[i*3  ], correcMatrixRow1[i*3+1], correcMatrixRow1[i*3+2]);
		//		printf("\n Y %e %e %e ", correcMatrixRow2[i*3  ], correcMatrixRow2[i*3+1], correcMatrixRow2[i*3+2]);
		//		printf("\n Z %e %e %e \n", correcMatrixRow3[i*3  ], correcMatrixRow3[i*3+1], correcMatrixRow3[i*3+2]);
			//}
			acc[i*3  ]=(relaxPress*coeffPressGrad/RHO[i])*(accX*correcMatrixRow1[i*3] + accY*correcMatrixRow1[i*3+1] + accZ*correcMatrixRow1[i*3+2]);
			acc[i*3+1]=(relaxPress*coeffPressGrad/RHO[i])*(accX*correcMatrixRow2[i*3] + accY*correcMatrixRow2[i*3+1] + accZ*correcMatrixRow2[i*3+2]);
			acc[i*3+2]=(relaxPress*coeffPressGrad/RHO[i])*(accX*correcMatrixRow3[i*3] + accY*correcMatrixRow3[i*3+1] + accZ*correcMatrixRow3[i*3+2]);
		}
	}
}

// Acceleration due to pressure gradient (Polygon wall)
void MpsParticle::calcWallPressGradient() {
	//int nPartNearMesh = partNearMesh.size();
	double VolumeForce = pow(partDist,dim);
	// Maximum velocity is the minimum of the computed and expected maximum velocities
	double maxVelocity = min(velMax, expectMaxVelocity);
	double velMax2 = maxVelocity*maxVelocity;
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {

	//particleNearWall[i]=true; // Only to show particles near polygon
	
	//if(particleType[i] == fluid) {
	if(particleType[i] == fluid && particleNearWall[i] == true) {
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double pressMin = press[i];
		double Pi = press[i];
		double ni = pndi[i];
		// Wall gradient Mitsume`s model
		double Rref_i[9], normaliw[3], normaliwSqrt;
		// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
		normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);

		if(normaliwSqrt > 1.0e-8) {
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
		Rai[0] = Rref_i[0]*accStar[i*3] + Rref_i[1]*accStar[i*3+1] + Rref_i[2]*accStar[i*3+2];
		Rai[1] = Rref_i[3]*accStar[i*3] + Rref_i[4]*accStar[i*3+1] + Rref_i[5]*accStar[i*3+2];
		Rai[2] = Rref_i[6]*accStar[i*3] + Rref_i[7]*accStar[i*3+1] + Rref_i[8]*accStar[i*3+2];

		// if(i == 16107)
		// {
		// 	printf("\ni:%5d timeCurrent: %lf / Rref_i: ", i, timeCurrent);
		// 	for(int rr=0; rr<9; rr++)
		// 		printf("%lf ", Rref_i[rr]);
		// 	printf("\ni:%5d timeCurrent: %lf / Rai: ", i, timeCurrent);
		// 	for(int rr=0; rr<3; rr++)
		// 		printf("%lf ", Rai[rr]);
		// 	printf("\ni:%5d timeCurrent: %lf / ai: %lf, %lf, %lf", i, timeCurrent, accStar[i*3], accStar[i*3+1], accStar[i*3+2]);
			
		// }
			
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		if(gradientType == 0 || gradientType == 2) {
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i) {
					if(pressMin > press[j]) pressMin = press[j];
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}}
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i) {

					double dst = sqrt(dstimj2);
					double wS = weightGradient(dst, reS, weightType);

					// Taylor pressure Pj
					double Pj = Pi + RHO[i]*(Rai[0]*v0imj + Rai[1]*v1imj + Rai[2]*v2imj);
					Pj = press[j];

					// if(i == 16107)
					// 	printf("\ni:%5d j:%5d timeCurrent: %lf / Pj: %lf / Pi: %lf / Zj: %lf / Zi: %lf", i, j, timeCurrent, Pj, Pi, pos[j*3+2], posMirrorZi);

					if(gradientType == 0)
						wS *= (Pj - pressMin)/dstimj2;//(press[j] - pressMin)/dstimj2
					else if(gradientType == 1)
						wS *= (Pj + Pi)/dstimj2;//(press[j] + press[i])/dstimj2
					else if(gradientType == 2)
						wS *= (Pj + Pi - 2.0*pressMin)/dstimj2;//(press[j] + press[i] - 2.0*pressMin)/dstimj2
					else if(gradientType == 3) {
						double nj = pndi[j];
						if(ni > 1.0e-8 && nj > 1.0e-8)
							wS *= (ni*Pj/nj + nj*Pi/ni)/dstimj2;
					}
					accX += v0imj*wS;	accY += v1imj*wS;	accZ += v2imj*wS;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// Add "i" contribution ("i" is a neighbor of "mirror i")
		double v0imi, v1imi, v2imi, dstimi2;
		sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
	  	
		if(dstimi2 < reS2) {

			double dst = sqrt(dstimi2);
			double wS = weightGradient(dst, reS, weightType);

			// Taylor pressure Pj
			double Pj = Pi + RHO[i]*(Rai[0]*v0imi + Rai[1]*v1imi + Rai[2]*v2imi);
			Pj = press[i];

			//if(i == 16107)
			//	printf("\ni:%5d timeCurrent: %lf / Pj: %lf / Pi: %lf / Zj: %lf / Zi: %lf", i, timeCurrent, Pj, Pi, posZi, posMirrorZi);

			if(gradientType == 0)
				wS *= (Pj - pressMin)/dstimi2;//(press[i] - pressMin)/dstimi2
			else if(gradientType == 1)
				wS *= (Pj + Pi)/dstimi2;//(press[i] + press[i])/dstimi2;
			else if(gradientType == 2)
				wS *= (Pj + Pi - 2.0*pressMin)/dstimi2;//(press[i] + press[i] - 2.0*pressMin)/dstimi2;
			else if(gradientType == 3) {
				double nj = pndi[i];
				if(ni > 1.0e-8 && nj > 1.0e-8)
					wS *= (ni*Pj/nj + nj*Pi/ni)/dstimi2;
			}
			accX += v0imi*wS;	accY += v1imi*wS;	accZ += v2imi*wS;
	  	}

		// Repulsive force
	  	double rpsForce[3];
	  	rpsForce[0]=rpsForce[1]=rpsForce[2] = 0.0;

	  	if(repulsiveForceType == repForceType::HARADA) {
			// Parallel analysis system for free-surface flow using MPS method with explicitly represented polygon wall boundary model
			// https://doi.org/10.1007/s40571-019-00269-6
			if(normaliwSqrt < reRepulsiveForce && normaliwSqrt > 1.0e-8) {
				double wijRep = RHO[i]/(timeStep*timeStep)*(reRepulsiveForce-normaliwSqrt);
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
			}
		}
		else if(repulsiveForceType == repForceType::MITSUME) {
			// Explicitly represented polygon wall boundary model for the explicit MPS method
			// https://doi.org/10.1007/s40571-015-0037-8
			if(normaliwSqrt < reRepulsiveForce && normaliwSqrt > 1.0e-8) {
				double wijRep = repForceCoefMitsume*weightGradient(normaliwSqrt, reRepulsiveForce, weightType);
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
			}
		}
		else if(repulsiveForceType == repForceType::LENNARD_JONES) {
			// Simulating Free Surface Flows with SPH
			// https://doi.org/10.1006/jcph.1994.1034
			if(normaliwSqrt < reRepulsiveForce && normaliwSqrt > 1.0e-8) {
				double R1 = (reRepulsiveForce/normaliwSqrt)*(reRepulsiveForce/normaliwSqrt);
				double R2 = R1*R1;
				double wijRep = (repForceCoefLennardJones*velMax2/normaliwSqrt)*(R2-R1)*RHO[i];
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
			}
		}
		else {
			// SPH particle boundary forces for arbitrary boundaries 
			// https://doi.org/10.1016/j.cpc.2009.05.008
			if(normaliwSqrt < reRepulsiveForce && normaliwSqrt > 1.0e-8) {
				double W1 = (1.0+3.0/2.0*normaliwSqrt/(reRepulsiveForce));
				double W2 = (1.0-normaliwSqrt/(reRepulsiveForce))*(1.0-normaliwSqrt/(reRepulsiveForce))*(1.0-normaliwSqrt/(reRepulsiveForce));
				double wijRep = (repForceCoefMonaghanKajtar*velMax2/(normaliwSqrt - 0.0*partDist))*(1.0/8.0)*(W1)*(W2)*RHO[i];
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
			}
	  	}
		// coeffPressGrad is a negative cte (-dim/noGrad)
		// Original
//		acc[i*3  ] += (relaxPress*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*coeffPressGrad - rpsForce[0])*invDns[partType::FLUID];
//		acc[i*3+1] += (relaxPress*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*coeffPressGrad - rpsForce[1])*invDns[partType::FLUID];
//		acc[i*3+2] += (relaxPress*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*coeffPressGrad - rpsForce[2])*invDns[partType::FLUID];
		// Modified
		acc[i*3  ] += (relaxPress*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*coeffPressGrad - rpsForce[0])/RHO[i];
		acc[i*3+1] += (relaxPress*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*coeffPressGrad - rpsForce[1])/RHO[i];
		acc[i*3+2] += (relaxPress*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*coeffPressGrad - rpsForce[2])/RHO[i];

		// FSI
		// Force on wall
		forceWall[i*3  ] += - (relaxPress*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*coeffPressGrad - rpsForce[0])*VolumeForce;
		forceWall[i*3+1] += - (relaxPress*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*coeffPressGrad - rpsForce[1])*VolumeForce;
		forceWall[i*3+2] += - (relaxPress*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*coeffPressGrad - rpsForce[2])*VolumeForce;
	}}
}

// Calculation of the volume of fraction if phase II in the mixture
void MpsParticle::calcVolumeFraction()
{
	if(Fraction_method == 1) {   //Linear distribution
#pragma omp parallel for schedule(dynamic,64)
		for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
			double sum1 = 0, sum2 = 0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			
			int ix, iy, iz;
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, dstij2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
					
					if(dstij2 < reS2) {
					if(j != i && particleType[j] == fluid) {
						sum1 = sum1 + 1;
						if(PTYPE[j] >= 2) sum2 = sum2 + 1;
						}}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			if(sum1 == 0)
				Cv[i] = 0.0;
			else 
				Cv[i] = sum2 / sum1;
		}}
	}
	else if(Fraction_method == 2) {   //Non linear :  Smoothed using the weight funtion
#pragma omp parallel for schedule(dynamic,64)
		for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
			double sum1 = 0, sum2 = 0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			
			int ix, iy, iz;
			bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
			int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
			for(int jz=minZ;jz<=maxZ;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*numBucketsXY + jy*numBucketsX + jx;
				int j = firstParticleInBucket[jb];
				if(j == -1) continue;
				while(true) {
					double v0ij, v1ij, v2ij, dstij2;
					
					// Particle distance r_ij = Xj - Xi_temporary_position
					sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);

					if(dstij2 < reS2) {
					if(j != i && particleType[j] == fluid) {
						double dst = sqrt(dstij2);
						double wS = weight(dst, reS, weightType);
						sum1 = sum1 + wS;
						if(PTYPE[j] >= 2) sum2 = sum2 + wS;
						}}
					j = nextParticleInSameBucket[j];
					if(j == -1) break;
				}
			}}}
			if(sum1 == 0)
				Cv[i] = 0.0;
			else 
				Cv[i] = sum2 / sum1;
		}}
	}
}

//void NonNwtVscTrm_omp(double *x_vel, double *y_vel, double *z_vel) {
// Viscosity interaction values for "real" fluid particles
void MpsParticle::calcViscosityInteractionVal() {

	//double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;//,grain_VF, *p_smooth;
	double d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
//	Xmin = 0.0 - partDist*3; Xmax = 1.65 + partDist*3;
//	Ymin = 0.0 - partDist*3; Ymax = 0.15 + partDist*3;
	Zmin = 0.0 - partDist*3; //Zmax = 0.7 + partDist*30;
	// damErosion3D
	Xmin = 0.0 - partDist*3; Xmax = 2.00 + partDist*3;
	Ymin = 0.0 - partDist*3; Ymax = 0.10 + partDist*3;
	// Changed !!!

	// Search free-surface particles for each interval of aa = 2 particles in wall
	int aa = 2, kx, ky;
	int kx_max = int((Xmax - Xmin) / aa * invPartDist) + 1;
	int ky_max = int((Ymax - Ymin) / aa * invPartDist) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
/*
	S11 = new double[numParticles + 1];
	S22 = new double[numParticles + 1];
	S33 = new double[numParticles + 1];
	S12 = new double[numParticles + 1];
	S13 = new double[numParticles + 1];
	S23 = new double[numParticles + 1];
*/
	//p_smooth = new double[numParticles + 1];
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
	for(kx = 1; kx <= kx_max; kx++)
	{
		for(ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		
			kx = int((posXi - Xmin) / aa * invPartDist) + 1;
			ky = int((posYi - Ymin) / aa * invPartDist) + 1;


			//if(posZi > BL[kx][ky] && Cv[i] > 0.5) { BL[kx][ky] = posZi; PS[kx][ky] = pnew[i]; }
			if(posZi > BL[kx][ky] && Cv[i] > 0.5) { BL[kx][ky] = posZi; PS[kx][ky] = press[i]; }
			if(posZi > WL[kx][ky] && PTYPE[i] == 1) { WL[kx][ky] = posZi; }
		}
	}

	// Strain rate calculation
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
	if(particleType[i] == fluid) {
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;

//		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					double dst = sqrt(dstij2);
					double wS = weight(dst, reS, weightType);
					double vec_ijx = vel[j*3  ]-velXi;	
					double vec_ijy = vel[j*3+1]-velYi;	
					double vec_ijz = vel[j*3+2]-velZi;
					double invDstij2 = 1.0/dstij2;

					sum1 += vec_ijx*v0ij*wS*invDstij2;
					sum2 += vec_ijx*v1ij*wS*invDstij2;
					sum3 += vec_ijx*v2ij*wS*invDstij2;
					
					sum4 += vec_ijy*v0ij*wS*invDstij2;
					sum5 += vec_ijy*v1ij*wS*invDstij2;
					sum6 += vec_ijy*v2ij*wS*invDstij2;
					
					sum7 += vec_ijz*v0ij*wS*invDstij2;
					sum8 += vec_ijz*v1ij*wS*invDstij2;
					sum9 += vec_ijz*v2ij*wS*invDstij2;
					
					sum10 += press[j]*wS;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// coeffPressGrad is a negative cte (-dim/pndGradientZero)
		Uxx = -coeffPressGrad*sum1; Uxy = -coeffPressGrad*sum2; Uxz = -coeffPressGrad*sum3;
		Uyx = -coeffPressGrad*sum4; Uyy = -coeffPressGrad*sum5; Uyz = -coeffPressGrad*sum6;
		Uzx = -coeffPressGrad*sum7; Uzy = -coeffPressGrad*sum8; Uzz = -coeffPressGrad*sum9;

		p_smooth[i] = sum10 / pndGradientZero;
		if(p_smooth[i] < 0) p_smooth[i] = 0;

		S11[i] = 0.5*(Uxx + Uxx);
		S12[i] = 0.5*(Uxy + Uyx);
		S13[i] = 0.5*(Uxz + Uzx);
		S22[i] = 0.5*(Uyy + Uyy);
		S23[i] = 0.5*(Uyz + Uzy);
		S33[i] = 0.5*(Uzz + Uzz);

		//II[i] = 0.5*Uxx*Uxx + 0.5*Uyy*Uyy + 0.25*(Uxy + Uyx)*(Uxy + Uyx);
		//II[i] = 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		II[i] = 0.5*(S11[i] * S11[i] + 2 * S12[i] * S12[i] + 2 * S13[i] * S13[i] + S22[i] * S22[i] + 2 * S23[i] * S23[i] + S33[i] * S33[i]);
		//II[i]= S11[i]*S22[i] +S22[i]*S33[i]+ S11[i]*S33[i] - S12[i]*S12[i] -S13[i]*S13[i]- S23[i]*S23[i];
		if(II[i] < 0 || II[i] * 0 != 0) II[i] = 0;
		//II=fabs(S11[i]*S22[i]-S12[i]*S12[i]);
		
		//std::cout << " II: " << II[i] << std::endl;
	}}

	// Newtonian viscosity
	if(fluidType == viscType::NEWTONIAN)
	{
#pragma omp parallel for
		for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
			if(PTYPE[i] <= 1)MEU[i] = KNM_VS1 * DNS_FL1;
			if(PTYPE[i] != 1)MEU[i] = KNM_VS2 * DNS_FL2;
		}}

//		if(TURB > 0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2.0*sqrt(II[i]);

//			if(NEUt[i] * 0 != 0)  NEUt[i] = 0;
//			if(NEUt[i] > 1)     NEUt[i] = 1;
//		}
	}

	// Granular Fluid
	if(fluidType == viscType::NON_NEWTONIAN)
	{
	// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
//#pragma omp parallel for
		for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {

//			double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];

			if(PTYPE[i] == 1) {
				MEU[i] = KNM_VS1 * DNS_FL1;
			}
			else if(PTYPE[i] == 2) {
				// phi: internal friction angle
				// phi2: maximum friction angle
				phi = (Cv[i] - 0.25)*PHI / (1.0 - 0.25);
				phi2 = (Cv[i] - 0.25)*PHI_2 / (1.0 - 0.25);
				if(Cv[i] <= 0.25) { phi = 1.0e-8; phi2 = 1.0e-8; } // phi close to zero
				if(PTYPE[i] <= 0) phi = PHI_BED; // ghost

				// normal stress calculation (mechanical pressure)
				p_rheo_new[i] = p_smooth[i];

				kx = int((posXi - Xmin) / aa * invPartDist) + 1;
				ky = int((posYi - Ymin) / aa * invPartDist) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress=(BL[k]-posYi+DL/2.0)*(DNS_FL2)*9.81;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - posZi + partDist / 2.0)*(DNS_FL2 - DNS_FL1)*9.81 - (velXi*velXi + velYi*velYi + velZi*velZi)*(DNS_FL2 - DNS_FL1) / 2.0;	// normal_stress= Gama.H

				if(p_smooth[i] - (WL[kx][ky] - posZi)*DNS_FL1*9.81<0) p_smooth[i] = (WL[kx][ky] - posZi)*DNS_FL1*9.81;
				if(timeCurrent <= 1) normal_stress = 1.0*(1.0 - timeCurrent)*(p_smooth[i] - (WL[kx][ky] - posZi)*DNS_FL1*9.81) + 1.0*(timeCurrent)*normal_stress;

//				normal_stress = p_smooth[i] - (WL[kx][ky] - posZi)*DNS_FL1*9.81;
//				normal_stress = p_smooth[i];
				//normal_stress=normal_stress*0.61*1500/DNS_FL2;
				if(normal_stress < 1 || Cv[i] < 0.5) normal_stress = 1;

				p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/DNS_SDT);		// Free-fall (dry granular material)
				Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/(DNS_FL1*Cd));	// Grain inertia (submerged)
				//Inertia[i] = sqrt(II[i])*(KNM_VS1*DNS_FL1)/normal_stress ;	// Viscous regime

//				Inertia[i] = 1.0;
				// VF_max VF_min
				VF[i] = VF_max - (VF_max - VF_min)*Inertia[i];
				if(VF[i] < VF_min) VF[i] = VF_min;
				RHO[i] = DNS_SDT * VF[i] + (1.0-VF[i])*DNS_FL1;
				phi = phi * VF[i] / VF_max;

				// Mohr-Coulomb
				double yield_stress = cohes * cos(phi) + normal_stress * sin(phi);

				if(yield_stress < 0) yield_stress = 0;

				double visc_max = (yield_stress*mm*0.5 + MEU0);

				if(II[i] > 0)
					MEU_Y[i] = yield_stress * (1.0 - exp(-mm * sqrt(II[i]))) / 2.0 / sqrt(II[i]);
				else
					MEU_Y[i] = visc_max;

				// H-B rheology

				//meu_0 = MEU0;

				// Non-linear Meu(I) rheology
				//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/DNS_FL2)+sqrt(II[i])*DG);			//free fall
				meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/(DNS_FL1*Cd))+sqrt(II[i])*DG);		//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(KNM_VS1*DNS_FL1)/(I0*normal_stress+sqrt(II[i])*(KNM_VS1*DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL2)/I0;		//free fall
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL1*Cd)/I0;	//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*(KNM_VS1*DNS_FL1)/I0;					//viscous

				if(II[i] <= 0 || (meu_0 * 0) != 0) meu_0 = MEU0;

				visc_max = (yield_stress*mm*0.5 + meu_0);

				//if(isnan(II[i]) || isinf(II[i])) {
					//std::cout << " viscmax: " << II[i] << std::endl;
				//	assert(visc_max >= 0 || visc_max <= 0);
				//}
				
				// Herschel bulkley papanastasiou
				MEU[i] = MEU_Y[i] + MEU0 * pow(4 * II[i], (N - 1) / 2);

				// MEU_Y rheological model
				//MEU[i] = MEU_Y[i] + meu_0;
				
				if(II[i] == 0 || MEU[i]>visc_max) {
					//std::cout << " MEU>viscmax: " << yield_stress*mm*0.5 << " meu0: " << meu_0 << " II: " << II[i] << std::endl;
					MEU[i] = visc_max;
				}
				if(PTYPE[i] <= 0) MEU[i] = MEU[i] * Cv[i] + DNS_FL1*KNM_VS1*(1.0 - Cv[i]);

				//if(MEU[i]/RHO[i] > maxVIS) maxVIS = MEU[i]/RHO[i];
			}
			
			if(PTYPE[i] >= 2) {
				if(Cv[i] > 0.5) RHO[i] = DNS_FL2;
				else RHO[i] = Cv[i] * DNS_FL2 + (1.0 - Cv[i]) * DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if(stress_cal_method == 2)
//		{
//			for(i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
//				for(l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if(i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2 * MEU[i] * MEU[j] / (MEU[i] + MEU[j]);
//						if((NEUt[i] + NEUt[j])>0) meuij = meuij + 2 * NEUt[i] * RHO[i] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j]);

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

//				Tau_xx[i] = (dim / n0) * 2 * sum1;
//				Tau_yy[i] = (dim / n0) * 2 * sum5;
//				Tau_zz[i] = (dim / n0) * 2 * sum9;

//				Tau_xy[i] = (dim / n0)*(sum2 + sum4);
//				Tau_xz[i] = (dim / n0)*(sum3 + sum7);
//				Tau_yz[i] = (dim / n0)*(sum6 + sum8);
//			}
//		}

	} // if(fluidType == viscType::NON_NEWTONIAN)

	//---------------------------------------------------------------

//	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS; //delete[]p_smooth;
//	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL; //p_smooth = NULL;

	delete[]BL; delete[]WL; delete[]PS;
	BL = NULL; WL = NULL; PS = NULL;
}

// Free-slip condition. Viscosity interaction values
void MpsParticle::calcWallSlipViscosityInteractionVal() {

	//double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;//, grain_VF, *p_smooth;
	double d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
//	Xmin = 0.0 - partDist*3; Xmax = 1.65 + partDist*3;
//	Ymin = 0.0 - partDist*3; Ymax = 0.15 + partDist*3;
	Zmin = 0.0 - partDist*3; //Zmax = 0.7 + partDist*30;
	// damErosion3D
	Xmin = 0.0 - partDist*3; Xmax = 2.00 + partDist*3;
	Ymin = 0.0 - partDist*3; Ymax = 0.10 + partDist*3;
	// Changed !!!

	// Search free-surface particles for each interval of aa = 2 particles in wall
	int aa = 2, kx, ky;
	int kx_max = int((Xmax - Xmin) / aa * invPartDist) + 1;
	int ky_max = int((Ymax - Ymin) / aa * invPartDist) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
	double aUxx, aUxy, aUxz, aUyx, aUyy, aUyz, aUzx, aUzy, aUzz;
/*
	S11 = new double[numParticles + 1];
	S22 = new double[numParticles + 1];
	S33 = new double[numParticles + 1];
	S12 = new double[numParticles + 1];
	S13 = new double[numParticles + 1];
	S23 = new double[numParticles + 1];
*/
	//p_smooth = new double[numParticles + 1];
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
	for(kx = 1; kx <= kx_max; kx++)
	{
		for(ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			
			kx = int((posXi - Xmin) / aa * invPartDist) + 1;
			ky = int((posYi - Ymin) / aa * invPartDist) + 1;

			//if(posZi>BL[kx][ky] && Cv[i]>0.5) { BL[kx][ky] = posZi; PS[kx][ky] = pnew[i]; }
			if(posZi>BL[kx][ky] && Cv[i]>0.5) { BL[kx][ky] = posZi; PS[kx][ky] = press[i]; }
			if(posZi>WL[kx][ky] && PTYPE[i] == 1) { WL[kx][ky] = posZi; }
		}
	}

	// Strain rate calculation
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
	//if(particleType[i] == fluid) {
	if(particleType[i] == fluid && particleNearWall[i] == true) {
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;

//		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2

		// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
		double Rref_i[9], normaliw[3], normalMod2;
		// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
		normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

		if(normalMod2 > 1.0e-8) {
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
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i) {
					double dst = sqrt(dstimj2);
					double wS = weight(dst, reS, weightType);

					double vec_mijx = vel[j*3  ]-velMirrorXi;	
					double vec_mijy = vel[j*3+1]-velMirrorYi;	
					double vec_mijz = vel[j*3+2]-velMirrorZi;

					sum1 += vec_mijx*v0imj*wS/dstimj2;
					sum2 += vec_mijx*v1imj*wS/dstimj2;
					sum3 += vec_mijx*v2imj*wS/dstimj2;
					
					sum4 += vec_mijy*v0imj*wS/dstimj2;
					sum5 += vec_mijy*v1imj*wS/dstimj2;
					sum6 += vec_mijy*v2imj*wS/dstimj2;
					
					sum7 += vec_mijz*v0imj*wS/dstimj2;
					sum8 += vec_mijz*v1imj*wS/dstimj2;
					sum9 += vec_mijz*v2imj*wS/dstimj2;
					
					sum10 += press[j]*wS;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Rref_i * gradU
		// coeffPressGrad is a negative cte (-dim/pndGradientZero)
		aUxx = -coeffPressGrad*(Rref_i[0]*sum1 + Rref_i[1]*sum4 + Rref_i[2]*sum7);
		aUxy = -coeffPressGrad*(Rref_i[0]*sum2 + Rref_i[1]*sum5 + Rref_i[2]*sum8);
		aUxz = -coeffPressGrad*(Rref_i[0]*sum3 + Rref_i[1]*sum6 + Rref_i[2]*sum9);

		aUyx = -coeffPressGrad*(Rref_i[3]*sum1 + Rref_i[4]*sum4 + Rref_i[5]*sum7);
		aUyy = -coeffPressGrad*(Rref_i[3]*sum2 + Rref_i[4]*sum5 + Rref_i[5]*sum8);
		aUyz = -coeffPressGrad*(Rref_i[3]*sum3 + Rref_i[4]*sum6 + Rref_i[5]*sum9);

		aUzx = -coeffPressGrad*(Rref_i[6]*sum1 + Rref_i[7]*sum4 + Rref_i[8]*sum7);
		aUzy = -coeffPressGrad*(Rref_i[6]*sum2 + Rref_i[7]*sum5 + Rref_i[8]*sum8);
		aUzz = -coeffPressGrad*(Rref_i[6]*sum3 + Rref_i[7]*sum6 + Rref_i[8]*sum9);

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
		p_smooth[i] += sum10 / pndGradientZero;
		if(p_smooth[i]<0) p_smooth[i] = 0;

		// 0.5*(gradU + gradUt)
		S11[i] = 0.5*(Uxx + Uxx);
		S12[i] = 0.5*(Uxy + Uyx);
		S13[i] = 0.5*(Uxz + Uzx);
		S22[i] = 0.5*(Uyy + Uyy);
		S23[i] = 0.5*(Uyz + Uzy);
		S33[i] = 0.5*(Uzz + Uzz);

		//II[i] = 0.5*Uxx*Uxx + 0.5*Uyy*Uyy + 0.25*(Uxy + Uyx)*(Uxy + Uyx);

		// Addition of II for particles near mesh
		//II[i] += 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		II[i] += 0.5*(S11[i] * S11[i] + 2 * S12[i] * S12[i] + 2 * S13[i] * S13[i] + S22[i] * S22[i] + 2 * S23[i] * S23[i] + S33[i] * S33[i]);
//		II[i] = 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		//II[i]= S11[i]*S22[i] +S22[i]*S33[i]+ S11[i]*S33[i] - S12[i]*S12[i] -S13[i]*S13[i]- S23[i]*S23[i] ;
		if(II[i]<0 || II[i] * 0 != 0) II[i] = 0;
		//II=fabs(S11[i]*S22[i]-S12[i]*S12[i]);
	}}
	
	// Newtonian viscosity
	if(fluidType == viscType::NEWTONIAN)
	{
	// Loop only for particles near mesh
#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<numParticles; i++) {
		//if(particleType[i] == fluid) {
		if(particleType[i] == fluid && particleNearWall[i] == true) {
			if(PTYPE[i] <= 1)MEU[i] = KNM_VS1 * DNS_FL1;
			if(PTYPE[i] != 1)MEU[i] = KNM_VS2 * DNS_FL2;
		}}

//		if(TURB>0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2.0*sqrt(II[i]);

//			if(NEUt[i] * 0 != 0)  NEUt[i] = 0;
//			if(NEUt[i]>1)     NEUt[i] = 1;
//		}
	}

	// Granular Fluid
	if(fluidType == viscType::NON_NEWTONIAN)
	{
		// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
		// Loop only for particles near mesh
//#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<numParticles; i++) {
		//if(particleType[i] == fluid) {
		if(particleType[i] == fluid && particleNearWall[i] == true) {
//			double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];

			if(PTYPE[i] == 1) {
				MEU[i] = KNM_VS1 * DNS_FL1;
			}
			else if(PTYPE[i] == 2) {
				phi = (Cv[i] - 0.25)*PHI / (1.0 - 0.25);
				phi2 = (Cv[i] - 0.25)*PHI_2 / (1.0 - 0.25);
				if(Cv[i] <= 0.25) { phi = 1.0e-8; phi2 = 1.0e-8; } // phi close to zero
				if(PTYPE[i] <= 0) phi = PHI_BED;

				// normal stress calculation (mehcanical pressure)
				p_rheo_new[i] = p_smooth[i];

				kx = int((posXi - Xmin) / aa * invPartDist) + 1;
				ky = int((posYi - Ymin) / aa * invPartDist) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress=(BL[k]-posYi+DL/2.0)*(DNS_FL2)*9.81;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - posZi + partDist / 2.0)*(DNS_FL2 - DNS_FL1)*9.81 - (velXi*velXi + velYi*velYi + velZi*velZi)*(DNS_FL2 - DNS_FL1) / 2.0;	// normal_stress= Gama.H

				if(p_smooth[i] - (WL[kx][ky] - posZi)*DNS_FL1*9.81<0) p_smooth[i] = (WL[kx][ky] - posZi)*DNS_FL1*9.81;
				if(timeCurrent <= 1) normal_stress = 1.0*(1.0 - timeCurrent)*(p_smooth[i] - (WL[kx][ky] - posZi)*DNS_FL1*9.81) + 1.0*(timeCurrent)*normal_stress;

//				normal_stress = p_smooth[i] - (WL[kx][ky] - posZi)*DNS_FL1*9.81;
//				normal_stress = p_smooth[i];
				//normal_stress=normal_stress*0.61*1500/DNS_FL2;
				if(normal_stress < 1 || Cv[i] < 0.5) normal_stress = 1;

				p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/DNS_SDT);		// Free-fall (dry granular material)
				Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/(DNS_FL1*Cd));	// Grain inertia (submerged)
				//Inertia[i] = sqrt(II[i])*(KNM_VS1*DNS_FL1)/normal_stress ;	// Viscous regime

				// VF_max VF_min
				VF[i] = VF_max - (VF_max - VF_min)*Inertia[i];
				if(VF[i] < VF_min) VF[i] = VF_min;
				RHO[i] = DNS_SDT * VF[i] + (1.0-VF[i])*DNS_FL1;
				phi = phi * VF[i] / VF_max;

				double yield_stress = cohes * cos(phi) + normal_stress * sin(phi);

				if(yield_stress < 0) yield_stress = 0;

				double visc_max = (yield_stress*mm*0.5 + MEU0);

				if(II[i]>0)
					MEU_Y[i] = yield_stress * (1.0 - exp(-mm * sqrt(II[i]))) / 2.0 / sqrt(II[i]);
				else
					MEU_Y[i] = visc_max;

				// H-B rheology

				//meu_0 = MEU0;

				// Non-linear Meu(I) rheology
				//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/DNS_FL2)+sqrt(II[i])*DG);			//free fall
				meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/(DNS_FL1*Cd))+sqrt(II[i])*DG);		//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(KNM_VS1*DNS_FL1)/(I0*normal_stress+sqrt(II[i])*(KNM_VS1*DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL2)/I0;		//free fall
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL1*Cd)/I0;	//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*(KNM_VS1*DNS_FL1)/I0;					//viscous

				if(II[i] <= 0 || (meu_0 * 0) != 0) meu_0 = MEU0;

				visc_max = (yield_stress*mm*0.5 + meu_0);

				// Herschel bulkley papanastasiou
				MEU[i] = MEU_Y[i] + MEU0 * pow(4 * II[i], (N - 1) / 2);

				// MEU_Y rheological model
				//MEU[i] = MEU_Y[i] + meu_0;
				
				if(II[i] == 0 || MEU[i]>visc_max) MEU[i] = visc_max;
				if(PTYPE[i] <= 0) MEU[i] = MEU[i] * Cv[i] + DNS_FL1*KNM_VS1*(1.0 - Cv[i]);
			}
			
			if(PTYPE[i] >= 2) {
				if(Cv[i] > 0.5) RHO[i] = DNS_FL2;
				else RHO[i] = Cv[i] * DNS_FL2 + (1.0 - Cv[i]) * DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if(stress_cal_method == 2)
//		{
//			for(i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
//				for(l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if(i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2 * MEU[i] * MEU[j] / (MEU[i] + MEU[j]);
//						if((NEUt[i] + NEUt[j])>0) meuij = meuij + 2 * NEUt[i] * RHO[i] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j]);

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

//				Tau_xx[i] = (dim / n0) * 2 * sum1;
//				Tau_yy[i] = (dim / n0) * 2 * sum5;
//				Tau_zz[i] = (dim / n0) * 2 * sum9;

//				Tau_xy[i] = (dim / n0)*(sum2 + sum4);
//				Tau_xz[i] = (dim / n0)*(sum3 + sum7);
//				Tau_yz[i] = (dim / n0)*(sum6 + sum8);
//			}
//		}

	} // if(fluidType == viscType::NON_NEWTONIAN)

	//---------------------------------------------------------------

//	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS;// delete[]p_smooth;
//	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL;// p_smooth = NULL;

	delete[]BL; delete[]WL; delete[]PS;
	BL = NULL; WL = NULL; PS = NULL;
}

// No-Slip condition. Viscosity interaction values
void MpsParticle::calcWallNoSlipViscosityInteractionVal() {

	//double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;//, grain_VF, *p_smooth;
	double d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
//	Xmin = 0.0 - partDist*3; Xmax = 1.65 + partDist*3;
//	Ymin = 0.0 - partDist*3; Ymax = 0.15 + partDist*3;
	Zmin = 0.0 - partDist*3; //Zmax = 0.7 + partDist*30;
	// damErosion3D
	Xmin = 0.0 - partDist*3; Xmax = 2.00 + partDist*3;
	Ymin = 0.0 - partDist*3; Ymax = 0.10 + partDist*3;
	// Changed !!!

	// Search free-surface particles for each interval of aa = 2 particles in wall
	int aa = 2, kx, ky;
	int kx_max = int((Xmax - Xmin) / aa * invPartDist) + 1;
	int ky_max = int((Ymax - Ymin) / aa * invPartDist) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
/*
	S11 = new double[numParticles + 1];
	S22 = new double[numParticles + 1];
	S33 = new double[numParticles + 1];
	S12 = new double[numParticles + 1];
	S13 = new double[numParticles + 1];
	S23 = new double[numParticles + 1];
*/
	//p_smooth = new double[numParticles + 1];
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
	for(kx = 1; kx <= kx_max; kx++)
	{
		for(ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		
			kx = int((posXi - Xmin) / aa * invPartDist) + 1;
			ky = int((posYi - Ymin) / aa * invPartDist) + 1;

			//if(posZi>BL[kx][ky] && Cv[i]>0.5) { BL[kx][ky] = posZi; PS[kx][ky] = pnew[i]; }
			if(posZi>BL[kx][ky] && Cv[i]>0.5) { BL[kx][ky] = posZi; PS[kx][ky] = press[i]; }
			if(posZi>WL[kx][ky] && PTYPE[i] == 1) { WL[kx][ky] = posZi; }
		}
	}

	// Strain rate calculation
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
	//if(particleType[i] == fluid) {
	if(particleType[i] == fluid && particleNearWall[i] == true) {
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;

//		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2];

		// Transformation matrix R_i = I
		double Rref_i[9], Rinv_i[9], normaliw[3], normalMod2;
	    // normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
	    normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
	    normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

	    if(normalMod2 > 1.0e-8) {
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

		if(nearMeshType[i] == meshType::FORCED) {
			viwall[0] = velVWall[0];
			viwall[1] = velVWall[1];
			viwall[2] = velVWall[2];
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
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i) {
					double dst = sqrt(dstimj2);
					double wS = weight(dst, reS, weightType);

					double vec_mijx = vel[j*3  ]-velMirrorXi;	
					double vec_mijy = vel[j*3+1]-velMirrorYi;	
					double vec_mijz = vel[j*3+2]-velMirrorZi;

					sum1 += vec_mijx*v0imj*wS/dstimj2;
					sum2 += vec_mijx*v1imj*wS/dstimj2;
					sum3 += vec_mijx*v2imj*wS/dstimj2;
					
					sum4 += vec_mijy*v0imj*wS/dstimj2;
					sum5 += vec_mijy*v1imj*wS/dstimj2;
					sum6 += vec_mijy*v2imj*wS/dstimj2;
					
					sum7 += vec_mijz*v0imj*wS/dstimj2;
					sum8 += vec_mijz*v1imj*wS/dstimj2;
					sum9 += vec_mijz*v2imj*wS/dstimj2;
					
					sum10 += press[j]*wS;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Rinv_i * gradU * Rref_i = - gradU * Rref_i
		// coeffPressGrad is a negative cte (-dim/pndGradientZero)
		Uxx = coeffPressGrad*(sum1*Rref_i[0] + sum2*Rref_i[3] + sum3*Rref_i[6]);
		Uxy = coeffPressGrad*(sum1*Rref_i[1] + sum2*Rref_i[4] + sum3*Rref_i[7]);
		Uxz = coeffPressGrad*(sum1*Rref_i[2] + sum2*Rref_i[5] + sum3*Rref_i[8]);

		Uyx = coeffPressGrad*(sum4*Rref_i[0] + sum5*Rref_i[3] + sum6*Rref_i[6]);
		Uyy = coeffPressGrad*(sum4*Rref_i[1] + sum5*Rref_i[4] + sum6*Rref_i[7]);
		Uyz = coeffPressGrad*(sum4*Rref_i[2] + sum5*Rref_i[5] + sum6*Rref_i[8]);

		Uzx = coeffPressGrad*(sum7*Rref_i[0] + sum8*Rref_i[3] + sum9*Rref_i[6]);
		Uzy = coeffPressGrad*(sum7*Rref_i[1] + sum8*Rref_i[4] + sum9*Rref_i[7]);
		Uzz = coeffPressGrad*(sum7*Rref_i[2] + sum8*Rref_i[5] + sum9*Rref_i[8]);

		// Addition of smoothed pressure for particles near mesh
		p_smooth[i] += sum10 / pndGradientZero;
		if(p_smooth[i]<0) p_smooth[i] = 0;

		// - (Rref_i * gradU) - (Rref_i * gradU)t
		S11[i] = 0.5*(Uxx + Uxx);
		S12[i] = 0.5*(Uxy + Uyx);
		S13[i] = 0.5*(Uxz + Uzx);
		S22[i] = 0.5*(Uyy + Uyy);
		S23[i] = 0.5*(Uyz + Uzy);
		S33[i] = 0.5*(Uzz + Uzz);

		//II[i] = 0.5*Uxx*Uxx + 0.5*Uyy*Uyy + 0.25*(Uxy + Uyx)*(Uxy + Uyx);
		
		// Addition of II for particles near mesh
		//II[i] += 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		II[i] += 0.5*(S11[i] * S11[i] + 2 * S12[i] * S12[i] + 2 * S13[i] * S13[i] + S22[i] * S22[i] + 2 * S23[i] * S23[i] + S33[i] * S33[i]);
//		II[i] = 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		//II[i]= S11[i]*S22[i] +S22[i]*S33[i]+ S11[i]*S33[i] - S12[i]*S12[i] -S13[i]*S13[i]- S23[i]*S23[i] ;
		if(II[i]<0 || II[i] * 0 != 0) II[i] = 0;
		//II=fabs(S11[i]*S22[i]-S12[i]*S12[i]);
	}}
	
	// Newtonian viscosity
	if(fluidType == viscType::NEWTONIAN)
	{
	// Loop only for particles near mesh
#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<numParticles; i++) {
		//if(particleType[i] == fluid) {
		if(particleType[i] == fluid && particleNearWall[i] == true) {
			if(PTYPE[i] <= 1)MEU[i] = KNM_VS1 * DNS_FL1;
			if(PTYPE[i] != 1)MEU[i] = KNM_VS2 * DNS_FL2;
		}}

//		if(TURB>0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2.0*sqrt(II[i]);

//			if(NEUt[i] * 0 != 0)  NEUt[i] = 0;
//			if(NEUt[i]>1)     NEUt[i] = 1;
//		}
	}

	// Granular Fluid
	if(fluidType == viscType::NON_NEWTONIAN)
	{
		// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
		// Loop only for particles near mesh
//#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<numParticles; i++) {
		//if(particleType[i] == fluid) {
		if(particleType[i] == fluid && particleNearWall[i] == true) {
//			double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
			double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
			double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];

			if(PTYPE[i] == 1) 
				MEU[i] = KNM_VS1 * DNS_FL1;
			else if(PTYPE[i] == 2) {
				phi = (Cv[i] - 0.25)*PHI / (1.0 - 0.25);
				phi2 = (Cv[i] - 0.25)*PHI_2 / (1.0 - 0.25);
				if(Cv[i] <= 0.25) { phi = 1.0e-8; phi2 = 1.0e-8; } // phi close to zero
				if(PTYPE[i] <= 0) phi = PHI_BED;

				// normal stress calculation (mechanical pressure)
				p_rheo_new[i] = p_smooth[i];

				kx = int((posXi - Xmin) / aa * invPartDist) + 1;
				ky = int((posYi - Ymin) / aa * invPartDist) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress=(BL[k]-posYi+DL/2.0)*(DNS_FL2)*9.81;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - posZi + partDist / 2.0)*(DNS_FL2 - DNS_FL1)*9.81 - (velXi*velXi + velYi*velYi + velZi*velZi)*(DNS_FL2 - DNS_FL1) / 2.0;	// normal_stress= Gama.H

				if(p_smooth[i] - (WL[kx][ky] - posZi)*DNS_FL1*9.81<0) p_smooth[i] = (WL[kx][ky] - posZi)*DNS_FL1*9.81;
				if(timeCurrent <= 1) normal_stress = 1.0*(1.0 - timeCurrent)*(p_smooth[i] - (WL[kx][ky] - posZi)*DNS_FL1*9.81) + 1.0*(timeCurrent)*normal_stress;

//				normal_stress = p_smooth[i] - (WL[kx][ky] - posZi)*DNS_FL1*9.81;
//				normal_stress = p_smooth[i];
				//normal_stress=normal_stress*0.61*1500/DNS_FL2;
				if(normal_stress < 1 || Cv[i] < 0.5) normal_stress = 1;

				p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/DNS_SDT);		// Free-fall (dry granular material)
				Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/(DNS_FL1*Cd));	// Grain inertia (submerged)
				//Inertia[i] = sqrt(II[i])*(KNM_VS1*DNS_FL1)/normal_stress ;	// Viscous regime

				// VF_max VF_min
				VF[i] = VF_max - (VF_max - VF_min)*Inertia[i];
				if(VF[i] < VF_min) VF[i] = VF_min;
				RHO[i] = DNS_SDT * VF[i] + (1.0-VF[i])*DNS_FL1;
				phi = phi * VF[i] / VF_max;

				double yield_stress = cohes * cos(phi) + normal_stress * sin(phi);

				if(yield_stress < 0) yield_stress = 0;

				double visc_max = (yield_stress*mm*0.5 + MEU0);

				if(II[i]>0)
					MEU_Y[i] = yield_stress * (1.0 - exp(-mm * sqrt(II[i]))) / 2.0 / sqrt(II[i]);
				else
					MEU_Y[i] = visc_max;

				// H-B rheology

				//meu_0 = MEU0;

				// Non-linear Meu(I) rheology
				//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/DNS_FL2)+sqrt(II[i])*DG);			//free fall
				meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/(DNS_FL1*Cd))+sqrt(II[i])*DG);		//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(KNM_VS1*DNS_FL1)/(I0*normal_stress+sqrt(II[i])*(KNM_VS1*DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL2)/I0;		//free fall
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL1*Cd)/I0;	//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*(KNM_VS1*DNS_FL1)/I0;					//viscous

				if(II[i] <= 0 || (meu_0 * 0) != 0) meu_0 = MEU0;

				visc_max = (yield_stress*mm*0.5 + meu_0);

				// Herschel bulkley papanastasiou
				MEU[i] = MEU_Y[i] + MEU0 * pow(4 * II[i], (N - 1) / 2);

				// MEU_Y rheological model
				//MEU[i] = MEU_Y[i] + meu_0;
				
				if(II[i] == 0 || MEU[i]>visc_max) MEU[i] = visc_max;
				if(PTYPE[i] <= 0) MEU[i] = MEU[i] * Cv[i] + DNS_FL1*KNM_VS1*(1.0 - Cv[i]);
			}
			
			if(PTYPE[i] >= 2) {
				if(Cv[i] > 0.5) RHO[i] = DNS_FL2;
				else RHO[i] = Cv[i] * DNS_FL2 + (1.0 - Cv[i]) * DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if(stress_cal_method == 2)
//		{
//			for(i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
//				for(l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if(i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2 * MEU[i] * MEU[j] / (MEU[i] + MEU[j]);
//						if((NEUt[i] + NEUt[j])>0) meuij = meuij + 2 * NEUt[i] * RHO[i] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j]);

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

//				Tau_xx[i] = (dim / n0) * 2 * sum1;
//				Tau_yy[i] = (dim / n0) * 2 * sum5;
//				Tau_zz[i] = (dim / n0) * 2 * sum9;

//				Tau_xy[i] = (dim / n0)*(sum2 + sum4);
//				Tau_xz[i] = (dim / n0)*(sum3 + sum7);
//				Tau_yz[i] = (dim / n0)*(sum6 + sum8);
//			}
//		}

	} // if(fluidType == viscType::NON_NEWTONIAN)

	//---------------------------------------------------------------

//	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS; //delete[]p_smooth;
//	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL;// p_smooth = NULL;
	
	delete[]BL; delete[]WL; delete[]PS;
	BL = NULL; WL = NULL; PS = NULL;
}

// Free-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
void MpsParticle::calcWallSlipViscosity() {
	//int nPartNearMesh = partNearMesh.size();
	double VolumeForce = pow(partDist,dim);
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
	//if(particleType[i] == fluid) {
	if(particleType[i] == fluid && particleNearWall[i] == true) {
		
		double meu_i = MEU[i];
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2

		// Transformation matrix Rref_i = I - 2.0*normal_iwall*normal_iwall
		double Rref_i[9], normaliw[3], normalMod2;
		// Normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
		normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

		if(normalMod2 > 1.0e-8) {
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
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && dstimj2 < reL2 && dstij2 < dstimj2) {
				if(j != i) {
					double dst = sqrt(dstimj2);
					double wL = weight(dst, reL, weightType);

					if((meu_i + MEU[j]) > 1.0e-8)
						NEU = 2 * meu_i * MEU[j] / (meu_i + MEU[j]);
					else
						NEU = 0.0;

//NEU = KNM_VS2 * DNS_FL2;

					if(PTYPE[i] == 1) NEU = NEU/DNS_FL1;
					else NEU = NEU/DNS_FL2;

					//if((NEUt[i] + NEUt[j])>0) NEU = NEU + (2 * NEUt[i] * RHO[j] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j])) / RHO[i];

					// Original
//					accX +=(vel[j*3  ]-velMirrorXi)*w;
//					accY +=(vel[j*3+1]-velMirrorYi)*w;
//					accZ +=(vel[j*3+2]-velMirrorZi)*w;
					// Modified
					accX +=(vel[j*3  ]-velMirrorXi)*wL*NEU;
					accY +=(vel[j*3+1]-velMirrorYi)*wL*NEU;
					accZ +=(vel[j*3+2]-velMirrorZi)*wL*NEU;

					//accX +=(Velk[j*3  ]-velMirrorXi)*w;
					//accY +=(Velk[j*3+1]-velMirrorYi)*w;
					//accZ +=(Velk[j*3+2]-velMirrorZi)*w;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Add "i" contribution ("i" is a neighbor of "mirror i")
		double v0imi, v1imi, v2imi, dstimi2;
		sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
		
		if(dstimi2 < reL2) {
			double dst = sqrt(dstimi2);
			double wL = weight(dst, reL, weightType);

			if(meu_i > 1.0e-8)
				NEU = 2 * meu_i * meu_i / (meu_i + meu_i);
			else
				NEU = 0.0;

//NEU = KNM_VS2 * DNS_FL2;

			if(PTYPE[i] == 1) NEU = NEU/DNS_FL1;
			else NEU = NEU/DNS_FL2;

			// Original
//			accX +=(velXi-velMirrorXi)*w;
//			accY +=(velYi-velMirrorYi)*w;
//			accZ +=(velZi-velMirrorZi)*w;

			// Modified
			accX +=(velXi-velMirrorXi)*wL*NEU;
			accY +=(velYi-velMirrorYi)*wL*NEU;
			accZ +=(velZi-velMirrorZi)*wL*NEU;
		}

		// Wall laplacian Mitsume`s model
		// Correction of velocity
		// Original
//      acc[i*3  ] += (Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*coeffViscosity;
//		acc[i*3+1] += (Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*coeffViscosity;
//		acc[i*3+2] += (Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*coeffViscosity;

		// coeffViscMultiphase = 2.0*dim/(pndLargeZero*lambdaZero);
		// Modified
		acc[i*3  ] += (Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*coeffViscMultiphase;
		acc[i*3+1] += (Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*coeffViscMultiphase;
		acc[i*3+2] += (Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*coeffViscMultiphase;

		
		// FSI
		// Force on wall
		forceWall[i*3  ] += - (Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*coeffViscMultiphase*VolumeForce*RHO[i];
		forceWall[i*3+1] += - (Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*coeffViscMultiphase*VolumeForce*RHO[i];
		forceWall[i*3+2] += - (Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*coeffViscMultiphase*VolumeForce*RHO[i];
	}}
}

// No-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
void MpsParticle::calcWallNoSlipViscosity() {
	//int nPartNearMesh = partNearMesh.size();
	double VolumeForce = pow(partDist,dim);
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
	//if(particleType[i] == fluid) {
	if(particleType[i] == fluid && particleNearWall[i] == true) {

		double meu_i = MEU[i];
		double accX = 0.0;			double accY = 0.0;			double accZ = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
//		double velXi = Velk[i*3  ];	double velYi = Velk[i*3+1];	double velZi = Velk[i*3+2];

		// Inverse matrix Rinv_i = - I
		double Rinv_i[9], normaliw[3], normalMod2;
		// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
		normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

		if(normalMod2 > 1.0e-8) {
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

		if(nearMeshType[i] == meshType::FORCED) {
			viwall[0] = velVWall[0];
			viwall[1] = velVWall[1];
			viwall[2] = velVWall[2];
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
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && dstimj2 < reL2 && dstij2 < dstimj2) {
				if(j != i) {
					double dst = sqrt(dstimj2);
					double wL = weight(dst, reL, weightType);

					if((meu_i + MEU[j]) > 1.0e-8)
						NEU = 2 * meu_i * MEU[j] / (meu_i + MEU[j]);
					else
						NEU = 0.0;

//NEU = KNM_VS2 * DNS_FL2;


					if(PTYPE[i] == 1) NEU = NEU/DNS_FL1;
					else NEU = NEU/DNS_FL2;
					
					//if((NEUt[i] + NEUt[j])>0) NEU = NEU + (2 * NEUt[i] * RHO[j] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j])) / RHO[i];

					// Original
//					accX +=(vel[j*3  ]-velMirrorXi)*w;
//					accY +=(vel[j*3+1]-velMirrorYi)*w;
//					accZ +=(vel[j*3+2]-velMirrorZi)*w;
					// Modified
					accX +=(vel[j*3  ]-velMirrorXi)*wL*NEU;
					accY +=(vel[j*3+1]-velMirrorYi)*wL*NEU;
					accZ +=(vel[j*3+2]-velMirrorZi)*wL*NEU;
					
					//accX +=(Velk[j*3  ]-velMirrorXi)*w;
					//accY +=(Velk[j*3+1]-velMirrorYi)*w;
					//accZ +=(Velk[j*3+2]-velMirrorZi)*w;

					//if(i==2817) {
					//	Fwall[i*3  ] += 1;//AA[0];
					//	Fwall[i*3+1] += 1;//AA[1];
					//	Fwall[i*3+2] += 1;//AA[2];
					//	std::cout << j << " " << Velk[j*3] << " " << Velk[j*3+1] << " " << Velk[j*3+2] << " Pj " << press[j] << std::endl;
					//}
					//accX += (Rinv_i[0]*(Velk[j*3  ]-velMirrorXi)+ Rinv_i[1]*(Velk[j*3+1]-velMirrorYi) + Rinv_i[2]*(Velk[j*3+2]-velMirrorZi))*w;
					//accY += (Rinv_i[3]*(Velk[j*3  ]-velMirrorXi)+ Rinv_i[4]*(Velk[j*3+1]-velMirrorYi) + Rinv_i[5]*(Velk[j*3+2]-velMirrorZi))*w;
					//accZ += (Rinv_i[6]*(Velk[j*3  ]-velMirrorXi)+ Rinv_i[7]*(Velk[j*3+1]-velMirrorYi) + Rinv_i[8]*(Velk[j*3+2]-velMirrorZi))*w;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Add "i" contribution ("i" is a neighbor of "mirror i")
		double v0imi, v1imi, v2imi, dstimi2;
		sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
		
		if(dstimi2 < reL2) {
			double dst = sqrt(dstimi2);
			double wL = weight(dst, reL, weightType);

			if(meu_i > 1.0e-8)
				NEU = 2 * meu_i * meu_i / (meu_i + meu_i);
			else
				NEU = 0.0;

//NEU = KNM_VS2 * DNS_FL2;

			if(PTYPE[i] == 1) NEU = NEU/DNS_FL1;
			else NEU = NEU/DNS_FL2;

			// Original
//			accX +=(velXi-velMirrorXi)*w;
//			accY +=(velYi-velMirrorYi)*w;
//			accZ +=(velZi-velMirrorZi)*w;

			// Modified
			accX +=(velXi-velMirrorXi)*wL*NEU;
			accY +=(velYi-velMirrorYi)*wL*NEU;
			accZ +=(velZi-velMirrorZi)*wL*NEU;

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
		// coeffViscMultiphase = 2.0*dim/(pndLargeZero*lambdaZero);
		// Original
//     	acc[i*3  ] += (Rinv_i[0]*accX + Rinv_i[1]*accY + Rinv_i[2]*accZ)*coeffViscosity;
//		acc[i*3+1] += (Rinv_i[3]*accX + Rinv_i[4]*accY + Rinv_i[5]*accZ)*coeffViscosity;
//		acc[i*3+2] += (Rinv_i[6]*accX + Rinv_i[7]*accY + Rinv_i[8]*accZ)*coeffViscosity;
		// Modified
		acc[i*3  ] += (Rinv_i[0]*accX + Rinv_i[1]*accY + Rinv_i[2]*accZ)*coeffViscMultiphase;
		acc[i*3+1] += (Rinv_i[3]*accX + Rinv_i[4]*accY + Rinv_i[5]*accZ)*coeffViscMultiphase;
		acc[i*3+2] += (Rinv_i[6]*accX + Rinv_i[7]*accY + Rinv_i[8]*accZ)*coeffViscMultiphase;
		//Acv[i*3  ] = (Rinv_i[0]*accX + Rinv_i[1]*accY + Rinv_i[2]*accZ)*coeffViscosity;
		//Acv[i*3+1] = (Rinv_i[3]*accX + Rinv_i[4]*accY + Rinv_i[5]*accZ)*coeffViscosity;
		//Acv[i*3+2] = (Rinv_i[6]*accX + Rinv_i[7]*accY + Rinv_i[8]*accZ)*coeffViscosity;

		//double AA[3];
		
		//AA[0] = (Rinv_i[0]*accX + Rinv_i[1]*accY + Rinv_i[2]*accZ)*coeffViscosity;
		//AA[1] = (Rinv_i[3]*accX + Rinv_i[4]*accY + Rinv_i[5]*accZ)*coeffViscosity;
		//AA[2] = (Rinv_i[6]*accX + Rinv_i[7]*accY + Rinv_i[8]*accZ)*coeffViscosity;

		//wallParticleForce1[i*3  ] = AA[0];
		//wallParticleForce1[i*3+1] = AA[1];
		//wallParticleForce1[i*3+2] = AA[2];

		//if(i==2817) {
			//Fwall[i*3  ] = AA[0];
			//Fwall[i*3+1] = AA[1];
			//Fwall[i*3+2] = AA[2];
			//std::cout << "t: " << timeCurrent << std::endl;
			//std::cout << "Fwall " << AA[0] << " " << AA[1] << " " << AA[2] << std::endl;
			//std::cout << "Veli " << velXi << " " << velYi << " " << velZi << std::endl;
			//std::cout << "Posmi " << posMirrorXi << " " << posMirrorYi << " " << posMirrorZi << std::endl;
			//std::cout << "pndi " << pndi[i] << " Pi " << press[i] << std::endl;
			//printf("Accx:%lf Accy:%lf Accz:%lf \n", AA[0],AA[1],AA[2]);
		//}
		//if(i==6) {
			//printf("Accx:%lf Accy:%lf Accz:%lf \n", acc[i*3],acc[i*3+1],acc[i*3+2]);
		//	printf("Time:%e\n", timeCurrent);
		//	printf("Xi:%e %e %e Xm:%e %e %e\n", posXi,posYi,posZi,posMirrorXi,posMirrorYi,posMirrorZi);
		//	printf("Vi:%e %e %e Vm:%e %e %e\n", velXi,velYi,velZi,velMirrorXi,velMirrorYi,velMirrorZi);
		//	printf("acc:%e %e %e\n", AA[0],AA[1],AA[2]);
		//}


		// FSI
		// Force on wall
		forceWall[i*3  ] += - (Rinv_i[0]*accX + Rinv_i[1]*accY + Rinv_i[2]*accZ)*coeffViscMultiphase*VolumeForce*RHO[i];
		forceWall[i*3+1] += - (Rinv_i[3]*accX + Rinv_i[4]*accY + Rinv_i[5]*accZ)*coeffViscMultiphase*VolumeForce*RHO[i];
		forceWall[i*3+2] += - (Rinv_i[6]*accX + Rinv_i[7]*accY + Rinv_i[8]*accZ)*coeffViscMultiphase*VolumeForce*RHO[i];
	}}
}

// 3D triangle to xy plane
// https://math.stackexchange.com/questions/856666/how-can-i-transform-a-3d-triangle-to-xy-plane
void MpsParticle::transformMatrix(double *V1, double *V2, double *V3, double *RM) {
	double A[3],B[3],C[3],U[3],V[3],W[3];
	// Translate the vertex "V1" to origin 0,0,0
	for(unsigned int i = 0; i < 3; i++) {
		//A[i] = V1[i] - V1[i];
		B[i] = V2[i] - V1[i];
		C[i] = V3[i] - V1[i];
	}
	// Define the vector U pointing from A to B, and normalise it.
	double normB = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
	if(normB > 1.0e-8) {
		for(unsigned int i = 0; i < 3; i++)
			U[i] = B[i]/normB;
	}
	else {
		for(unsigned int i = 0; i < 3; i++)
			U[i] = 0.0;
	}
	// Take the cross product which is at right angles to the triangle.
	W[0] = U[1] * C[2] - U[2] * C[1];
	W[1] = U[2] * C[0] - U[0] * C[2];
	W[2] = U[0] * C[1] - U[1] * C[0];
	// Normalise it to give the unit vector
	double normW = sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
	if(normW > 1.0e-8) {
		for(unsigned int i = 0; i < 3; i++)
			W[i] = W[i]/normW;
	}
	else {
		for(unsigned int i = 0; i < 3; i++)
			W[i] = 0.0;
	}
	// Find the cross-product automatically a unit vector.
	V[0] = U[1] * W[2] - U[2] * W[1];
	V[1] = U[2] * W[0] - U[0] * W[2];
	V[2] = U[0] * W[1] - U[1] * W[0];
	// In coordinates corresponding to the basis {U, V, W }, the triangle lies. The rotation matrix that carries the usual
	// (1,0,0) to U, the usual (0,1,0) to V, and the usual (,0,0,1) to W has U, V and W as its columns.
	// You want the other direction, so take the inverse - which for a rotation matrix is just the transpose.
	for(unsigned int i = 0; i < 3; i++) {
		RM[i  ] = U[i];
		RM[i+3] = V[i];
		RM[i+6] = W[i];
	}

	/*
	// Projection of P represented in triangle coordinate system 3D (X,Y,Z) -> 2D (U,V,W).
	proj_pn = RM*proj_p';
	// Square approximation support in triangle coordinate system 3D (X,Y,Z) -> 2D (U,V,W). Here, we assume the
	// value in W as 0.
	ls = sqrt(re^2-p_pb^2);
	squareSupport = [[-ls+proj_pn(1);ls+proj_pn(1);ls+proj_pn(1);-ls+proj_pn(1)],...
	[-ls+proj_pn(2);-ls+proj_pn(2);ls+proj_pn(2);ls+proj_pn(2)]];
	// Triangle represented in triangle coordinate system 3D
	// (X,Y,Z) -> 2D (U,V,W). Here we assume the value in W as 0
	triangle2D = [[an(1);bn(1);cn(1)],[an(2);bn(2);cn(2)]];
	// Sutherland Hodgman clipping - 2D (U,V,W). Here, the value in W is assumed 0 for all points.
	clippedPolygon2D = sutherlandHodgman(triangle2D,squareSupport);
	*/
}

// 3D -> 2D
void MpsParticle::transform3Dto2D(double *P1, double *RM) {
	double A[3];
	// Point represented in triangle coordinate system (U,V,W). The value in W is the same for all points (3D -> 2D).
	// Trasnpose of RM
	for(unsigned int i = 0; i < 3; i++) {
		A[i] = RM[3*i]*P1[0]+RM[3*i+1]*P1[1]+RM[3*i+2]*P1[2];
	}
	for(unsigned int i = 0; i < 3; i++) {
		P1[i] = A[i];
	}
}

// 2D -> 3D
void MpsParticle::transform2Dto3D(double *P1, double *RM) {
	double A[3];
	// Point represented in coordinate system (X,Y,Z) (2D -> 3D).
	for(unsigned int i = 0; i < 3; i++) {
		A[i] = RM[i]*P1[0]+RM[i+3]*P1[1]+RM[i+6]*P1[2];
	}
	for(unsigned int i = 0; i < 3; i++) {
		P1[i] = A[i];
	}
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
// 			if(normalMod2 > 1.0e-8) {
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
// 		transformMatrix(V1, V2, V3, RM);

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

// Update velocity and positions
void MpsParticle::updateVelocityPosition2nd() {
	velMax = 0.0;						// Maximum flow velocity
	double auxiliar[5] = {1.2, -3.3, 4.3, -0.3, 5.6};

	// https://stackoverflow.com/questions/39989473/use-openmp-in-c11-to-find-the-maximum-of-the-calculated-values
#pragma omp parallel
{
	double local_vMax = 0.0;
#pragma omp for
	for(int i=0; i<numParticles; i++) {
		if(particleType[i] == fluid) {
			vel[i*3  ]+=acc[i*3  ]*timeStep;	vel[i*3+1]+=acc[i*3+1]*timeStep;	vel[i*3+2]+=acc[i*3+2]*timeStep;
			//if(particleType[i] == fluid) {
			pos[i*3  ]+=acc[i*3  ]*timeStep*timeStep;	pos[i*3+1]+=acc[i*3+1]*timeStep*timeStep;	pos[i*3+2]+=acc[i*3+2]*timeStep*timeStep;
			//}
			acc[i*3]=acc[i*3+1]=acc[i*3+2]=0.0;

			//pos[i*3  ]=Posk[i*3 ]+vel[i*3  ]*timeStep;	pos[i*3+1]=Posk[i*3+1]+vel[i*3+1]*timeStep;	pos[i*3+2]=Posk[i*3+2]+vel[i*3+2]*timeStep;
			//Posk[i*3  ]=pos[i*3  ];	Posk[i*3+1]=pos[i*3+1];	Posk[i*3+2]=pos[i*3+2];
			//Velk[i*3  ]=vel[i*3  ];	Velk[i*3+1]=vel[i*3+1];	Velk[i*3+2]=vel[i*3+2];

			//wallParticleForce1[i*3  ]=acc[i*3  ];	wallParticleForce1[i*3+1]=acc[i*3+1];	wallParticleForce1[i*3+2]=acc[i*3+2];
			//wallParticleForce2[i*3  ]=Acv[i*3  ];	wallParticleForce2[i*3+1]=Acv[i*3+1];	wallParticleForce2[i*3+2]=Acv[i*3+2];
			//acc[i*3]=acc[i*3+1]=acc[i*3+2]=0.0;
			//Acv[i*3]=Acv[i*3+1]=Acv[i*3+2]=0.0;

			double vMod2 = vel[i*3  ]*vel[i*3  ] + vel[i*3+1]*vel[i*3+1] + vel[i*3+2]*vel[i*3+2];
			if(vMod2 > local_vMax*local_vMax)
				local_vMax = sqrt(vMod2);
		}
	}

#pragma omp critical
	{
		if (local_vMax > velMax)
			velMax = local_vMax;
	}
}
	CFLcurrent = timeStep*velMax/partDist;
}

// Shifting technique
// Improvements for accuracy and stability in a weakly-compressible particle method
// https://www.sciencedirect.com/science/article/pii/S0045793016302250
void MpsParticle::calcShifting() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
	if(particleType[i] == fluid) {
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double duXi = 0.0;	double duYi = 0.0;	double duZi = 0.0;
		//double ni = 0.0;
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					double dst = sqrt(dstij2);
					double dw = delWeight(dst, reS, weightType);
					duXi += dw*(vel[j*3  ]-velXi);
					duYi += dw*(vel[j*3+1]-velYi);
					duZi += dw*(vel[j*3+2]-velZi);
					//double w = weight(dst, r, weightType);
					//ni += w;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

//		if(wallType == boundaryWallType::POLYGON) {
//			if(pndi[i] > pndThreshold*pndSmallZero)
//		}
//		else if(wallType == boundaryWallType::PARTICLE) 
//			if(pndi[i] > pndThreshold*pndSmallZero || numNeigh[i] > neighThreshold*numNeighZero)
//		}
//		if(pndSmall[i] > pndThreshold*pndSmallZero || numNeigh[i] > neighThreshold*numNeighZero)
		if(particleBC[i] == inner) {
			vel[i*3  ] -= coeffShifting1*duXi;
			vel[i*3+1] -= coeffShifting1*duYi;
			vel[i*3+2] -= coeffShifting1*duZi;
		}
		// else {
		// 	double Inn[9], duAux[3];
		// 	// I - nxn
		// 	Inn[0] = 1.0 - normal[i*3  ]*normal[i*3  ]; Inn[1] = 0.0 - normal[i*3  ]*normal[i*3+1]; Inn[2] = 0.0 - normal[i*3  ]*normal[i*3+2];
		// 	Inn[3] = 0.0 - normal[i*3+1]*normal[i*3  ]; Inn[4] = 1.0 - normal[i*3+1]*normal[i*3+1]; Inn[5] = 0.0 - normal[i*3+1]*normal[i*3+2];
		// 	Inn[6] = 0.0 - normal[i*3+2]*normal[i*3  ]; Inn[7] = 0.0 - normal[i*3+2]*normal[i*3+1]; Inn[8] = 1.0 - normal[i*3+2]*normal[i*3+2];
		// 	// (I - nxn)dr
		// 	duAux[0] = Inn[0]*duXi + Inn[1]*duYi + Inn[2]*duZi;
		// 	duAux[1] = Inn[3]*duXi + Inn[4]*duYi + Inn[5]*duZi;
		// 	duAux[2] = Inn[6]*duXi + Inn[7]*duYi + Inn[8]*duZi;
		// 	vel[i*3  ] -= coeffShifting1*duAux[0];
		// 	vel[i*3+1] -= coeffShifting1*duAux[1];
		// 	vel[i*3+2] -= coeffShifting1*duAux[2];
		// }
	}}
}

// normal vector on the fluid
// An accurate and stable multiphase moving particle semi-implicit method based on a corrective matrix for all particle interaction models
// https://onlinelibrary.wiley.com/doi/full/10.1002/nme.5844
void MpsParticle::calcNormalParticles() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double dr_ix = 0.0;	double dr_iy = 0.0;	double dr_iz = 0.0;
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					double dst = sqrt(dstij2);
					double wS = weight(dst, reS, weightType);
					dr_ix += v0ij*wS/dst;
					dr_iy += v1ij*wS/dst;
					dr_iz += v2ij*wS/dst;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		normal[i*3  ] = - dr_ix/pndSmallZero;
		normal[i*3+1] = - dr_iy/pndSmallZero;
		normal[i*3+2] = - dr_iz/pndSmallZero;
		
		/*
		if(wallType == boundaryWallType::PARTICLE)
		{
			// Normalize
			double norm2 = dr_ix*dr_ix + dr_iy*dr_iy + dr_iz*dr_iz;
			if(norm2 > 0.0) {
				double norm = sqrt(norm2);
				normal[i*3  ] /= norm;
				normal[i*3+1] /= norm;
				normal[i*3+2] /= norm;
			}
			else {
				normal[i*3  ] = 0.0;
				normal[i*3+1] = 0.0;
				normal[i*3+2] = 0.0;
			}
		}
		*/
	}
}

// normal vector on the fluid (Polygon wall)
void MpsParticle::calcWallNormalParticles() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
	//if(particleType[i] == fluid) {
	if(particleType[i] == fluid && particleNearWall[i] == true) {
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double dr_ix = 0.0;	double dr_iy = 0.0;	double dr_iz = 0.0;
		// Wall gradient Mitsume`s model
		double Rref_i[9], normaliw[3], normaliwSqrt;
		// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
		normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);

		if(normaliwSqrt > 1.0e-8) {
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
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i) {
					double dst = sqrt(dstimj2);
					double wS = weight(dst, reS, weightType);
					dr_ix += v0imj*wS/dst;
					dr_iy += v1imj*wS/dst;
					dr_iz += v2imj*wS/dst;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Add "i" contribution ("i" is a neighbor of "mirror i")
		double v0imi, v1imi, v2imi, dstimi2;
		sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
	  	
		if(dstimi2 < reS2) {
			double dst = sqrt(dstimi2);
			double wS = weight(dst, reS, weightType);
			dr_ix += v0imi*wS/dst;
			dr_iy += v1imi*wS/dst;
			dr_iz += v2imi*wS/dst;
		}

		double drx = Rref_i[0]*dr_ix + Rref_i[1]*dr_iy + Rref_i[2]*dr_iz;
		double dry = Rref_i[3]*dr_ix + Rref_i[4]*dr_iy + Rref_i[5]*dr_iz;
		double drz = Rref_i[6]*dr_ix + Rref_i[7]*dr_iy + Rref_i[8]*dr_iz;
		
		normal[i*3  ] += - drx/pndSmallZero;
		normal[i*3+1] += - dry/pndSmallZero;
		normal[i*3+2] += - drz/pndSmallZero;

/*
		// Normalize
		double norm2 = normal[i*3  ]*normal[i*3  ] + normal[i*3+1]*normal[i*3+1] + normal[i*3+2]*normal[i*3+2];
		if(norm2 > 0.0) {
			double norm = sqrt(norm2);
			normal[i*3  ] /= norm;
			normal[i*3+1] /= norm;
			normal[i*3+2] /= norm;
		}
		else {
			normal[i*3  ] = 0.0;
			normal[i*3+1] = 0.0;
			normal[i*3+2] = 0.0;
		}
*/	
	}}
}


// Shifting technique (Polygon wall)
void MpsParticle::calcWallShifting() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
	//if(particleType[i] == fluid) {
	if(particleType[i] == fluid && particleNearWall[i] == true) {
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double duXi = 0.0;	double duYi = 0.0;	double duZi = 0.0;

		// No-slip

		// Inverse matrix Rinv_i = - I
		double Rinv_i[9], normaliw[3], normalMod2;
		// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
		normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

		if(normalMod2 > 1.0e-8) {
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

		if(nearMeshType[i] == meshType::FORCED) {
			viwall[0] = velVWall[0];
			viwall[1] = velVWall[1];
			viwall[2] = velVWall[2];
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
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i) {
					double dst = sqrt(dstimj2);
					double dw = delWeight(dst, reS, weightType);
					duXi += dw*(vel[j*3  ]-velMirrorXi);
					duYi += dw*(vel[j*3+1]-velMirrorYi);
					duZi += dw*(vel[j*3+2]-velMirrorZi);
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// Add "i" contribution ("i" is a neighbor of "mirror i")
		double v0imi, v1imi, v2imi, dstimi2;
		sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
	  	
		if(dstimi2 < reS2) {
			double dst = sqrt(dstimi2);
			double dw = delWeight(dst, reS, weightType);
			duXi += dw*(velXi-velMirrorXi);
			duYi += dw*(velYi-velMirrorYi);
			duZi += dw*(velZi-velMirrorZi);
		}
	
//		if(wallType == boundaryWallType::POLYGON) {
//			if(pndi[i] > pndThreshold*pndSmallZero)
//		}
//		else if(wallType == boundaryWallType::PARTICLE) 
//			if(pndi[i] > pndThreshold*pndSmallZero || numNeigh[i] > neighThreshold*numNeighZero)
//		}
//		if(pndSmall[i] > pndThreshold*pndSmallZero || numNeigh[i] > neighThreshold*numNeighZero)
		if(particleBC[i] == inner) {
			double dux = Rinv_i[0]*duXi + Rinv_i[1]*duYi + Rinv_i[2]*duZi;
			double duy = Rinv_i[3]*duXi + Rinv_i[4]*duYi + Rinv_i[5]*duZi;
			double duz = Rinv_i[6]*duXi + Rinv_i[7]*duYi + Rinv_i[8]*duZi;
			vel[i*3  ] -= coeffShifting1*dux;
			vel[i*3+1] -= coeffShifting1*duy;
			vel[i*3+2] -= coeffShifting1*duz;
		}
		// else {
		// 	double Inn[9], duAux[3];
		// 	// I - nxn
		// 	Inn[0] = 1.0 - normal[i*3  ]*normal[i*3  ]; Inn[1] = 0.0 - normal[i*3  ]*normal[i*3+1]; Inn[2] = 0.0 - normal[i*3  ]*normal[i*3+2];
		// 	Inn[3] = 0.0 - normal[i*3+1]*normal[i*3  ]; Inn[4] = 1.0 - normal[i*3+1]*normal[i*3+1]; Inn[5] = 0.0 - normal[i*3+1]*normal[i*3+2];
		// 	Inn[6] = 0.0 - normal[i*3+2]*normal[i*3  ]; Inn[7] = 0.0 - normal[i*3+2]*normal[i*3+1]; Inn[8] = 1.0 - normal[i*3+2]*normal[i*3+2];
		// 	double dux = Rinv_i[0]*duXi + Rinv_i[1]*duYi + Rinv_i[2]*duZi;
		//	double duy = Rinv_i[3]*duXi + Rinv_i[4]*duYi + Rinv_i[5]*duZi;
		//	double duz = Rinv_i[6]*duXi + Rinv_i[7]*duYi + Rinv_i[8]*duZi;
		// 	// (I - nxn)dr
		// 	duAux[0] = Inn[0]*dux + Inn[1]*duy + Inn[2]*duz;
		// 	duAux[1] = Inn[3]*dux + Inn[4]*duy + Inn[5]*duz;
		// 	duAux[2] = Inn[6]*dux + Inn[7]*duy + Inn[8]*duz;
		// 	vel[i*3  ] -= coeffShifting1*duAux[0];
		// 	vel[i*3+1] -= coeffShifting1*duAux[1];
		// 	vel[i*3+2] -= coeffShifting1*duAux[2];
		// }
	}}
}

// Concentration and Gradient of concentration
void MpsParticle::calcConcAndConcGradient() {
	// Concentration
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
		concentration[i] = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					double dst = sqrt(dstij2);
					double wS = weight(dst, reS, weightType);
					concentration[i] += wS;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// Add PND due Wall polygon
		concentration[i] += pndWallContribution[i];

		concentration[i] /= pndSmallZero;
	}
	// Gradient of concentration
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
	if(particleType[i] == fluid) {
		gradConcentration[i*3  ] = 0.0;	gradConcentration[i*3+1] = 0.0;	gradConcentration[i*3+2] = 0.0;
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double duXi = 0.0;	double duYi = 0.0;	double duZi = 0.0;
		double conc_i = concentration[i];
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && (dstij2 < dstimj2 || wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					double dst = sqrt(dstij2);
					double wS = weight(dst, reS, weightType);
					gradConcentration[i*3  ] += (conc_i + concentration[j])*v0ij*wS/dstij2;
					gradConcentration[i*3+1] += (conc_i + concentration[j])*v1ij*wS/dstij2;
					gradConcentration[i*3+2] += (conc_i + concentration[j])*v2ij*wS/dstij2;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		//coeffPressGrad = -dim/pndGradientZero
		gradConcentration[3*i  ] *= -coeffPressGrad;
		gradConcentration[3*i+1] *= -coeffPressGrad;
		gradConcentration[3*i+2] *= -coeffPressGrad;

//		if(wallType == boundaryWallType::POLYGON) {
//			if(pndi[i] > pndThreshold*pndSmallZero)
//		}
//		else if(wallType == boundaryWallType::PARTICLE) 
//			if(pndi[i] > pndThreshold*pndSmallZero || numNeigh[i] > neighThreshold*numNeighZero)
//		}
//		if(pndSmall[i] > pndThreshold*pndSmallZero || numNeigh[i] > neighThreshold*numNeighZero)
		if(particleBC[i] == inner) {
			// coeffShifting2 = coefA*partDist*partDist*cflNumber*machNumber;	// Coefficient used to adjust velocity
			pos[i*3  ] -= coeffShifting2*gradConcentration[3*i  ];
			pos[i*3+1] -= coeffShifting2*gradConcentration[3*i+1];
			pos[i*3+2] -= coeffShifting2*gradConcentration[3*i+2];
		}
		/*
		else
		{
			double Inn[9], drAux[3];
			// I - nxn
			Inn[0] = 1.0 - normal[i*3  ]*normal[i*3  ]; Inn[1] = 0.0 - normal[i*3  ]*normal[i*3+1]; Inn[2] = 0.0 - normal[i*3  ]*normal[i*3+2];
			Inn[3] = 0.0 - normal[i*3+1]*normal[i*3  ]; Inn[4] = 1.0 - normal[i*3+1]*normal[i*3+1]; Inn[5] = 0.0 - normal[i*3+1]*normal[i*3+2];
			Inn[6] = 0.0 - normal[i*3+2]*normal[i*3  ]; Inn[7] = 0.0 - normal[i*3+2]*normal[i*3+1]; Inn[8] = 1.0 - normal[i*3+2]*normal[i*3+2];
			// (I - nxn)dr
			drAux[0] = Inn[0]*gradConcentration[3*i  ] + Inn[1]*gradConcentration[3*i+1] + Inn[2]*gradConcentration[3*i+2];
			drAux[1] = Inn[3]*gradConcentration[3*i  ] + Inn[4]*gradConcentration[3*i+1] + Inn[5]*gradConcentration[3*i+2];
			drAux[2] = Inn[6]*gradConcentration[3*i  ] + Inn[7]*gradConcentration[3*i+1] + Inn[8]*gradConcentration[3*i+2];
			// coeffShifting2 = coefA*partDist*partDist*cflNumber*machNumber;	// Coefficient used to adjust Velocity
			pos[i*3  ] -= coeffShifting2*drAux[0];
			pos[i*3+1] -= coeffShifting2*drAux[1];
			pos[i*3+2] -= coeffShifting2*drAux[2];
		}
		*/
		// else {
		// 	double Inn[9], duAux[3];
		// 	// I - nxn
		// 	Inn[0] = 1.0 - normal[i*3  ]*normal[i*3  ]; Inn[1] = 0.0 - normal[i*3  ]*normal[i*3+1]; Inn[2] = 0.0 - normal[i*3  ]*normal[i*3+2];
		// 	Inn[3] = 0.0 - normal[i*3+1]*normal[i*3  ]; Inn[4] = 1.0 - normal[i*3+1]*normal[i*3+1]; Inn[5] = 0.0 - normal[i*3+1]*normal[i*3+2];
		// 	Inn[6] = 0.0 - normal[i*3+2]*normal[i*3  ]; Inn[7] = 0.0 - normal[i*3+2]*normal[i*3+1]; Inn[8] = 1.0 - normal[i*3+2]*normal[i*3+2];
		// 	// (I - nxn)dr
		// 	duAux[0] = Inn[0]*duXi + Inn[1]*duYi + Inn[2]*duZi;
		// 	duAux[1] = Inn[3]*duXi + Inn[4]*duYi + Inn[5]*duZi;
		// 	duAux[2] = Inn[6]*duXi + Inn[7]*duYi + Inn[8]*duZi;
		// 	vel[i*3  ] -= coeffShifting1*duAux[0];
		// 	vel[i*3+1] -= coeffShifting1*duAux[1];
		// 	vel[i*3+2] -= coeffShifting1*duAux[2];
		// }
	}}
}

// Concentration and Gradient of concentration (Polygon wall)
void MpsParticle::calcWallConcAndConcGradient() {
	// Gradient of concentration due Polygon wall
	//int nPartNearMesh = partNearMesh.size();
	double VolumeForce = pow(partDist,dim);
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<numParticles; i++) {
	//if(particleType[i] == fluid) {
	if(particleType[i] == fluid && particleNearWall[i] == true) {
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double posMirrorXi = mirrorParticlePos[i*3  ];	double posMirrorYi = mirrorParticlePos[i*3+1];	double posMirrorZi = mirrorParticlePos[i*3+2];
		double drX = 0.0;			double drY = 0.0;			double drZ = 0.0;
		double gradCiWallX = 0.0;			double gradCiWallY = 0.0;			double gradCiWallZ = 0.0;
		double ni = pndi[i];
		double conc_i = concentration[i];

		// Wall gradient Mitsume`s model
		double Rref_i[9], normaliw[3], normaliwSqrt;
		// normal fluid-wall particle = 0.5*(normal fluid-mirror particle)
		normaliw[0] = 0.5*(posXi - posMirrorXi); normaliw[1] = 0.5*(posYi - posMirrorYi); normaliw[2] = 0.5*(posZi - posMirrorZi);
		normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);

		if(normaliwSqrt > 1.0e-8) {
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
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);
				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2);

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i) {
					double dst = sqrt(dstimj2);
					double wS = weightGradient(dst, reS, weightType);

					drX += (conc_i + concentration[j])*v0imj*wS/dstimj2;
					drY += (conc_i + concentration[j])*v1imj*wS/dstimj2;
					drZ += (conc_i + concentration[j])*v2imj*wS/dstimj2;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		// Add "i" contribution ("i" is a neighbor of "mirror i")
		double v0imi, v1imi, v2imi, dstimi2;
		sqrDistBetweenParticles(i, posMirrorXi, posMirrorYi, posMirrorZi, v0imi, v1imi, v2imi, dstimi2);
		
		if(dstimi2 < reS2) {
			double dst = sqrt(dstimi2);
			double wS = weightGradient(dst, reS, weightType);

			drX += (conc_i + conc_i)*v0imi*wS/dstimi2;
			drY += (conc_i + conc_i)*v1imi*wS/dstimi2;
			drZ += (conc_i + conc_i)*v2imi*wS/dstimi2;
		}

		// coeffPressGrad is a negative cte (-dim/noGrad)
		// Original
//		acc[i*3  ] += (relaxPress*(Rref_i[0]*accX + Rref_i[1]*accY + Rref_i[2]*accZ)*coeffPressGrad - rpsForce[0])*invDns[partType::FLUID];
//		acc[i*3+1] += (relaxPress*(Rref_i[3]*accX + Rref_i[4]*accY + Rref_i[5]*accZ)*coeffPressGrad - rpsForce[1])*invDns[partType::FLUID];
//		acc[i*3+2] += (relaxPress*(Rref_i[6]*accX + Rref_i[7]*accY + Rref_i[8]*accZ)*coeffPressGrad - rpsForce[2])*invDns[partType::FLUID];
		// Modified
		gradCiWallX = -(Rref_i[0]*drX + Rref_i[1]*drY + Rref_i[2]*drZ)*coeffPressGrad;
		gradCiWallY = -(Rref_i[3]*drX + Rref_i[4]*drY + Rref_i[5]*drZ)*coeffPressGrad;
		gradCiWallZ = -(Rref_i[6]*drX + Rref_i[7]*drY + Rref_i[8]*drZ)*coeffPressGrad;

		gradConcentration[i*3  ] += gradCiWallX;
		gradConcentration[i*3+1] += gradCiWallY;
		gradConcentration[i*3+2] += gradCiWallZ;

		if(particleBC[i] == inner) {
			// coeffShifting2 = coefA*partDist*partDist*cflNumber*machNumber;	// Coefficient used to adjust velocity
			pos[i*3  ] -= coeffShifting2*gradCiWallX;
			pos[i*3+1] -= coeffShifting2*gradCiWallY;
			pos[i*3+2] -= coeffShifting2*gradCiWallZ;
		}
	}}
}

// normal vector on the fluid
void MpsParticle::calcNormalConcentration() {
#pragma omp parallel for
	for(int i=0; i<numParticles; i++) {
		double norm2GradCi = gradConcentration[3*i]*gradConcentration[3*i] + gradConcentration[3*i+1]*gradConcentration[3*i+1] + gradConcentration[3*i+2]*gradConcentration[3*i+2];
		
		if(norm2GradCi > 0.0) {
			double norm = sqrt(norm2GradCi);
			normal[i*3  ] = -gradConcentration[i*3  ]/norm;
			normal[i*3+1] = -gradConcentration[i*3+1]/norm;
			normal[i*3+2] = -gradConcentration[i*3+2]/norm;
		}
		else {
			normal[i*3  ] = 0.0;
			normal[i*3+1] = 0.0;
			normal[i*3+2] = 0.0;
		}
	}
}

// Update velocity at wall and dummy particles
void MpsParticle::updateVelocityParticlesWallDummy() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<numParticles; i++) {
	if(particleType[i] == wall) {
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		double velXi = vel[i*3  ];	double velYi = vel[i*3+1];	double velZi = vel[i*3+2];
		double duXi = 0.0;	double duYi = 0.0;	double duZi = 0.0;
		double ni = 0.0;
		
		int ix, iy, iz;
		bucketCoordinates(ix, iy, iz, posXi, posYi, posZi);
		int minZ = (iz-1)*((int)(dim-2.0)); int maxZ = (iz+1)*((int)(dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*numBucketsXY + jy*numBucketsX + jx;
			int j = firstParticleInBucket[jb];
			if(j == -1) continue;
			while(true) {
				double v0ij, v1ij, v2ij, dstij2;
				
				// Particle distance r_ij = Xj - Xi_temporary_position
				sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2);

				if(dstij2 < reS2) {
//				if(j != i) {
				if(j != i && particleType[j] == fluid) {
					double dst = sqrt(dstij2);
					double wS = weight(dst, reS, weightType);
					ni += wS;
					duXi += vel[j*3  ]*wS;
					duYi += vel[j*3+1]*wS;
					duZi += vel[j*3+2]*wS;
				}}
				j = nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}
		if(ni > 0.0) {
			vel[i*3  ] = 2.0*vel[i*3  ] - duXi/ni;
			vel[i*3+1] = 2.0*vel[i*3+1] - duYi/ni;
			vel[i*3+2] = 2.0*vel[i*3+2] - duZi/ni;
		}
		else {
			vel[i*3  ] = 2.0*vel[i*3  ] - duXi;
			vel[i*3+1] = 2.0*vel[i*3+1] - duYi;
			vel[i*3+2] = 2.0*vel[i*3+2] - duZi;
		}
	}}
}

// Pressure sensors
void MpsParticle::writePressSensors() {

	double P1, P2, P3, P4, posP1[3], posP2[3], posP3[3], posP4[3];
	double pndP1, pndP2, pndP3, pndP4, P1wij, P2wij, P3wij, P4wij;
	double riP1, riP2, riP3, riP4;

	double xPrs1Min, xPrs1Max, xPrs2Min, xPrs2Max;
	double yPrs1Min, yPrs1Max, yPrs2Min, yPrs2Max;
	double zPrs1min, zPrs2min, zPrs3min, zPrs4min;
	double zPrs1max, zPrs2max, zPrs3max, zPrs4max;

	// Pressure at sensor
	P1=P2=P3=P4=0.0;
	// PND for each sensor
	pndP1=pndP2=pndP3=pndP4=0.0;
	// Pressure * wij
	P1wij=P2wij=P3wij=P4wij=0.0;
	// Square distance between the fluid particle and the sensor
	riP1=riP2=riP3=riP4=10.0e10;

	// Sensor limits
	// Dam 1610
	/*
	xPrs1Min = 1.61 - 2.0*partDist; xPrs1Max = 1.61 + partDist;
	yPrs1Min = 0.075 - 2.0*partDist; yPrs1Max = 0.075 + 2.0*partDist;
	zPrs1min = 0.003 - partDist; zPrs2min = 0.015 - partDist; zPrs3min = 0.030 - partDist; zPrs4min = 0.080 - partDist;
	zPrs1max = 0.003 + partDist; zPrs2max = 0.015 + partDist; zPrs3max = 0.030 + partDist; zPrs4max = 0.080 + partDist;
*/
	// Hydrostatic
	xPrs1Min = 0.1 - 2.0*partDist; xPrs1Max = 0.1 + partDist; xPrs2Min = 0.0 - 2.0*partDist; xPrs2Max = 0.0 + partDist;
	yPrs1Min = 0.1 - 2.0*partDist; yPrs1Max = 0.1 + 2.0*partDist; yPrs2Min = 0.1 - 2.0*partDist; yPrs2Max = 0.1 + 2.0*partDist;
	zPrs1min = 0.0 - partDist; zPrs2min = 0.1 - partDist;
	zPrs1max = 0.0 + partDist; zPrs2max = 0.1 + partDist;

	// Sensor positions
	// Dam 1610
	/*
	posP1[0] = 1.61; posP1[1] = 0.075; posP1[2] = 0.003;
	posP2[0] = 1.61; posP2[1] = 0.075; posP2[2] = 0.015;
	posP3[0] = 1.61; posP3[1] = 0.075; posP3[2] = 0.030;
	posP4[0] = 1.61; posP4[1] = 0.075; posP4[2] = 0.080;
*/
	// Hydrostatic
	posP1[0] = 0.1; posP1[1] = 0.1; posP1[2] = 0.0;
	posP2[0] = 0.0; posP2[1] = 0.1; posP2[2] = 0.1;

#pragma omp parallel for reduction(+: pndP1, pndP2, pndP3, pndP4, P1wij, P2wij, P3wij, P4wij)
	for(int i=0; i<numParticles; i++) {
		double posXi = pos[i*3  ];	double posYi = pos[i*3+1];	double posZi = pos[i*3+2];
		
		/*
		// Pressure at a specific particle close to the sensor
		// Dam 1610
		if(posXi >= xPrs1Min && posXi <= xPrs1Max && posYi >= yPrs1Min && posYi <= yPrs1Max) {
			// Sensor 1
			if(posZi >= zPrs1min && posZi <= zPrs1max) {
				if(press[i] > 0.0) {
					double v0 = posP1[0] - posXi;
					double v1 = posP1[1] - posYi;
					double v2 = posP1[2] - posZi;
					double dst2 = v0*v0+v1*v1+v2*v2;
					// Closest fluid particle
					if(dst2 < riP1) {
						P1 = press[i];
						riP1 = dst2;
					}
				}
			}
			// Sensor 2
			if(posZi >= zPrs2min && posZi <= zPrs2max) {
				if(press[i] > 0.0) {
					double v0 = posP2[0] - posXi;
					double v1 = posP2[1] - posYi;
					double v2 = posP2[2] - posZi;
					double dst2 = v0*v0+v1*v1+v2*v2;
					// Closest fluid particle
					if(dst2 < riP2) {
						P2 = press[i];
						riP2 = dst2;
					}
				}
			}
			// Sensor 3
			if(posZi >= zPrs3min && posZi <= zPrs3max) {
				if(press[i] > 0.0) {
					double v0 = posP3[0] - posXi;
					double v1 = posP3[1] - posYi;
					double v2 = posP3[2] - posZi;
					double dst2 = v0*v0+v1*v1+v2*v2;
					// Closest fluid particle
					if(dst2 < riP3) {
						P3 = press[i];
						riP3 = dst2;
					}
				}
			}
			// Sensor 4
			if(posZi >= zPrs4min && posZi <= zPrs4max) {
				if(press[i] > 0.0) {
					double v0 = posP4[0] - posXi;
					double v1 = posP4[1] - posYi;
					double v2 = posP4[2] - posZi;
					double dst2 = v0*v0+v1*v1+v2*v2;
					// Closest fluid particle
					if(dst2 < riP4) {
						P4 = press[i];
						riP4 = dst2;
					}
				}
			}
		}

		// Weighted average pressure
		double v0 = posP1[0] - posXi;
		double v1 = posP1[1] - posYi;
		double v2 = posP1[2] - posZi;
		double dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2 < reS2) {
			double dst = sqrt(dst2);
			double wS = weight(dst, reS, weightType);
			if(press[i] > 0.0){
				pndP1 += wS;
				P1wij += press[i]*wS;
			}
		}
		
		v0 = posP2[0] - posXi;
		v1 = posP2[1] - posYi;
		v2 = posP2[2] - posZi;
		dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2 < reS2) {
			double dst = sqrt(dst2);
			double wS = weight(dst, reS, weightType);
			if(press[i] > 0.0){
				pndP2 += wS;
				P2wij += press[i]*wS;
			}
		}
		
		v0 = posP3[0] - posXi;
		v1 = posP3[1] - posYi;
		v2 = posP3[2] - posZi;
		dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2 < reS2) {
			double dst = sqrt(dst2);
			double wS = weight(dst, reS, weightType);
			if(press[i] > 0.0){
				pndP3 += wS;
				P3wij += press[i]*wS;
			}
		}
		
		v0 = posP4[0] - posXi;
		v1 = posP4[1] - posYi;
		v2 = posP4[2] - posZi;
		dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2 < reS2) {
			double dst = sqrt(dst2);
			double wS = weight(dst, reS, weightType);
			if(press[i] > 0.0){
				pndP4 += wS;
				P4wij += press[i]*wS;
			}
		}
		*/
		// Pressure at a specific particle close to the sensor
		// Hydrostatic
		if((posXi >= xPrs1Min && posXi <= xPrs1Max && posYi >= yPrs1Min && posYi <= yPrs1Max) || (posXi >= xPrs2Min && posXi <= xPrs2Max && posYi >= yPrs2Min && posYi <= yPrs2Max)) {
			// Sensor 1
			if(posZi >= zPrs1min && posZi <= zPrs1max) {
				if(press[i] > 0.0) {
					double v0 = posP1[0] - posXi;
					double v1 = posP1[1] - posYi;
					double v2 = posP1[2] - posZi;
					double dst2 = v0*v0+v1*v1+v2*v2;
					// Closest fluid particle
					if(dst2 < riP1) {
						P1 = press[i];
						riP1 = dst2;
					}
				}
			}
			// Sensor 2
			if(posZi >= zPrs2min && posZi <= zPrs2max) {
				if(press[i] > 0.0) {
					double v0 = posP2[0] - posXi;
					double v1 = posP2[1] - posYi;
					double v2 = posP2[2] - posZi;
					double dst2 = v0*v0+v1*v1+v2*v2;
					// Closest fluid particle
					if(dst2 < riP2) {
						P2 = press[i];
						riP2 = dst2;
					}
				}
			}
			
		}

		// Weighted average pressure
		double v0 = posP1[0] - posXi;
		double v1 = posP1[1] - posYi;
		double v2 = posP1[2] - posZi;
		double dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2 < reS2) {
			double dst = sqrt(dst2);
			double wS = weight(dst, reS, weightType);
			if(press[i] > 0.0){
				pndP1 += wS;
				P1wij += press[i]*wS;
			}
		}
		
		v0 = posP2[0] - posXi;
		v1 = posP2[1] - posYi;
		v2 = posP2[2] - posZi;
		dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2 < reS2) {
			double dst = sqrt(dst2);
			double wS = weight(dst, reS, weightType);
			if(press[i] > 0.0){
				pndP2 += wS;
				P2wij += press[i]*wS;
			}
		}

	}

	char *pressTxtFilenameChar = new char[pressTxtFilename.length()+1];
	strcpy (pressTxtFilenameChar, pressTxtFilename.c_str());
	// pressTxtFilenameChar now contains a c-string copy of pressTxtFilename

	// Open the File to write
	pressTxtFile = fopen(pressTxtFilenameChar, "a");
	if(pressTxtFile == NULL) perror ("Error opening press txt file");

	delete[] pressTxtFilenameChar;
	pressTxtFilenameChar = NULL;
	/*
	// dam 1610
	fprintf(pressTxtFile,"\n%lf\t%lf\t%lf\t%lf\t%lf",timeCurrent, P1,P2,P3,P4);
	double Pmean;
	
	if(pndP1 > 0.0)
		Pmean = P1wij/pndP1;
	else
		Pmean = 0.0;
	fprintf(pressTxtFile,"\t%lf",Pmean);
	
	if(pndP2 > 0.0)
		Pmean = P2wij/pndP2;
	else
		Pmean = 0.0;
	fprintf(pressTxtFile,"\t%lf",Pmean);
	
	if(pndP3 > 0.0)
		Pmean = P3wij/pndP3;
	else
		Pmean = 0.0;
	fprintf(pressTxtFile,"\t%lf",Pmean);
	
	if(pndP4 > 0.0)
		Pmean = P4wij/pndP4;
	else
		Pmean = 0.0;
	fprintf(pressTxtFile,"\t%lf",Pmean);
	*/

	// Hydrostatic
	fprintf(pressTxtFile,"\n%lf\t%lf\t%lf",timeCurrent, P1,P2);
	double Pmean;
	
	if(pndP1 > 0.0)
		Pmean = P1wij/pndP1;
	else
		Pmean = 0.0;
	fprintf(pressTxtFile,"\t%lf",Pmean);
	
	if(pndP2 > 0.0)
		Pmean = P2wij/pndP2;
	else
		Pmean = 0.0;
	fprintf(pressTxtFile,"\t%lf",Pmean);
	

	// Close press file
	fclose(pressTxtFile);

}

void MpsParticle::writeHeaderTxtFiles() {
	if(txtForce ==  true) {
		char *forceTxtFilenameChar = new char[forceTxtFilename.length()+1];
		strcpy (forceTxtFilenameChar, forceTxtFilename.c_str());
		// forceTxtFilenameChar now contains a c-string copy of forceTxtFilename
		
		// Open the File to write
		forceTxtFile = fopen(forceTxtFilenameChar, "w");
		if(forceTxtFile == NULL) perror ("Error opening force txt file");
		fprintf(forceTxtFile,"Time(s)\tForce(N) X\tY\tZ");
		// Close force file
		fclose(forceTxtFile);

		delete[] forceTxtFilenameChar;
		forceTxtFilenameChar = NULL;
	}
	if(txtPress ==  true) {
		char *pressTxtFilenameChar = new char[pressTxtFilename.length()+1];
		strcpy (pressTxtFilenameChar, pressTxtFilename.c_str());
		// pressTxtFilenameChar now contains a c-string copy of pressTxtFilename

		// Open the File to write
		pressTxtFile = fopen(pressTxtFilenameChar, "w");
		if(pressTxtFile == NULL) perror ("Error opening press txt file");

		// dam 1610
		// fprintf(pressTxtFile,"Time(s)\tP1(Pa)\tP2(Pa)\tP3(Pa)\tP4(Pa)\tPm1(Pa)\tPm2(Pa)\tPm3(Pa)\tPm4(Pa)");

		// Hydrostatic
		fprintf(pressTxtFile,"Time(s)\tP1(Pa)\tP2(Pa)\tPm1(Pa)\tPm2(Pa)");

		// Close press file
		fclose(pressTxtFile);

		delete[] pressTxtFilenameChar;
		pressTxtFilenameChar = NULL;
	}
}

// Call functions to write output files
void MpsParticle::writeOutputFiles() {
	// writeProfAscii();
	// Write particle data (VTU files)
	if(vtuType == 0) {
		if(freeSurfWall == true) {
			writeVtuAsciiFreeSurface();
		}
		else {
			writeVtuAscii();
		}
	}
	else {
		if(freeSurfWall == true) {
			writeVtuBinaryFreeSurface();
		}
		else {
			writeVtuBinary();
		}
	}
}

// Write data. Format .prof
void MpsParticle::writeProfAscii()
{
	char output_filename[256];
	sprintf(output_filename, "output%05d.prof",fileNumber);
	fp = fopen(output_filename, "w");
	fprintf(fp,"%d\n",numParticles);
	for(int i=0; i<numParticles; i++)
	{
		int a[2];
		double b[9];
		a[0]=i;	a[1]=particleType[i];
		b[0]=pos[i*3];	b[1]=pos[i*3+1];	b[2]=pos[i*3+2];
		b[3]=vel[i*3];	b[4]=vel[i*3+1];	b[5]=vel[i*3+2];
		b[6]=press[i];		b[7]=pressAverage[i]/iterOutput;
		b[8]=pndi[i];
		fprintf(fp," %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",a[0],a[1],b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[8]);
		pressAverage[i]=0.0;
	}
	// Close .prof file
	fclose(fp);
}

// https://stackoverflow.com/questions/105252
// https://stackoverflow.com/questions/10913666/error-writing-binary-vtk-files
// https://stackoverflow.com/questions/55829282/write-vtk-file-in-binary-format
template <typename T> void MpsParticle::SwapEnd(T& var)
{
	char* varArray = reinterpret_cast<char*>(&var);
	for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
		swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}

// Write data. Format .vtu (Paraview)
void MpsParticle::writeVtuBinary()
{
	char output_filename[256];
	char *output_folder_char = new char[vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, vtuOutputFoldername.c_str());
	// output_folder_char now contains a c-string copy of vtuOutputFoldername

	//char output_folder_char[vtuOutputFoldername.size() + 1];
    //strcpy(output_folder_char, vtuOutputFoldername.c_str());    // or pass &vtuOutputFoldername[0]

    //string aux_filename = path + "/mesh" + mesh_IDstr + "_" + numstr + ".stl";

	sprintf(output_filename, "%s/output%05d.vtu",output_folder_char,fileNumber);

	delete[] output_folder_char;
	output_folder_char = NULL;

	// BINARY FILE
	ofstream file;
	file.open(output_filename, ios::out | ios::binary);

	int nParticles = numParticles;

	// Header
	//file << "<?xml version='1.0' encoding='UTF-8'?>" << endl;
	file << "<VTKFile type='UnstructuredGrid' version='1.0' byte_order='LittleEndian' header_type='UInt64'>" << endl;
	file << "  <UnstructuredGrid>" << endl;
	file << "    <Piece NumberOfPoints='" << nParticles << "' NumberOfCells='" << nParticles << "'>" << endl;
	
	// Point data
	file << "      <PointData>" << endl;
	file << "        <DataArray type='Float32' Name='Velocity' NumberOfComponents='3' format='binary'>" << endl;
	for(size_t i=0; i<numParticles; i++)
	{
		{
			float ptx = (float)vel[i*3  ];
			float pty = (float)vel[i*3+1];
			float ptz = (float)vel[i*3+2];

			SwapEnd(ptx);
			SwapEnd(pty);
			SwapEnd(ptz);

			file.write(reinterpret_cast<char*>(&ptx), sizeof(float));
			file.write(reinterpret_cast<char*>(&pty), sizeof(float));
			file.write(reinterpret_cast<char*>(&ptz), sizeof(float));
		}
	}
	file << endl << "        </DataArray>" << endl;
	
	file << "        <DataArray type='Float32' Name='PreSmallsure' format='binary'>" << endl;
	for(size_t i=0; i<numParticles; i++)
	{
		{
			float ptv = (float)press[i];
			SwapEnd(ptv);
			file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
		}
	}
	file << endl << "        </DataArray>" << endl;
	
	file << "<        DataArray type='Int32' Name='BC' format='binary'>" << endl;
	for(size_t i=0; i<numParticles; i++)
	{
		{
			int32_t type = particleBC[i];
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
		}
	}
	file << endl << "        </DataArray>" << endl;

	if(outputPnd) 
	{
		file << "        <DataArray type='Float32' Name='pnd' format='binary'>" << endl;
		for(size_t i=0; i<numParticles; i++)
		{
			{
				float ptv = pndi[i];
				SwapEnd(ptv);
				file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
			}
		}
		file << endl << "        </DataArray>" << endl;
		
		file << "        <DataArray type='Float32' Name='pndSmall' format='binary'>" << endl;
		for(size_t i=0; i<numParticles; i++)
		{
			{
				float ptv = pndSmall[i];
				SwapEnd(ptv);
				file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
			}
		}
		file << endl << "        </DataArray>" << endl;
	}

	if(outputNeigh)
	{
		file << "        <DataArray type='Int32' Name='nNeigh' format='binary'>" << endl;
		for(size_t i=0; i<numParticles; i++)
		{
			{
				int32_t type = numNeigh[i];
				int type_i = static_cast<int>(type);
				SwapEnd(type_i);
				file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
			}
		}
		file << endl << "        </DataArray>" << endl;
	}

	file << "      </PointData>" << endl;

	// Cell data
	file << "      <CellData>" << endl;
	file << "      </CellData>" << endl;

	//Points
	file << "      <Points>" << endl;
	file << "        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='binary'>" << endl;
	for(size_t i=0; i<numParticles; i++)
	{
		{
			float ptx = (float)pos[i*3  ];
			float pty = (float)pos[i*3+1];
			float ptz = (float)pos[i*3+2];

			SwapEnd(ptx);
			SwapEnd(pty);
			SwapEnd(ptz);

			file.write(reinterpret_cast<char*>(&ptx), sizeof(float));
			file.write(reinterpret_cast<char*>(&pty), sizeof(float));
			file.write(reinterpret_cast<char*>(&ptz), sizeof(float));
		}
	}
	file << endl <<"        </DataArray>" << endl;
	file << "      </Points>" << endl;

	// Cells
	file << "      <Cells>" << endl;
	file << "        <DataArray type='Int64' Name='connectivity' format='binary'>" << endl;
	for(size_t i=0, ii=0; i<numParticles; i++)
	{
		{
			int64_t type = ii;
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
			ii++;
		}
	}
	file << endl << "        </DataArray>" << endl;
	
	file << "        <DataArray type='Int64' Name='offsets' format='binary'>" << endl;
	for(size_t i=0, ii=0; i<numParticles; i++)
	{
		{
			int64_t type = ii+1;
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
			ii++;
		}
	}
	file << endl << "        </DataArray>" << endl;

	file << "        <DataArray type='UInt8' Name='types' format='binary'>" << endl;
	for(size_t i=0; i<numParticles; i++)
	{
		{
			uint8_t type = 1;
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
		}
	}
	file << endl << "        </DataArray>" << endl;
	file << "      </Cells>" << endl;
	file << "    </Piece>" << endl;
	file << "  </UnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;

	file.close();
}

// Write data. Format .vtu (Paraview). Print only free-surface particles
void MpsParticle::writeVtuBinaryFreeSurface()
{

	char output_filename[256];
	char *output_folder_char = new char[vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, vtuOutputFoldername.c_str());
	// output_folder_char now contains a c-string copy of vtuOutputFoldername

	sprintf(output_filename, "%s/output%05d.vtu",output_folder_char,fileNumber);

	delete[] output_folder_char;
	output_folder_char = NULL;

	int nParticles = 0;

	for(int i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
			nParticles++;
	}

	// BINARY FILE
	ofstream file;
	file.open(output_filename, ios::out | ios::binary);

	// Header
	//file << "<?xml version='1.0' encoding='UTF-8'?>" << endl;
	file << "<VTKFile type='UnstructuredGrid' version='1.0' byte_order='LittleEndian' header_type='UInt64'>" << endl;
	file << "  <UnstructuredGrid>" << endl;
	file << "    <Piece NumberOfPoints='" << nParticles << "' NumberOfCells='" << nParticles << "'>" << endl;
	
	// Point data
	file << "      <PointData>" << endl;
	file << "        <DataArray type='Float32' Name='Velocity' NumberOfComponents='3' format='binary'>" << endl;
	for(size_t i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			float ptx = (float)vel[i*3  ];
			float pty = (float)vel[i*3+1];
			float ptz = (float)vel[i*3+2];

			SwapEnd(ptx);
			SwapEnd(pty);
			SwapEnd(ptz);

			file.write(reinterpret_cast<char*>(&ptx), sizeof(float));
			file.write(reinterpret_cast<char*>(&pty), sizeof(float));
			file.write(reinterpret_cast<char*>(&ptz), sizeof(float));
		}
	}
	file << endl << "        </DataArray>" << endl;
	
	file << "        <DataArray type='Float32' Name='PreSmallsure' format='binary'>" << endl;
	for(size_t i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			float ptv = (float)press[i];
			SwapEnd(ptv);
			file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
		}
	}
	file << endl << "        </DataArray>" << endl;
	
	file << "<        DataArray type='Int32' Name='BC' format='binary'>" << endl;
	for(size_t i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			int32_t type = particleBC[i];
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
		}
	}
	file << endl << "        </DataArray>" << endl;

	if(outputPnd)
	{
		file << "        <DataArray type='Float32' Name='pnd' format='binary'>" << endl;
		for(size_t i=0; i<numParticles; i++)
		{
			if(particleBC[i] == surface || particleType[i] == wall)
			{
				float ptv = pndi[i];
				SwapEnd(ptv);
				file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
			}
		}
		file << endl << "        </DataArray>" << endl;
		
		file << "        <DataArray type='Float32' Name='pndSmall' format='binary'>" << endl;
		for(size_t i=0; i<numParticles; i++)
		{
			if(particleBC[i] == surface || particleType[i] == wall)
			{
				float ptv = pndSmall[i];
				SwapEnd(ptv);
				file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
			}
		}
		file << endl << "        </DataArray>" << endl;
		
	}

	if(outputNeigh)
	{
		file << "        <DataArray type='Int32' Name='nNeigh' format='binary'>" << endl;
		for(size_t i=0; i<numParticles; i++)
		{
			if(particleBC[i] == surface || particleType[i] == wall)
			{
				int32_t type = numNeigh[i];
				int type_i = static_cast<int>(type);
				SwapEnd(type_i);
				file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
			}
		}
		file << endl << "        </DataArray>" << endl;
		
	}

	file << "      </PointData>" << endl;

	// Cell data
	file << "      <CellData>" << endl;
	file << "      </CellData>" << endl;

	//Points
	file << "      <Points>" << endl;
	file << "        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='binary'>" << endl;
	for(size_t i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			float ptx = (float)pos[i*3  ];
			float pty = (float)pos[i*3+1];
			float ptz = (float)pos[i*3+2];

			SwapEnd(ptx);
			SwapEnd(pty);
			SwapEnd(ptz);

			file.write(reinterpret_cast<char*>(&ptx), sizeof(float));
			file.write(reinterpret_cast<char*>(&pty), sizeof(float));
			file.write(reinterpret_cast<char*>(&ptz), sizeof(float));
		}
	}
	file << endl <<"        </DataArray>" << endl;
	file << "      </Points>" << endl;

	// Cells
	file << "      <Cells>" << endl;
	file << "        <DataArray type='Int64' Name='connectivity' format='binary'>" << endl;
	for(size_t i=0, ii=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			int64_t type = ii;
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
			ii++;
		}
	}
	file << endl << "        </DataArray>" << endl;
	
	file << "        <DataArray type='Int64' Name='offsets' format='binary'>" << endl;
	for(size_t i=0, ii=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			int64_t type = ii+1;
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
			ii++;
		}
	}
	file << endl << "        </DataArray>" << endl;

	file << "        <DataArray type='UInt8' Name='types' format='binary'>" << endl;
	for(size_t i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			uint8_t type = 1;
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
		}
	}
	file << endl << "        </DataArray>" << endl;
	file << "      </Cells>" << endl;
	file << "    </Piece>" << endl;
	file << "  </UnstructuredGrid>" << endl;
	file << "</VTKFile>" << endl;

	file.close();
}

// Write data. Format .vtu (Paraview)
void MpsParticle::writeVtuAscii()
{
	char output_filename[256];
	char *output_folder_char = new char[vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, vtuOutputFoldername.c_str());
	// output_folder_char now contains a c-string copy of vtuOutputFoldername

	sprintf(output_filename, "%s/output%05d.vtu",output_folder_char,fileNumber);

	delete[] output_folder_char;
	output_folder_char = NULL;

	int nParticles = numParticles;

	// ASCII FILE
	fp = fopen(output_filename, "w");

	// Header
	fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
	fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
	fprintf(fp,"  <UnstructuredGrid>\n");
	fprintf(fp,"    <Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",nParticles,nParticles);

	// Points
	fprintf(fp,"      <Points>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='ascii'>\n          ");
	for(int i=0; i<numParticles; i++)
	{
		fprintf(fp,"%lf %lf %lf ",pos[i*3],pos[i*3+1],pos[i*3+2]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"      </Points>\n");

	// Point data
	fprintf(fp,"      <PointData>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='Velocity' NumberOfComponents='3' format='ascii'>\n");
	for(int i=0; i<numParticles; i++)
	{
		//double val=sqrt(vel[i*3]*vel[i*3]+vel[i*3+1]*vel[i*3+1]+vel[i*3+2]*vel[i*3+2]);
		fprintf(fp,"%f %f %f ",(float)vel[i*3],(float)vel[i*3+1],(float)vel[i*3+2]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	//fprintf(fp,"        <DataArray type='Float32' Name='preSmallsave' format='ascii'>\n");
	//for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)(pressAverage[i]/iterOutput));}
	//fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='pressure' format='ascii'>\n");
	for(int i=0; i<numParticles; i++)
	{
		fprintf(fp,"%f ",(float)press[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='BC' format='ascii'>\n");
	for(int i=0; i<numParticles; i++)
	{
		fprintf(fp,"%d ",particleBC[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	
	if(outputPnd)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='pnd' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			fprintf(fp,"%f ",(float)pndi[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='pndSmall' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			fprintf(fp,"%f ",(float)pndSmall[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='pndk' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++)
//		{
//			fprintf(fp,"%f ",(float)pndki[i]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='pndsk' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++)
//		{
//			fprintf(fp,"%f ",(float)pndski[i]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");
	}
	if(outputNeigh) {
		fprintf(fp,"        <DataArray type='Int32' Name='nNeigh' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			fprintf(fp,"%d ",numNeigh[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}

	if(outputDeviation)
	{
//		fprintf(fp,"        <DataArray type='Float32' Name='devSquare' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)npcdDeviation2[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='npcdDeviation' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			fprintf(fp,"%f %f %f ",(float)npcdDeviation[i*3],(float)npcdDeviation[i*3+1],(float)npcdDeviation[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='polygonNormal' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			fprintf(fp,"%f %f %f ",(float)polygonNormal[i*3],(float)polygonNormal[i*3+1],(float)polygonNormal[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='deviationDotPolygonNormal' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)deviationDotPolygonNormal[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
	}

	if(outputConcentration)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='concentration' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			fprintf(fp,"%f ",(float)concentration[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='GradCi' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			fprintf(fp,"%f %f %f ",(float)gradConcentration[i*3],(float)gradConcentration[i*3+1],(float)gradConcentration[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}
	
//	fprintf(fp,"        <DataArray type='Float32' Name='wallParticleForce1' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<numParticles; i++) {
//		fprintf(fp,"%f %f %f ",(float)wallParticleForce1[i*3],(float)wallParticleForce1[i*3+1],(float)wallParticleForce1[i*3+2]);
//	}
//	fprintf(fp,"\n        </DataArray>\n");

//	fprintf(fp,"        <DataArray type='Float32' Name='wallParticleForce2' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<numParticles; i++) {
//		fprintf(fp,"%f %f %f ",(float)wallParticleForce2[i*3],(float)wallParticleForce2[i*3+1],(float)wallParticleForce2[i*3+2]);
//	}
//	fprintf(fp,"\n        </DataArray>\n");

//	fprintf(fp,"        <DataArray type='Float32' Name='Normal' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<numParticles; i++) {
//		fprintf(fp,"%f %f %f ",(float)normal[i*3],(float)normal[i*3+1],(float)normal[i*3+2]);
//	}
//	fprintf(fp,"\n        </DataArray>\n");
	
	if(outputNonNewtonian)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='RHO' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)RHO[i]);}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='ConcVol' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)Cv[i]);}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='MEU' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)MEU[i]);}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='MEUy' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)MEU_Y[i]);}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='II' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)II[i]);}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Int32' Name='PTYPE' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {	fprintf(fp,"%d ",PTYPE[i]);}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='p_smooth' format='ascii'>\n");
		for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)p_smooth[i]);}
		fprintf(fp,"\n        </DataArray>\n");
	}
	
	if(outputAuxiliar)
	{
		fprintf(fp,"        <DataArray type='Int32' Name='ParticleType' format='ascii'>\n");
		for(int i=0; i<numParticles; i++) {
			fprintf(fp,"%d ",particleType[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");

//		fprintf(fp,"        <DataArray type='Float32' Name='mirrorParticlePos' NumberOfComponents='3' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {
//			fprintf(fp,"%f %f %f ",(float)mirrorParticlePos[i*3],(float)mirrorParticlePos[i*3+1],(float)mirrorParticlePos[i*3+2]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

//		fprintf(fp,"        <DataArray type='Int32' Name='NearWall' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {
//			fprintf(fp,"%d ",particleNearWall[i]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

		//fprintf(fp,"        <DataArray type='Float32' Name='Fwall' NumberOfComponents='3' format='ascii'>\n");
		//for(int i=0; i<numParticles; i++) {
		//	fprintf(fp,"%f %f %f ",(float)Fwall[i*3],(float)Fwall[i*3+1],(float)Fwall[i*3+2]);
		//}
		//fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='Fwall' NumberOfComponents='3' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {
//			fprintf(fp,"%f %f %f ",(float)forceWall[i*3],(float)forceWall[i*3+1],(float)forceWall[i*3+2]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='DIV' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)velDivergence[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='Di' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)diffusiveTerm[i]);}
//		fprintf(fp,"\n        </DataArray>\n");

		fprintf(fp,"        <DataArray type='Int32' Name='nearMeshType' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			fprintf(fp,"%d ",nearMeshType[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");

		//fprintf(fp,"        <DataArray type='Float32' Name='numNeighborsSurfaceParticles' format='ascii'>\n");
		//for(int i=0; i<numParticles; i++) {fprintf(fp,"%f ",(float)numNeighborsSurfaceParticles[i]);}
		//fprintf(fp,"\n        </DataArray>\n");
	}

	fprintf(fp,"      </PointData>\n");

	// Cells
	fprintf(fp,"      <Cells>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='connectivity' format='ascii'>\n");
	for(int i=0, ii = 0; i<numParticles; i++)
	{
		fprintf(fp,"%d ",ii);
		ii++;
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='offsets' format='ascii'>\n");
	for(int i=0, ii=0; i<numParticles; i++)
	{
		fprintf(fp,"%d ",ii+1);
		ii++;
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='UInt8' Name='types' format='ascii'>\n");
	for(int i=0; i<numParticles; i++)
	{
		fprintf(fp,"1 ");
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"      </Cells>\n");
	fprintf(fp,"    </Piece>\n");
	fprintf(fp,"  </UnstructuredGrid>\n");
	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
}

// Write data. Format .vtu (Paraview). Print only free-surface particles
void MpsParticle::writeVtuAsciiFreeSurface()
{
	char output_filename[256];
	char *output_folder_char = new char[vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, vtuOutputFoldername.c_str());
	// output_folder_char now contains a c-string copy of vtuOutputFoldername

	sprintf(output_filename, "%s/output%05d.vtu",output_folder_char,fileNumber);

	delete[] output_folder_char;
	output_folder_char = NULL;

	// ASCII FILE
	fp = fopen(output_filename, "w");

	int nParticles = 0;
	for(int i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
			nParticles++;
	}

	// Header
	fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
	fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
	fprintf(fp,"  <UnstructuredGrid>\n");
	fprintf(fp,"    <Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",nParticles,nParticles);

	// Points
	fprintf(fp,"      <Points>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='ascii'>\n          ");
	for(int i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			fprintf(fp,"%lf %lf %lf ",pos[i*3],pos[i*3+1],pos[i*3+2]);
		}
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"      </Points>\n");

	// Point data
	fprintf(fp,"      <PointData>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='Velocity' NumberOfComponents='3' format='ascii'>\n");
	for(int i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
			fprintf(fp,"%f %f %f ",(float)vel[i*3],(float)vel[i*3+1],(float)vel[i*3+2]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	//fprintf(fp,"        <DataArray type='Float32' Name='preSmallsave' format='ascii'>\n");
	//for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)(pressAverage[i]/iterOutput));}
	//fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='pressure' format='ascii'>\n");
	for(int i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
			fprintf(fp,"%f ",(float)press[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='BC' format='ascii'>\n");
	for(int i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
			fprintf(fp,"%d ",particleBC[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	
	if(outputPnd)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='pnd' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			if(particleBC[i] == surface || particleType[i] == wall)
				fprintf(fp,"%f ",(float)pndi[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='pndSmall' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			if(particleBC[i] == surface || particleType[i] == wall)
				fprintf(fp,"%f ",(float)pndSmall[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}

	if(outputNeigh)
	{
		fprintf(fp,"        <DataArray type='Int32' Name='nNeigh' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			if(particleBC[i] == surface || particleType[i] == wall)
				fprintf(fp,"%d ",numNeigh[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}

	if(outputDeviation)
	{
//		fprintf(fp,"        <DataArray type='Float32' Name='devSquare' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)npcdDeviation2[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='npcdDeviation' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			if(particleBC[i] == surface || particleType[i] == wall)
				fprintf(fp,"%f %f %f ",(float)npcdDeviation[i*3],(float)npcdDeviation[i*3+1],(float)npcdDeviation[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='polygonNormal' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			fprintf(fp,"%f %f %f ",(float)polygonNormal[i*3],(float)polygonNormal[i*3+1],(float)polygonNormal[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='deviationDotPolygonNormal' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)deviationDotPolygonNormal[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
	}

	if(outputConcentration)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='concentration' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			if(particleBC[i] == surface || particleType[i] == wall)
				fprintf(fp,"%f ",(float)concentration[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='GradCi' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			if(particleBC[i] == surface || particleType[i] == wall)
				fprintf(fp,"%f %f %f ",(float)gradConcentration[i*3],(float)gradConcentration[i*3+1],(float)gradConcentration[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}
	
//	fprintf(fp,"        <DataArray type='Float32' Name='wallParticleForce1' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<numParticles; i++) {
//		fprintf(fp,"%f %f %f ",(float)wallParticleForce1[i*3],(float)wallParticleForce1[i*3+1],(float)wallParticleForce1[i*3+2]);
//	}
//	fprintf(fp,"\n        </DataArray>\n");

//	fprintf(fp,"        <DataArray type='Float32' Name='wallParticleForce2' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<numParticles; i++) {
//		fprintf(fp,"%f %f %f ",(float)wallParticleForce2[i*3],(float)wallParticleForce2[i*3+1],(float)wallParticleForce2[i*3+2]);
//	}
//	fprintf(fp,"\n        </DataArray>\n");
	
	if(outputNonNewtonian)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='RHO' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {
				if(particleBC[i] == surface || particleType[i] == wall)
					fprintf(fp,"%f ",(float)RHO[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='ConcVol' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {
				if(particleBC[i] == surface || particleType[i] == wall)
					fprintf(fp,"%f ",(float)Cv[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='MEU' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {
				if(particleBC[i] == surface || particleType[i] == wall)
					fprintf(fp,"%f ",(float)MEU[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='MEUy' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {
				if(particleBC[i] == surface || particleType[i] == wall)
					fprintf(fp,"%f ",(float)MEU_Y[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='II' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {
				if(particleBC[i] == surface || particleType[i] == wall)
					fprintf(fp,"%f ",(float)II[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Int32' Name='PTYPE' format='ascii'>\n");
			for(int i=0; i<numParticles; i++) {
				if(particleBC[i] == surface || particleType[i] == wall)
					fprintf(fp,"%d ",PTYPE[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='p_smooth' format='ascii'>\n");
		for(int i=0; i<numParticles; i++) {
			if(particleBC[i] == surface || particleType[i] == wall)
				fprintf(fp,"%f ",(float)p_smooth[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}
	
	if(outputAuxiliar)
	{

//		fprintf(fp,"        <DataArray type='Int32' Name='ParticleType' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {
//			if(particleBC[i] == surface || particleType[i] == wall)
//				fprintf(fp,"%d ",particleType[i]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

//		fprintf(fp,"        <DataArray type='Float32' Name='mirrorParticlePos' NumberOfComponents='3' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {
//			if(particleBC[i] == surface || particleType[i] == wall)
//			fprintf(fp,"%f %f %f ",(float)mirrorParticlePos[i*3],(float)mirrorParticlePos[i*3+1],(float)mirrorParticlePos[i*3+2]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

//		fprintf(fp,"        <DataArray type='Int32' Name='NearWall' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {
//			if(particleBC[i] == surface || particleType[i] == wall)
//				fprintf(fp,"%d ",particleNearWall[i]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

		//fprintf(fp,"        <DataArray type='Float32' Name='Fwall' NumberOfComponents='3' format='ascii'>\n");
		//for(int i=0; i<numParticles; i++) {
		//	fprintf(fp,"%f %f %f ",(float)Fwall[i*3],(float)Fwall[i*3+1],(float)Fwall[i*3+2]);
		//}
		//fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='Fwall' NumberOfComponents='3' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {
//			fprintf(fp,"%f %f %f ",(float)forceWall[i*3],(float)forceWall[i*3+1],(float)forceWall[i*3+2]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='DIV' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)velDivergence[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='Di' format='ascii'>\n");
//		for(int i=0; i<numParticles; i++) {	fprintf(fp,"%f ",(float)diffusiveTerm[i]);}
//		fprintf(fp,"\n        </DataArray>\n");

		fprintf(fp,"        <DataArray type='Int32' Name='nearMeshType' format='ascii'>\n");
		for(int i=0; i<numParticles; i++)
		{
			if(particleBC[i] == surface || particleType[i] == wall)
				fprintf(fp,"%d ",nearMeshType[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");

		//fprintf(fp,"        <DataArray type='Float32' Name='numNeighborsSurfaceParticles' format='ascii'>\n");
		//for(int i=0; i<numParticles; i++) {fprintf(fp,"%f ",(float)numNeighborsSurfaceParticles[i]);}
		//fprintf(fp,"\n        </DataArray>\n");
	
	}

	fprintf(fp,"      </PointData>\n");

	// Cells
	fprintf(fp,"      <Cells>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='connectivity' format='ascii'>\n");
	for(int i=0, ii = 0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			fprintf(fp,"%d ",ii);
			ii++;
		}
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='offsets' format='ascii'>\n");
	for(int i=0, ii=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			fprintf(fp,"%d ",ii+1);
			ii++;
		}
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='UInt8' Name='types' format='ascii'>\n");
	for(int i=0; i<numParticles; i++)
	{
		if(particleBC[i] == surface || particleType[i] == wall)
		{
			fprintf(fp,"1 ");
		}
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"      </Cells>\n");
	fprintf(fp,"    </Piece>\n");
	fprintf(fp,"  </UnstructuredGrid>\n");
	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
}

// Write header for vtu files
void MpsParticle::writePvd()
{
	char pvd_filename[256];
	char *output_folder_char = new char[vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, vtuOutputFoldername.c_str());
	// output_folder_char now contains a c-string copy of vtuOutputFoldername

	string vtuOutputFoldername_copy;
	vtuOutputFoldername_copy = vtuOutputFoldername;
	size_t found = vtuOutputFoldername_copy.find('/');
	if(found != std::string::npos)
	{
		vtuOutputFoldername_copy.erase(0,found+1);
	}

	int mkdirOK;
	// Creating a directory
#if defined(_WIN32) || defined(WIN32) || defined(__MINGW32__) || defined(__BORLANDC__)
	mkdirOK = _mkdir(output_folder_char);
#else
	mkdirOK = mkdir(output_folder_char, 0777);
#endif

	if(mkdirOK == -1) {
		printf("Unable to create OUTPUT directory or it has already been created!\n");
	}
	else {
		printf("Directory created.\n");
	}

	sprintf(pvd_filename, "%s.pvd",output_folder_char);

	fp = fopen(pvd_filename, "w");
	int nIter = ceil(timeSimulation/timeStep);

	//fprintf(fp,"<VTKFile type=""Collection"" version=""0.1"" byte_order=""LittleEndian"">\n");
	fprintf(fp,"<VTKFile type='Collection' version='0.1' byte_order='LittleEndian'>\n");
	fprintf(fp,"  <Collection>\n");
	int j = 0;
	for(int i=0;i<nIter;i++)
	{
		if(i % iterOutput == 0)
		{
			double timePrint = timeStep*i;
 			fprintf(fp,"    <DataSet timestep='%.6f' group='A' part='0' file='%s/output%05d.vtu'/>\n",timePrint,vtuOutputFoldername_copy.c_str(),j);
 			j++;
		}
	}
	fprintf(fp,"  </Collection>\n");
 	fprintf(fp,"</VTKFile>\n");

	fclose(fp);

	delete[] output_folder_char;
	output_folder_char = NULL;
}
