/**
 * @defgroup   MPSINPUTOUTPUT Input/Output
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#if defined(_WIN32) || defined(WIN32) || defined(__MINGW32__) || defined(__BORLANDC__)
#include <direct.h>		///< mkdir for Windows
#endif
#include <experimental/filesystem> 	///< numeric_limits
#include <fstream>		///<> std::ofstream std::ifstream
#include <iostream>					///< cout
#include <sys/stat.h>	///< mkdir for Linux
#include <sys/time.h>	///< gettimeofday
#include "json.hpp"		///< json input file
#include "MpsInputOutput.h"

using namespace std;

// Constructor declaration
MpsInputOutput::MpsInputOutput()
{
}
// Destructor declaration
MpsInputOutput::~MpsInputOutput()
{
}

// Return time
double MpsInputOutput::getTime() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 1.0e-6);
}

void MpsInputOutput::displayInfo(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsPressure *PPress, const int intervalIter) {
	if(PSystem->numOfIterations%intervalIter == 0) {
		PSystem->timerEnd = getTime();
		int seconds, hours, minutes;
		seconds = int(PSystem->timerEnd - PSystem->timerStart);
		//seconds = int(timer_end - timer_sta);
		minutes = seconds / 60;
		hours = minutes / 60;
		printf("Iteration: %5dth Time: %lfsec Num. Particles: %d Max Velocity: %lfm/s Courant: %lf", 
			PSystem->numOfIterations, PSystem->timeCurrent, Particles->numParticles, PSystem->velMax, PSystem->CFLcurrent);
		if(PSystem->fluidType == viscType::NON_NEWTONIAN) {
			printf(" CFLvisc: %lf", PSystem->CFLvisc);
		}
		if(PSystem->mpsType == calcPressType::IMPLICIT_PND || PSystem->mpsType == calcPressType::IMPLICIT_PND_DIVU){
			printf(" Solver iterations: %3d Estimated error: %.2e", PPress->getSolverIter(), PPress->getSolverError());
		}
		printf(" RunTime: %02dh%02dm%02dsec\n", int(hours), int(minutes%60), int(seconds%60));
	}
}

// Read input data from file .json to class MpsParticle
void MpsInputOutput::readInputFile(MpsParticleSystem *PSystem, MpsParticle *Particles) {

	// Runtime start
	PSystem->timerStart = getTime();

	char json_folder[] = "input/";
	char json_file_char [1000];
	char json_path_char [1000];
	bool readOK = false;

	printf("   ____       ___                    ____    ____       							\n");
	printf("  /\\  _`\\    /\\_ \\           /'\\_/`\\/\\  _`\\ /\\  _`\\     				\n");
	printf("  \\ \\ \\L\\ \\__\\//\\ \\   __  __/\\      \\ \\ \\L\\ \\ \\,\\L\\_\\   		\n");
	printf("   \\ \\ ,__/ __`\\ \\ \\ /\\ \\/\\ \\ \\ \\__\\ \\ \\ ,__/\\/_\\__ \\   			\n");
	printf("    \\ \\ \\/\\ \\L\\ \\_\\ \\\\ \\ \\_\\ \\ \\ \\_/\\ \\ \\ \\/   /\\ \\L\\ \\ 	\n");
	printf("     \\ \\_\\ \\____/\\____\\/`____ \\ \\_\\\\ \\_\\ \\_\\   \\ `\\____\\ 		\n");
	printf("      \\/_/\\/___/\\/____/`/___/> \\/_/ \\/_/\\/_/    \\/_____/					\n");
	printf("                          /\\___/                        							\n");
	printf("                          \\/__/                         							\n");

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
		scanf("%999s", json_file_char);
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
	try
	{
		json je = json::parse(ifs,
						/* callback */ nullptr,
						/* allow exceptions */ true,
						/* skip_comments */ true);

		// Access the values
		// Types of Simulations
		try {
			const json jflags = je.at("flags");
			
			PSystem->wallType = jflags.value("wall_type", 1);
			PSystem->femOn = 	jflags.value("fem_MESH", false);
			PSystem->forcedOn = jflags.value("forced_MESH", false);
	
			// Types of Outputs
			try {
				const json joutVtu = jflags.at("output_VTU");
				
				vtuType = 				joutVtu.value("type", 0);
				freeSurfWall = 			joutVtu.value("only_freeSurface", false);
				outputPnd = 			joutVtu.value("pnd", false);
				outputNeigh = 			joutVtu.value("neigh", false);
				outputDeviation = 		joutVtu.value("deviation", false);
				outputConcentration = 	joutVtu.value("concentration", false);
				outputAuxiliar = 		joutVtu.value("auxiliar", false);
				outputNonNewtonian =	joutVtu.value("non_newtonian", false);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}
			
			try {
				const json joutTxt = jflags.at("output_TXT");
				txtPress = joutTxt.value("press", false);
				txtForce = joutTxt.value("force", false);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

		}
		catch (json::out_of_range& e) {
			cout << e.what() << '\n';
		}
	
		// Paths of input and output files/folders
		try {
			const json jpNames = je.at("pathNames");

			gridFilename = 			jpNames.value("particle_grid_file", "oops");
			meshRigidFilename = 	jpNames.value("mesh_rigid_file", "oops");
			meshDeformableFilename =jpNames.value("mesh_deformable_file", "oops");
			meshForcedFilename = 	jpNames.value("mesh_forced_file", "oops");
			vtuOutputFoldername = 	jpNames.value("vtu_output_folder", "oops");
			forceTxtFilename = 		jpNames.value("forceTxt_file", "oops");
			pressTxtFilename = 		jpNames.value("pressTxt_file", "oops");
		}
		catch (json::out_of_range& e) {
			cout << e.what() << '\n';
		}

		// Geometry dimension limits
		try {
			const json jdomain = je.at("domain");
			
			try {
				const json jdomainMin = jdomain.at("min");
				PSystem->domainMinX = jdomainMin.value("x", 0.0);
				PSystem->domainMinY = jdomainMin.value("y", 0.0);
				PSystem->domainMinZ = jdomainMin.value("z", 0.0);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}
			
			try {
				const json jdomainMax = jdomain.at("max");
				PSystem->domainMaxX = jdomainMax.value("x", 0.0);
				PSystem->domainMaxY = jdomainMax.value("y", 0.0);
				PSystem->domainMaxZ = jdomainMax.value("z", 0.0);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			// Domain Periodic boundary condition
			try {
				const json jdomainBoundary = jdomain.at("boundary");

				PSystem->domainTypeBC = jdomainBoundary.value("type", 0);
				PSystem->limitTypeBC = 	jdomainBoundary.value("limit", 0);
				PSystem->numBC = 1;

				try {
					const json jdomainBoundaryDir = jdomainBoundary.at("direction");

					PSystem->periodicDirectionX = jdomainBoundaryDir.value("x", false);
					PSystem->periodicDirectionY = jdomainBoundaryDir.value("y", false);
					PSystem->periodicDirectionZ = jdomainBoundaryDir.value("z", false);
				}
				catch (json::out_of_range& e) {
					cout << e.what() << '\n';

					PSystem->periodicDirectionX = false;
					PSystem->periodicDirectionY = false;
					PSystem->periodicDirectionZ = false;
				}
				
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

		}
		catch (json::out_of_range& e) {
			cout << e.what() << '\n';
		}

		// Physical parameters
		try {
			const json jphysical = je.at("physical");

			PSystem->densityFluid = jphysical.value("fluid_density", 1000.0);
			PSystem->densityWall = 	jphysical.value("wall_density", 1000.0);
			PSystem->KNM_VS1 = 		jphysical.value("kinematic_visc", 0.000001);
			PSystem->fluidType = 	jphysical.value("fluid_type", 0);

			try {
				const json jphysicalGrav = jphysical.at("gravity");

				PSystem->gravityX = jphysicalGrav.value("x", 0.0);
				PSystem->gravityY = jphysicalGrav.value("y", 0.0);
				PSystem->gravityZ = jphysicalGrav.value("z", -9.81);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			// Rheological parameters
			try {
				const json jphysicalRheo = jphysical.at("rheological");

				PSystem->KNM_VS2 = 	jphysicalRheo.value("kinematic_visc_phase_2", 0.000001);
				PSystem->DNS_FL1 = 	jphysicalRheo.value("fluid_density_phase_1", 1000.0);
				PSystem->DNS_FL2 = 	jphysicalRheo.value("fluid_density_phase_2", 1540.0);
				PSystem->DNS_SDT = 	jphysicalRheo.value("sediment_density", 1540.0);
				PSystem->N = 		jphysicalRheo.value("power_law_index", 1.2);
				PSystem->MEU0 = 	jphysicalRheo.value("consistency_index", 0.03);

				try {
					const json jphysicalRheoPhi = jphysicalRheo.at("phi");
					
					PSystem->PHI_1 = 	jphysicalRheoPhi.value("lower", 0.541);
					PSystem->PHI_WAL = 	jphysicalRheoPhi.value("wall", 0.541);
					PSystem->PHI_BED = 	jphysicalRheoPhi.value("bed", 0.541);
					PSystem->PHI_2 = 	jphysicalRheoPhi.value("second", 0.6);
				}
				catch (json::out_of_range& e) {
					cout << e.what() << '\n';
				}

				PSystem->cohes = 				jphysicalRheo.value("cohes_coeff", 0.0);
				PSystem->Fraction_method = 		jphysicalRheo.value("fraction_method", 2);
				// visc_max = jphysicalRheo.value("visc_max", 20);
				PSystem->DG = 					jphysicalRheo.value("grain_size", 0.0035);
				PSystem->I0 = 					jphysicalRheo.value("I0", 0.75);
				PSystem->mm = 					jphysicalRheo.value("mm", 100.0);
				PSystem->stress_calc_method = 	jphysicalRheo.value("stress_calc_method", 1);
				
				try {
					const json jphysicalRheoVisc = jphysicalRheo.at("viscosity");
					
					PSystem->visc_itr_num = jphysicalRheoVisc.value("iter_num", 1);
					PSystem->visc_error = 	jphysicalRheoVisc.value("error", 0.0);
					PSystem->visc_ave = 	jphysicalRheoVisc.value("average", 0.0);
				}
				catch (json::out_of_range& e) {
					cout << e.what() << '\n';
				}

				PSystem->Cd = jphysicalRheo.value("drag_coeff", 0.47);

				try {
					const json jphysicalRheoVol = jphysicalRheo.at("volume_fraction");
					
					PSystem->VF_min = jphysicalRheoVol.value("min", 0.25);
					PSystem->VF_max = jphysicalRheoVol.value("max", 0.65);
				}
				catch (json::out_of_range& e) {
					cout << e.what() << '\n';
				}

			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

		}
		catch (json::out_of_range& e) {
			cout << e.what() << '\n';
		}
		
		// Numerical parameters
		try {
			const json jnumerical = je.at("numerical");
			
			PSystem->dim = 				jnumerical.value("dimension", 3.0);
			PSystem->partDist = 		jnumerical.value("particle_dist", 0.01);
			PSystem->timeStep = 		jnumerical.value("time_step", 0.0005);
			PSystem->timeSimulation = 	jnumerical.value("final_time", 1.0);
			PSystem->iterOutput = 		jnumerical.value("iter_output", 80);
			PSystem->cflNumber = 		jnumerical.value("CFL_number", 0.2);
			PSystem->weightType = 		jnumerical.value("weight_type", 0);
			PSystem->slipCondition = 	jnumerical.value("slip_condition", 0);

			try {
				const json jnumericalRe = jnumerical.at("effective_radius");
				
				PSystem->reS = jnumericalRe.value("small", 2.1);
				PSystem->reL = jnumericalRe.value("large", 2.1);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			try {
				const json jnumericalGrad = jnumerical.at("gradient");
				
				PSystem->gradientType = 		jnumericalGrad.value("type", 0);
				PSystem->gradientCorrection = 	jnumericalGrad.value("correction", false);
				PSystem->relaxPress = 			jnumericalGrad.value("relax_fact", 1.0);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			PSystem->mpsType = jnumerical.value("mps_type", 1);

			try {
				const json jnumericalExp = jnumerical.at("explicit_mps");

				try {
					const json jnumericalExpEq = jnumericalExp.at("equation_state");

					PSystem->soundSpeed = 	jnumericalExpEq.value("speed_sound", 15.0);
					PSystem->gamma = 		jnumericalExpEq.value("gamma", 7.0);
				}
				catch (json::out_of_range& e) {
					cout << e.what() << '\n';
				}
				
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			try {
				const json jnumericalImp = jnumerical.at("semi_implicit_mps");

				PSystem->solverType = jnumericalImp.value("solver_type", 0);

				try {
					const json jnumericalImpWeak = jnumericalImp.at("weak_compressibility");

					PSystem->alphaCompressibility = jnumericalImpWeak.value("alpha", 0.000001);
				}
				catch (json::out_of_range& e) {
					cout << e.what() << '\n';
				}

				try {
					const json jnumericalImpSource = jnumericalImp.at("source_term");

					PSystem->relaxPND = jnumericalImpSource.value("relax_pnd", 0.001);
				}
				catch (json::out_of_range& e) {
					cout << e.what() << '\n';
				}
				
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			try {
				const json jnumericalShift = jnumerical.at("particle_shifting");

				PSystem->shiftingType = jnumericalShift.value("type", 2);
				PSystem->dri = 			jnumericalShift.value("DRI", 0.01);
				PSystem->coefA = 		jnumericalShift.value("coef_A", 2.0);
				PSystem->machNumber = 	jnumericalShift.value("mach_number", 0.1);
				PSystem->VEL_A = 		jnumericalShift.value("adj_vel_A", 0.1);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			try {
				const json jnumericalPnd = jnumerical.at("pnd");
				
				PSystem->pndType = 		jnumericalPnd.value("type", 0);
				PSystem->diffusiveCoef =jnumericalPnd.value("diffusive_coeff", 0.35);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			try {
				const json jnumericalWallForce = jnumerical.at("wall_repulsive_force");
				
				PSystem->repulsiveForceType = 	jnumericalWallForce.value("type", 2);
				PSystem->reRepulsiveForce = 	jnumericalWallForce.value("re", 0.5);
				PSystem->expectMaxVelocity = 	jnumericalWallForce.value("maxVel", 6.0);

				try {
					const json jnumericalWallForceCoef = jnumericalWallForce.at("coefficient");

					PSystem->repForceCoefMitsume = 			jnumericalWallForceCoef.value("Mitsume", 40000000.0);
					PSystem->repForceCoefLennardJones = 	jnumericalWallForceCoef.value("Lennard-Jones", 2.0);
					PSystem->repForceCoefMonaghanKajtar = 	jnumericalWallForceCoef.value("Monaghan-Kajtar", 1.0);
				}
				catch (json::out_of_range& e) {
					cout << e.what() << '\n';
				}

				PSystem->EPS_RE = jnumericalWallForce.value("eps_re", 0.01);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			try {
				const json jnumericalFree = jnumerical.at("free_surface_threshold");
				
				PSystem->freeSurfType = 	jnumericalFree.value("type", 0);
				PSystem->pndThreshold = 	jnumericalFree.value("pnd", 0.98);
				PSystem->neighThreshold = 	jnumericalFree.value("neigh", 0.85);
				PSystem->npcdThreshold = 	jnumericalFree.value("NPCD", 0.20);
				PSystem->thetaThreshold = 	jnumericalFree.value("ARC", 45.0);
				PSystem->normThreshold = 	jnumericalFree.value("normal", 0.1);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			try {
				const json jnumericalCol = jnumerical.at("particle_collision");
				
				PSystem->collisionType = 	jnumericalCol.value("type", 0);
				PSystem->collisionRatio = 	jnumericalCol.value("ratio", 0.20);
				PSystem->distLimitRatio = 	jnumericalCol.value("dist_limit_ratio", 0.85);
				PSystem->lambdaCollision = 	jnumericalCol.value("lambda", 0.20);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			try {
				const json jnumericalPtype = jnumerical.at("particle_type");

				PSystem->ghost = 	jnumericalPtype.value("ghost", -1);
				PSystem->fluid = 	jnumericalPtype.value("fluid",  0);
				PSystem->wall = 	jnumericalPtype.value("wall", 2);
				PSystem->dummyWall =jnumericalPtype.value("dummyWall", 3);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

			try {
				const json jnumericalBtype = jnumerical.at("boundary_type");
				
				PSystem->surface = 	jnumericalBtype.value("free_surface", 1);
				PSystem->inner = 	jnumericalBtype.value("inner", 0);
				PSystem->other = 	jnumericalBtype.value("other", -1);
			}
			catch (json::out_of_range& e) {
				cout << e.what() << '\n';
			}

		}
		catch (json::out_of_range& e) {
			cout << e.what() << '\n';
		}
		
	}
	catch (json::parse_error& e)
	{
		// output exception information
		cout 	<< '\n' << "JSON INPUT FILE ERROR !!!" << '\n'
				<< "message: " << e.what() << '\n'
				<< "exception id: " << e.id << '\n'
				<< "byte position of error: " << e.byte << '\n'
				<< "For more information, see https://json.nlohmann.me/home/exceptions" << endl;
		
		throw std::exception();
	}

	PSystem->numPartTypes = 3;

	printf("OK\n");

	printf("Reading GRID File... ");
	readMpsParticleFile(PSystem, Particles);
	printf("OK\n");

	// Error message
	if(PSystem->partDist*PSystem->partDist <= PSystem->epsilonZero) {
		printf("\nError: Operations using squared particle distance lower than the machine double precision."
			" Use a model with a higher particle distance.\n");
		throw std::exception();
	}

	if(PSystem->mpsType == calcPressType::IMPLICIT_PND || PSystem->mpsType == calcPressType::IMPLICIT_PND_DIVU) {
		if(PSystem->pndType == calcPNDType::DIFFUSIVE || PSystem->pndType == calcPNDType::MEAN_SUM_WIJ) {
			printf("\nError: Please, set 'pnd-type to 0' in the json file for incompressible MPS\n");
			throw std::exception();
		}
	}

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
	// cout << "PSystem->densityFluid: " << PSystem->densityFluid << " | ";
	// cout << "densityWall: " << densityWall << " | ";
	// cout << "PSystem->KNM_VS1-2: " << PSystem->KNM_VS1 << ": " << KNM_VS2 << " | ";
	// cout << "DNS_FL1-2-DST: " << DNS_FL1 << ": " << DNS_FL2 << ": " << DNS_SDT << endl;
	// cout << "fluidType: " << fluidType << " | ";
	// cout << "N: " << N << " | ";
	// cout << "MEU0: " << MEU0 << " | ";
	// cout << "PHI-FL-WAL-BED-2: " << PHI_1 << ": " << PHI_WAL << ": " << PHI_BED << ": " << PHI_2 << endl;
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
	// cout << "lo: " << PSystem->partDist << " | ";
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
	// cout << "Particles->numParticles: " << numPartTypes << endl;
		
	// // cout << endl;
	// Close .json file
	fclose(js);
}


// Read data from file .grid to class MpsParticle and allocates memory
void MpsInputOutput::readMpsParticleFile(MpsParticleSystem *PSystem, MpsParticle *Particles) {
	char *grid_file_char = new char[gridFilename.length()+1];
	strcpy (grid_file_char, gridFilename.c_str());
	// grid_file_char now contains a c-string copy of grid_file

	fp = fopen(grid_file_char, "r");
	if(fp == NULL) {
		perror ("Error opening grid file");
		throw std::exception();
	}

	delete[] grid_file_char;
	grid_file_char = NULL;

	int zeroZero;
	fscanf(fp,"%d",&zeroZero);
	fscanf(fp,"%d",&Particles->numParticles);									// Read number of particles
	Particles->numParticlesZero = Particles->numParticles;
	// printf("Number of particles: %d\n",Particles->numParticles);

	// Memory allocation
	// Scalars
	Particles->particleType = (int*)malloc(sizeof(int)*Particles->numParticles);		// Particle type
	Particles->particleBC = (int*)malloc(sizeof(int)*Particles->numParticles);			// BC particle type
	Particles->numNeigh = (int*)malloc(sizeof(int)*Particles->numParticles);			// Number of neighbors

	Particles->press = (double*)malloc(sizeof(double)*Particles->numParticles);			// Particle pressure
	Particles->pressAverage = (double*)malloc(sizeof(double)*Particles->numParticles);	// Time averaged particle pressure
	Particles->pndi = (double*)malloc(sizeof(double)*Particles->numParticles);			// PND
	Particles->pndki = (double*)malloc(sizeof(double)*Particles->numParticles);			// PND step k
	Particles->pndski = (double*)malloc(sizeof(double)*Particles->numParticles);		// Mean fluid neighbor PND step k
	Particles->pndSmall = (double*)malloc(sizeof(double)*Particles->numParticles);		// PND small = sum(wij)
	Particles->npcdDeviation2 = (double*)malloc(sizeof(double)*Particles->numParticles);// NPCD deviation modulus
	Particles->concentration = (double*)malloc(sizeof(double)*Particles->numParticles);	// Concentration
	Particles->velDivergence = (double*)malloc(sizeof(double)*Particles->numParticles);	// Divergence of velocity
	Particles->diffusiveTerm = (double*)malloc(sizeof(double)*Particles->numParticles);	// Diffusive term
	Particles->pndAuxVar1 = (double*)malloc(sizeof(double)*Particles->numParticles);	// Auxiliar PND variable
	Particles->pndAuxVar2 = (double*)malloc(sizeof(double)*Particles->numParticles);	// Auxiliar PND variable

	Particles->Dns = (double*)malloc(sizeof(double)*PSystem->numPartTypes);				// Density
	Particles->invDns = (double*)malloc(sizeof(double)*PSystem->numPartTypes);			// Inverse of Density
	

	// Vectors
	Particles->acc = (double*)malloc(sizeof(double)*Particles->numParticles*3);			// Particle acceleration
	Particles->accStar = (double*)malloc(sizeof(double)*Particles->numParticles*3);		// Particle acceleration due gravity and viscosity
	Particles->pos = (double*)malloc(sizeof(double)*Particles->numParticles*3);			// Particle position
	Particles->vel = (double*)malloc(sizeof(double)*Particles->numParticles*3);			// Particle velocity
	Particles->npcdDeviation = (double*)malloc(sizeof(double)*Particles->numParticles*3);		// NPCD deviation
	Particles->gradConcentration = (double*)malloc(sizeof(double)*Particles->numParticles*3);	// Gradient of concentration
	Particles->correcMatrixRow1 = (double*)malloc(sizeof(double)*Particles->numParticles*3);	// Correction matrix - Row 1
	Particles->correcMatrixRow2 = (double*)malloc(sizeof(double)*Particles->numParticles*3);	// Correction matrix - Row 2
	Particles->correcMatrixRow3 = (double*)malloc(sizeof(double)*Particles->numParticles*3);	// Correction matrix - Row 3
	Particles->normal = (double*)malloc(sizeof(double)*Particles->numParticles*3);				// Particle normal
	Particles->dvelCollision = (double*)malloc(sizeof(double)*Particles->numParticles*3);		// Variation of velocity due collision

	// Polygons
	// Scalars
	Particles->nearMeshType = (int*)malloc(sizeof(int)*Particles->numParticles);			// Type of mesh near particle
	Particles->particleNearWall = (bool*)malloc(sizeof(bool)*Particles->numParticles);		// Particle near polygon wall
	Particles->numNeighWallContribution = (int*)malloc(sizeof(int)*Particles->numParticles);// Number of neighbors due wall

	Particles->pndWallContribution = (double*)malloc(sizeof(double)*Particles->numParticles);		// PND wall
	Particles->deviationDotPolygonNormal = (double*)malloc(sizeof(double)*Particles->numParticles);	// Deviation vector X polygonal wall
	Particles->numNeighborsSurfaceParticles = (double*)malloc(sizeof(double)*Particles->numParticles);// Number of free-surface particle neighbors
	Particles->distParticleWall2 = (double*)malloc(sizeof(double)*Particles->numParticles);			// Squared distance of particle to triangle mesh
	// Vectors
	Particles->particleAtWallPos = (double*)malloc(sizeof(double)*Particles->numParticles*3);	// Particle at wall coordinate
	Particles->mirrorParticlePos = (double*)malloc(sizeof(double)*Particles->numParticles*3);	// Mirrored particle coordinate
	Particles->wallParticleForce1 = (double*)malloc(sizeof(double)*Particles->numParticles*3);	// Wall-Particle force
	Particles->wallParticleForce2 = (double*)malloc(sizeof(double)*Particles->numParticles*3);	// Wall-Particle force
	Particles->polygonNormal = (double*)malloc(sizeof(double)*Particles->numParticles*3);		// Polygon normal
	
//	Posk = (double*)malloc(sizeof(double)*Particles->numParticles*3);		// Particle coordinates
//	Velk = (double*)malloc(sizeof(double)*Particles->numParticles*3);		// Particle velocity
//	Acv = (double*)malloc(sizeof(double)*Particles->numParticles*3);		// Part

	// Non-Newtonian
	// Scalars
	Particles->PTYPE = (int*)malloc(sizeof(int)*Particles->numParticles);				// Type of fluid
	Particles->RHO = (double*)malloc(sizeof(double)*Particles->numParticles);			// Fluid density
	Particles->MEU = (double*)malloc(sizeof(double)*Particles->numParticles);			// Dynamic viscosity

	// Alocate memory only if Non-newtonian is used
	if(PSystem->fluidType == viscType::NON_NEWTONIAN) {
		Particles->Cv = (double*)malloc(sizeof(double)*Particles->numParticles);		// Concentration
		Particles->II = (double*)malloc(sizeof(double)*Particles->numParticles);		// Invariant
		Particles->MEU_Y = (double*)malloc(sizeof(double)*Particles->numParticles);		// Dynamic viscosity ??
		Particles->Inertia = (double*)malloc(sizeof(double)*Particles->numParticles);	//
		Particles->pnew = (double*)malloc(sizeof(double)*Particles->numParticles);		// New pressure
		Particles->p_rheo_new = (double*)malloc(sizeof(double)*Particles->numParticles);//
		Particles->p_smooth = (double*)malloc(sizeof(double)*Particles->numParticles);	//
		Particles->VF = (double*)malloc(sizeof(double)*Particles->numParticles);		//
		Particles->S12 = (double*)malloc(sizeof(double)*Particles->numParticles);		//
		Particles->S13 = (double*)malloc(sizeof(double)*Particles->numParticles);		//
		Particles->S23 = (double*)malloc(sizeof(double)*Particles->numParticles);		//
		Particles->S11 = (double*)malloc(sizeof(double)*Particles->numParticles);		//
		Particles->S22 = (double*)malloc(sizeof(double)*Particles->numParticles);		//
		Particles->S33 = (double*)malloc(sizeof(double)*Particles->numParticles);		//
	}

	// FSI
	// Scalars
	Particles->elementID = (int*)malloc(sizeof(int)*Particles->numParticles);			// Element ID
	// Vectors
	Particles->forceWall = (double*)malloc(sizeof(double)*Particles->numParticles*3);	// Force on wall

	// Solver PPE
	// Particles->pressurePPE = Eigen::VectorXd::Zero(Particles->numParticles);
	// Particles->sourceTerm = Eigen::VectorXd::Zero(Particles->numParticles);

	// Set values from .grid file
	for(int i=0; i<Particles->numParticles; i++) {
		int a[2];
		double b[8];

		// Uncomment here to read .prof file
		//fscanf(fp," %d %d %lf %lf %lf %lf %lf %lf %lf %lf",&a[0],&a[1],&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7]);
		// Uncomment here to read .grid file
		a[0] = 0;
		fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf",&a[1],&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7]);
		Particles->particleType[i]=a[1];
		Particles->pos[i*3]=b[0];	Particles->pos[i*3+1]=b[1];	Particles->pos[i*3+2]=b[2];
		Particles->vel[i*3]=b[3];	Particles->vel[i*3+1]=b[4];	Particles->vel[i*3+2]=b[5];
		//	|				|					|
		//	i*3 = x			i*3+1 = y			i*3+2 = z
		Particles->press[i]=b[6];	Particles->pressAverage[i]=b[7];

//		printf("X: %d %lf %lf %lf %lf %lf %lf %lf %lf\n",particleType[i],pos[3*i],pos[3*i+1],pos[3*i+2],vel[3*i],vel[3*i+1],vel[3*i+2],press[i],pressAverage[i]);
//		Posk[i*3]=b[0];	Posk[i*3+1]=b[1];	Posk[i*3+2]=b[2];
//		Velk[i*3]=b[3];	Velk[i*3+1]=b[4];	Velk[i*3+2]=b[5];
	}
	// Close .grid file
	fclose(fp);
	
	// Set vectors to zero
	for(int i=0;i<Particles->numParticles*3;i++) {
		Particles->acc[i]=0.0;Particles->accStar[i]=0.0;Particles->npcdDeviation[i]=0.0;
		Particles->gradConcentration[i]=0.0;Particles->correcMatrixRow1[i]=0.0;Particles->correcMatrixRow2[i]=0.0;
		Particles->correcMatrixRow3[i]=0.0;Particles->normal[i]=0.0;Particles->dvelCollision[i]=0.0;
		Particles->particleAtWallPos[i]=0.0;Particles->mirrorParticlePos[i]=0.0;Particles->wallParticleForce1[i]=0.0;
		Particles->wallParticleForce2[i]=0.0;Particles->polygonNormal[i]=0.0;Particles->forceWall[i]=0.0;
		//Acv[i]=0.0;
	}

	// Set scalars to zero or infinity(10e8)
	for(int i=0; i<Particles->numParticles; i++) {
		Particles->particleBC[i]=0;Particles->numNeigh[i]=0;Particles->numNeighWallContribution[i]=0;Particles->elementID[i]=0;
		Particles->particleNearWall[i]=false;
		Particles->nearMeshType[i]=meshType::FIXED;

		Particles->pndi[i]=0.0;Particles->pndki[i]=0.0;Particles->pndski[i]=0.0;Particles->pndSmall[i]=0.0;
		Particles->npcdDeviation2[i]=0.0;Particles->concentration[i]=0.0;Particles->velDivergence[i]=0.0;
		Particles->diffusiveTerm[i]=0.0;Particles->pndWallContribution[i]=0.0;Particles->deviationDotPolygonNormal[i]=0.0;
		Particles->numNeighborsSurfaceParticles[i]=0.0;Particles->pndAuxVar1[i]=0.0;Particles->pndAuxVar2[i]=0.0;

		Particles->distParticleWall2[i]=10e8*PSystem->partDist;
	}
	// Set Non-Newtonian scalars to zero
	if(PSystem->fluidType == viscType::NON_NEWTONIAN) {
		for(int i=0; i<Particles->numParticles; i++) {
			Particles->Cv[i]=0.0;Particles->II[i]=0.0;
			Particles->MEU_Y[i]=0.0;Particles->Inertia[i]=0.0;
			Particles->pnew[i]=0.0;Particles->p_rheo_new[i]=0.0;
			Particles->p_smooth[i]=0.0;Particles->VF[i]=0.0;
			Particles->S12[i]=0.0;Particles->S13[i]=0.0;
			Particles->S23[i]=0.0;Particles->S11[i]=0.0;
			Particles->S22[i]=0.0;Particles->S33[i]=0.0;
		}
	}
	// Assign type and density
	for(int i=0; i<Particles->numParticles; i++) {
		/*
		// Assign type and density
		if(pos[i*3+2] <= 0.3) {
			PTYPE[i]=2;
			Particles->RHO[i] = PSystem->DNS_FL2;
			// CHANGED Only at the first time step
			Particles->MEU[i] = PSystem->KNM_VS2 * PSystem->DNS_FL2;
		}
		else {
			PTYPE[i]=1;
			Particles->RHO[i] = PSystem->DNS_FL1;
			// CHANGED Only at the first time step
			Particles->MEU[i] = PSystem->KNM_VS1 * PSystem->DNS_FL1;
		}
		*/

		if(PSystem->fluidType == viscType::NEWTONIAN) {
			Particles->RHO[i] = PSystem->DNS_FL1;
			Particles->PTYPE[i] = 1;
			Particles->MEU[i] = PSystem->KNM_VS1 * PSystem->DNS_FL1;
		}
		// Multiphase simulations - Granular Fluid
		if(PSystem->fluidType == viscType::NON_NEWTONIAN) {
			// Assign type and density
			if(Particles->particleType[i] == 1) {
				Particles->particleType[i] = 0;
				Particles->PTYPE[i] = 2;
				Particles->RHO[i] = PSystem->DNS_FL2;
				// CHANGED Only at the first time step
				Particles->MEU[i] = PSystem->KNM_VS2 * PSystem->DNS_FL2;
			}
			else {
				//Particles->particleType[i] = 0;
				Particles->PTYPE[i] = 1;
				Particles->RHO[i] = PSystem->DNS_FL1;
				// CHANGED Only at the first time step
				Particles->MEU[i] = PSystem->KNM_VS1 * PSystem->DNS_FL1;
			}
		}
	}
}


// Call functions to write output files
void MpsInputOutput::writeOutputFiles(MpsParticleSystem *PSystem, MpsParticle *Particles) {
	// writeProfAscii(Particles);
	// Write particle data (VTU files)
	if(vtuType == 0) {
		if(freeSurfWall == true) {
			writeVtuAsciiFreeSurface(PSystem, Particles);
		}
		else {
			writeVtuAscii(PSystem, Particles);
		}
	}
	else {
		if(freeSurfWall == true) {
			writeVtuBinaryFreeSurface(PSystem, Particles);
		}
		else {
			writeVtuBinary(PSystem, Particles);
		}
	}
}

// Write data. Format .prof
void MpsInputOutput::writeProfAscii(MpsParticleSystem *PSystem, MpsParticle *Particles)
{
	char output_filename[256];
	sprintf(output_filename, "output%05d.prof",PSystem->fileNumber);
	fp = fopen(output_filename, "w");
	fprintf(fp,"%d\n",Particles->numParticles);
	for(int i=0; i<Particles->numParticles; i++)
	{
		int a[2];
		double b[9];
		a[0]=i;	a[1]=Particles->particleType[i];
		b[0]=Particles->pos[i*3];	b[1]=Particles->pos[i*3+1];	b[2]=Particles->pos[i*3+2];
		b[3]=Particles->vel[i*3];	b[4]=Particles->vel[i*3+1];	b[5]=Particles->vel[i*3+2];
		b[6]=Particles->press[i];	b[7]=Particles->pressAverage[i]/PSystem->iterOutput;
		b[8]=Particles->pndi[i];
		fprintf(fp," %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",a[0],a[1],b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[8]);
		Particles->pressAverage[i]=0.0;
	}
	// Close .prof file
	fclose(fp);
}

// https://stackoverflow.com/questions/105252
// https://stackoverflow.com/questions/10913666/error-writing-binary-vtk-files
// https://stackoverflow.com/questions/55829282/write-vtk-file-in-binary-format
template <typename T> void MpsInputOutput::SwapEnd(T& var)
{
	char* varArray = reinterpret_cast<char*>(&var);
	for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
		swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}

// Write data. Format .vtu (Paraview)
void MpsInputOutput::writeVtuBinary(MpsParticleSystem *PSystem, MpsParticle *Particles)
{
	char output_filename[256];
	char *output_folder_char = new char[vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, vtuOutputFoldername.c_str());
	// output_folder_char now contains a c-string copy of vtuOutputFoldername

	//char output_folder_char[vtuOutputFoldername.size() + 1];
    //strcpy(output_folder_char, vtuOutputFoldername.c_str());    // or pass &vtuOutputFoldername[0]

    //string aux_filename = path + "/mesh" + mesh_IDstr + "_" + numstr + ".stl";

	sprintf(output_filename, "%s/output%05d.vtu",output_folder_char,PSystem->fileNumber);

	delete[] output_folder_char;
	output_folder_char = NULL;

	// BINARY FILE
	ofstream file;
	file.open(output_filename, ios::out | ios::binary);

	int nParticles = Particles->numParticles;

	// Header
	//file << "<?xml version='1.0' encoding='UTF-8'?>" << endl;
	file << "<VTKFile type='UnstructuredGrid' version='1.0' byte_order='LittleEndian' header_type='UInt64'>" << endl;
	file << "  <UnstructuredGrid>" << endl;
	file << "    <Piece NumberOfPoints='" << nParticles << "' NumberOfCells='" << nParticles << "'>" << endl;
	
	// Point data
	file << "      <PointData>" << endl;
	file << "        <DataArray type='Float32' Name='Velocity' NumberOfComponents='3' format='binary'>" << endl;
	for(size_t i=0; i<Particles->numParticles; i++)
	{
		{
			float ptx = (float)Particles->vel[i*3  ];
			float pty = (float)Particles->vel[i*3+1];
			float ptz = (float)Particles->vel[i*3+2];

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
	for(size_t i=0; i<Particles->numParticles; i++)
	{
		{
			float ptv = (float)Particles->press[i];
			SwapEnd(ptv);
			file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
		}
	}
	file << endl << "        </DataArray>" << endl;
	
	file << "<        DataArray type='Int32' Name='BC' format='binary'>" << endl;
	for(size_t i=0; i<Particles->numParticles; i++)
	{
		{
			int32_t type = Particles->particleBC[i];
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
		}
	}
	file << endl << "        </DataArray>" << endl;

	if(outputPnd) 
	{
		file << "        <DataArray type='Float32' Name='pnd' format='binary'>" << endl;
		for(size_t i=0; i<Particles->numParticles; i++)
		{
			{
				float ptv = Particles->pndi[i];
				SwapEnd(ptv);
				file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
			}
		}
		file << endl << "        </DataArray>" << endl;
		
		file << "        <DataArray type='Float32' Name='pndSmall' format='binary'>" << endl;
		for(size_t i=0; i<Particles->numParticles; i++)
		{
			{
				float ptv = Particles->pndSmall[i];
				SwapEnd(ptv);
				file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
			}
		}
		file << endl << "        </DataArray>" << endl;
	}

	if(outputNeigh)
	{
		file << "        <DataArray type='Int32' Name='nNeigh' format='binary'>" << endl;
		for(size_t i=0; i<Particles->numParticles; i++)
		{
			{
				int32_t type = Particles->numNeigh[i];
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
	for(size_t i=0; i<Particles->numParticles; i++)
	{
		{
			float ptx = (float)Particles->pos[i*3  ];
			float pty = (float)Particles->pos[i*3+1];
			float ptz = (float)Particles->pos[i*3+2];

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
	for(size_t i=0, ii=0; i<Particles->numParticles; i++)
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
	for(size_t i=0, ii=0; i<Particles->numParticles; i++)
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
	for(size_t i=0; i<Particles->numParticles; i++)
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
void MpsInputOutput::writeVtuBinaryFreeSurface(MpsParticleSystem *PSystem, MpsParticle *Particles)
{

	char output_filename[256];
	char *output_folder_char = new char[vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, vtuOutputFoldername.c_str());
	// output_folder_char now contains a c-string copy of vtuOutputFoldername

	sprintf(output_filename, "%s/output%05d.vtu",output_folder_char,PSystem->fileNumber);

	delete[] output_folder_char;
	output_folder_char = NULL;

	int nParticles = 0;

	for(int i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
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
	for(size_t i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
		{
			float ptx = (float)Particles->vel[i*3  ];
			float pty = (float)Particles->vel[i*3+1];
			float ptz = (float)Particles->vel[i*3+2];

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
	for(size_t i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
		{
			float ptv = (float)Particles->press[i];
			SwapEnd(ptv);
			file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
		}
	}
	file << endl << "        </DataArray>" << endl;
	
	file << "<        DataArray type='Int32' Name='BC' format='binary'>" << endl;
	for(size_t i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
		{
			int32_t type = Particles->particleBC[i];
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
		}
	}
	file << endl << "        </DataArray>" << endl;

	if(outputPnd)
	{
		file << "        <DataArray type='Float32' Name='pnd' format='binary'>" << endl;
		for(size_t i=0; i<Particles->numParticles; i++)
		{
			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
			{
				float ptv = Particles->pndi[i];
				SwapEnd(ptv);
				file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
			}
		}
		file << endl << "        </DataArray>" << endl;
		
		file << "        <DataArray type='Float32' Name='pndSmall' format='binary'>" << endl;
		for(size_t i=0; i<Particles->numParticles; i++)
		{
			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
			{
				float ptv = Particles->pndSmall[i];
				SwapEnd(ptv);
				file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
			}
		}
		file << endl << "        </DataArray>" << endl;
		
	}

	if(outputNeigh)
	{
		file << "        <DataArray type='Int32' Name='nNeigh' format='binary'>" << endl;
		for(size_t i=0; i<Particles->numParticles; i++)
		{
			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
			{
				int32_t type = Particles->numNeigh[i];
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
	for(size_t i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
		{
			float ptx = (float)Particles->pos[i*3  ];
			float pty = (float)Particles->pos[i*3+1];
			float ptz = (float)Particles->pos[i*3+2];

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
	for(size_t i=0, ii=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
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
	for(size_t i=0, ii=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
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
	for(size_t i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
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
void MpsInputOutput::writeVtuAscii(MpsParticleSystem *PSystem, MpsParticle *Particles)
{
	char output_filename[256];
	char *output_folder_char = new char[vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, vtuOutputFoldername.c_str());
	// output_folder_char now contains a c-string copy of vtuOutputFoldername

	sprintf(output_filename, "%s/output%05d.vtu",output_folder_char,PSystem->fileNumber);

	delete[] output_folder_char;
	output_folder_char = NULL;

	int nParticles = Particles->numParticles;

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
	for(int i=0; i<Particles->numParticles; i++)
	{
		fprintf(fp,"%lf %lf %lf ",Particles->pos[i*3],Particles->pos[i*3+1],Particles->pos[i*3+2]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"      </Points>\n");

	// Point data
	fprintf(fp,"      <PointData>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='Velocity' NumberOfComponents='3' format='ascii'>\n");
	for(int i=0; i<Particles->numParticles; i++)
	{
		//double val=sqrt(Particles->vel[i*3]*Particles->vel[i*3]+Particles->vel[i*3+1]*Particles->vel[i*3+1]+Particles->vel[i*3+2]*Particles->vel[i*3+2]);
		fprintf(fp,"%f %f %f ",(float)Particles->vel[i*3],(float)Particles->vel[i*3+1],(float)Particles->vel[i*3+2]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	//fprintf(fp,"        <DataArray type='Float32' Name='preSmallsave' format='ascii'>\n");
	//for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)(Particles->pressAverage[i]/PSystem->iterOutput));}
	//fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='pressure' format='ascii'>\n");
	for(int i=0; i<Particles->numParticles; i++)
	{
		fprintf(fp,"%f ",(float)Particles->press[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='BC' format='ascii'>\n");
	for(int i=0; i<Particles->numParticles; i++)
	{
		fprintf(fp,"%d ",Particles->particleBC[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	
	if(outputPnd)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='pnd' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			fprintf(fp,"%f ",(float)Particles->pndi[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='pndSmall' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			fprintf(fp,"%f ",(float)Particles->pndSmall[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='pndk' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++)
//		{
//			fprintf(fp,"%f ",(float)pndki[i]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='pndsk' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++)
//		{
//			fprintf(fp,"%f ",(float)pndski[i]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");
	}
	if(outputNeigh) {
		fprintf(fp,"        <DataArray type='Int32' Name='nNeigh' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			fprintf(fp,"%d ",Particles->numNeigh[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}

	if(outputDeviation)
	{
//		fprintf(fp,"        <DataArray type='Float32' Name='devSquare' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)npcdDeviation2[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='npcdDeviation' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			fprintf(fp,"%f %f %f ",(float)Particles->npcdDeviation[i*3],(float)Particles->npcdDeviation[i*3+1],(float)Particles->npcdDeviation[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='polygonNormal' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			fprintf(fp,"%f %f %f ",(float)Particles->polygonNormal[i*3],(float)Particles->polygonNormal[i*3+1],(float)Particles->polygonNormal[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='deviationDotPolygonNormal' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)deviationDotPolygonNormal[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
	}

	if(outputConcentration)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='concentration' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			fprintf(fp,"%f ",(float)Particles->concentration[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='GradCi' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			fprintf(fp,"%f %f %f ",(float)Particles->gradConcentration[i*3],(float)Particles->gradConcentration[i*3+1],(float)Particles->gradConcentration[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}
	
//	fprintf(fp,"        <DataArray type='Float32' Name='wallParticleForce1' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<Particles->numParticles; i++) {
//		fprintf(fp,"%f %f %f ",(float)wallParticleForce1[i*3],(float)wallParticleForce1[i*3+1],(float)wallParticleForce1[i*3+2]);
//	}
//	fprintf(fp,"\n        </DataArray>\n");

//	fprintf(fp,"        <DataArray type='Float32' Name='wallParticleForce2' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<Particles->numParticles; i++) {
//		fprintf(fp,"%f %f %f ",(float)wallParticleForce2[i*3],(float)wallParticleForce2[i*3+1],(float)wallParticleForce2[i*3+2]);
//	}
//	fprintf(fp,"\n        </DataArray>\n");

//	fprintf(fp,"        <DataArray type='Float32' Name='Normal' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<Particles->numParticles; i++) {
//		fprintf(fp,"%f %f %f ",(float)normal[i*3],(float)normal[i*3+1],(float)normal[i*3+2]);
//	}
//	fprintf(fp,"\n        </DataArray>\n");
	
	if(outputNonNewtonian)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='RHO' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)Particles->RHO[i]);}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='ConcVol' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)Particles->Cv[i]);}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='MEU' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)Particles->MEU[i]);}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='MEUy' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)Particles->MEU_Y[i]);}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='II' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)Particles->II[i]);}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Int32' Name='PTYPE' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%d ",Particles->PTYPE[i]);}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='p_smooth' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)Particles->p_smooth[i]);}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='Inertia' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)Particles->Inertia[i]);}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='Normal_Stress' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)Particles->p_rheo_new[i]);}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='VF' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)Particles->VF[i]);}
		fprintf(fp,"\n        </DataArray>\n");
	}
	
	if(outputAuxiliar)
	{
		fprintf(fp,"        <DataArray type='Int32' Name='ParticleType' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++) {
			fprintf(fp,"%d ",Particles->particleType[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");

//		fprintf(fp,"        <DataArray type='Float32' Name='mirrorParticlePos' NumberOfComponents='3' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {
//			fprintf(fp,"%f %f %f ",(float)mirrorParticlePos[i*3],(float)mirrorParticlePos[i*3+1],(float)mirrorParticlePos[i*3+2]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

//		fprintf(fp,"        <DataArray type='Int32' Name='NearWall' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {
//			fprintf(fp,"%d ",particleNearWall[i]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

		//fprintf(fp,"        <DataArray type='Float32' Name='Fwall' NumberOfComponents='3' format='ascii'>\n");
		//for(int i=0; i<Particles->numParticles; i++) {
		//	fprintf(fp,"%f %f %f ",(float)Fwall[i*3],(float)Fwall[i*3+1],(float)Fwall[i*3+2]);
		//}
		//fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='Fwall' NumberOfComponents='3' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {
//			fprintf(fp,"%f %f %f ",(float)forceWall[i*3],(float)forceWall[i*3+1],(float)forceWall[i*3+2]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='DIV' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)velDivergence[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='Di' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)diffusiveTerm[i]);}
//		fprintf(fp,"\n        </DataArray>\n");

		if(PSystem->wallType == boundaryWallType::POLYGON)
		{
			fprintf(fp,"        <DataArray type='Int32' Name='nearMeshType' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++)
			{
				fprintf(fp,"%d ",Particles->nearMeshType[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		}

		//fprintf(fp,"        <DataArray type='Float32' Name='numNeighborsSurfaceParticles' format='ascii'>\n");
		//for(int i=0; i<Particles->numParticles; i++) {fprintf(fp,"%f ",(float)numNeighborsSurfaceParticles[i]);}
		//fprintf(fp,"\n        </DataArray>\n");
	}

	fprintf(fp,"      </PointData>\n");

	// Cells
	fprintf(fp,"      <Cells>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='connectivity' format='ascii'>\n");
	for(int i=0, ii = 0; i<Particles->numParticles; i++)
	{
		fprintf(fp,"%d ",ii);
		ii++;
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='offsets' format='ascii'>\n");
	for(int i=0, ii=0; i<Particles->numParticles; i++)
	{
		fprintf(fp,"%d ",ii+1);
		ii++;
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='UInt8' Name='types' format='ascii'>\n");
	for(int i=0; i<Particles->numParticles; i++)
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
void MpsInputOutput::writeVtuAsciiFreeSurface(MpsParticleSystem *PSystem, MpsParticle *Particles)
{
	char output_filename[256];
	char *output_folder_char = new char[vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, vtuOutputFoldername.c_str());
	// output_folder_char now contains a c-string copy of vtuOutputFoldername

	sprintf(output_filename, "%s/output%05d.vtu",output_folder_char,PSystem->fileNumber);

	delete[] output_folder_char;
	output_folder_char = NULL;

	// ASCII FILE
	fp = fopen(output_filename, "w");

	int nParticles = 0;
	for(int i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
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
	for(int i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
		{
			fprintf(fp,"%lf %lf %lf ",Particles->pos[i*3],Particles->pos[i*3+1],Particles->pos[i*3+2]);
		}
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"      </Points>\n");

	// Point data
	fprintf(fp,"      <PointData>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='Velocity' NumberOfComponents='3' format='ascii'>\n");
	for(int i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
			fprintf(fp,"%f %f %f ",(float)Particles->vel[i*3],(float)Particles->vel[i*3+1],(float)Particles->vel[i*3+2]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	//fprintf(fp,"        <DataArray type='Float32' Name='preSmallsave' format='ascii'>\n");
	//for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)(Particles->pressAverage[i]/PSystem->iterOutput));}
	//fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='pressure' format='ascii'>\n");
	for(int i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
			fprintf(fp,"%f ",(float)Particles->press[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='BC' format='ascii'>\n");
	for(int i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
			fprintf(fp,"%d ",Particles->particleBC[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	
	if(outputPnd)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='pnd' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
				fprintf(fp,"%f ",(float)Particles->pndi[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='pndSmall' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
				fprintf(fp,"%f ",(float)Particles->pndSmall[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}

	if(outputNeigh)
	{
		fprintf(fp,"        <DataArray type='Int32' Name='nNeigh' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
				fprintf(fp,"%d ",Particles->numNeigh[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}

	if(outputDeviation)
	{
//		fprintf(fp,"        <DataArray type='Float32' Name='devSquare' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)npcdDeviation2[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='npcdDeviation' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
				fprintf(fp,"%f %f %f ",(float)Particles->npcdDeviation[i*3],(float)Particles->npcdDeviation[i*3+1],(float)Particles->npcdDeviation[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='polygonNormal' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			fprintf(fp,"%f %f %f ",(float)Particles->polygonNormal[i*3],(float)Particles->polygonNormal[i*3+1],(float)Particles->polygonNormal[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='deviationDotPolygonNormal' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)deviationDotPolygonNormal[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
	}

	if(outputConcentration)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='concentration' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
				fprintf(fp,"%f ",(float)Particles->concentration[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='GradCi' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++)
		{
			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
				fprintf(fp,"%f %f %f ",(float)Particles->gradConcentration[i*3],(float)Particles->gradConcentration[i*3+1],(float)Particles->gradConcentration[i*3+2]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}
	
//	fprintf(fp,"        <DataArray type='Float32' Name='wallParticleForce1' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<Particles->numParticles; i++) {
//		fprintf(fp,"%f %f %f ",(float)wallParticleForce1[i*3],(float)wallParticleForce1[i*3+1],(float)wallParticleForce1[i*3+2]);
//	}
//	fprintf(fp,"\n        </DataArray>\n");

//	fprintf(fp,"        <DataArray type='Float32' Name='wallParticleForce2' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<Particles->numParticles; i++) {
//		fprintf(fp,"%f %f %f ",(float)wallParticleForce2[i*3],(float)wallParticleForce2[i*3+1],(float)wallParticleForce2[i*3+2]);
//	}
//	fprintf(fp,"\n        </DataArray>\n");
	
	if(outputNonNewtonian)
	{
		fprintf(fp,"        <DataArray type='Float32' Name='RHO' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {
				if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
					fprintf(fp,"%f ",(float)Particles->RHO[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='ConcVol' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {
				if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
					fprintf(fp,"%f ",(float)Particles->Cv[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='MEU' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {
				if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
					fprintf(fp,"%f ",(float)Particles->MEU[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='MEUy' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {
				if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
					fprintf(fp,"%f ",(float)Particles->MEU_Y[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='II' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {
				if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
					fprintf(fp,"%f ",(float)Particles->II[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Int32' Name='PTYPE' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++) {
				if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
					fprintf(fp,"%d ",Particles->PTYPE[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		fprintf(fp,"        <DataArray type='Float32' Name='p_smooth' format='ascii'>\n");
		for(int i=0; i<Particles->numParticles; i++) {
			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
				fprintf(fp,"%f ",(float)Particles->p_smooth[i]);
		}
		fprintf(fp,"\n        </DataArray>\n");
	}
	
	if(outputAuxiliar)
	{

//		fprintf(fp,"        <DataArray type='Int32' Name='ParticleType' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {
//			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
//				fprintf(fp,"%d ",Particles->particleType[i]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

//		fprintf(fp,"        <DataArray type='Float32' Name='mirrorParticlePos' NumberOfComponents='3' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {
//			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
//			fprintf(fp,"%f %f %f ",(float)mirrorParticlePos[i*3],(float)mirrorParticlePos[i*3+1],(float)mirrorParticlePos[i*3+2]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

//		fprintf(fp,"        <DataArray type='Int32' Name='NearWall' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {
//			if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
//				fprintf(fp,"%d ",particleNearWall[i]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

		//fprintf(fp,"        <DataArray type='Float32' Name='Fwall' NumberOfComponents='3' format='ascii'>\n");
		//for(int i=0; i<Particles->numParticles; i++) {
		//	fprintf(fp,"%f %f %f ",(float)Fwall[i*3],(float)Fwall[i*3+1],(float)Fwall[i*3+2]);
		//}
		//fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='Fwall' NumberOfComponents='3' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {
//			fprintf(fp,"%f %f %f ",(float)forceWall[i*3],(float)forceWall[i*3+1],(float)forceWall[i*3+2]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='DIV' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)velDivergence[i]);}
//		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='Di' format='ascii'>\n");
//		for(int i=0; i<Particles->numParticles; i++) {	fprintf(fp,"%f ",(float)diffusiveTerm[i]);}
//		fprintf(fp,"\n        </DataArray>\n");

		if(PSystem->wallType == boundaryWallType::POLYGON)
		{
			fprintf(fp,"        <DataArray type='Int32' Name='nearMeshType' format='ascii'>\n");
			for(int i=0; i<Particles->numParticles; i++)
			{
				if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
					fprintf(fp,"%d ",Particles->nearMeshType[i]);
			}
			fprintf(fp,"\n        </DataArray>\n");
		}

		//fprintf(fp,"        <DataArray type='Float32' Name='numNeighborsSurfaceParticles' format='ascii'>\n");
		//for(int i=0; i<Particles->numParticles; i++) {fprintf(fp,"%f ",(float)numNeighborsSurfaceParticles[i]);}
		//fprintf(fp,"\n        </DataArray>\n");
	
	}

	fprintf(fp,"      </PointData>\n");

	// Cells
	fprintf(fp,"      <Cells>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='connectivity' format='ascii'>\n");
	for(int i=0, ii = 0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
		{
			fprintf(fp,"%d ",ii);
			ii++;
		}
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='offsets' format='ascii'>\n");
	for(int i=0, ii=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
		{
			fprintf(fp,"%d ",ii+1);
			ii++;
		}
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='UInt8' Name='types' format='ascii'>\n");
	for(int i=0; i<Particles->numParticles; i++)
	{
		if(Particles->particleBC[i] == PSystem->surface || Particles->particleType[i] == PSystem->wall)
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
void MpsInputOutput::writePvd(MpsParticleSystem *PSystem, MpsParticle *Particles)
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
	int nIter = ceil(PSystem->timeSimulation/PSystem->timeStep);

	//fprintf(fp,"<VTKFile type=""Collection"" version=""0.1"" byte_order=""LittleEndian"">\n");
	fprintf(fp,"<VTKFile type='Collection' version='0.1' byte_order='LittleEndian'>\n");
	fprintf(fp,"  <Collection>\n");
	int j = 0;
	for(int i=0;i<nIter;i++)
	{
		if(i % PSystem->iterOutput == 0)
		{
			double timePrint = PSystem->timeStep*i;
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


// Pressure sensors
void MpsInputOutput::writePressSensors(MpsParticleSystem *PSystem, MpsParticle *Particles) {

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
	// xPrs1Min = 1.61 - 2.0*PSystem->partDist; xPrs1Max = 1.61 + PSystem->partDist;
	// yPrs1Min = 0.075 - 2.0*PSystem->partDist; yPrs1Max = 0.075 + 2.0*PSystem->partDist;
	// zPrs1min = 0.003 - PSystem->partDist; zPrs2min = 0.015 - PSystem->partDist; zPrs3min = 0.030 - PSystem->partDist; zPrs4min = 0.080 - PSystem->partDist;
	// zPrs1max = 0.003 + PSystem->partDist; zPrs2max = 0.015 + PSystem->partDist; zPrs3max = 0.030 + PSystem->partDist; zPrs4max = 0.080 + PSystem->partDist;

	// Hydrostatic
	xPrs1Min = 0.1 - 2.0*PSystem->partDist; xPrs1Max = 0.1 + PSystem->partDist; xPrs2Min = 0.0 - 2.0*PSystem->partDist; xPrs2Max = 0.0 + PSystem->partDist;
	yPrs1Min = 0.1 - 2.0*PSystem->partDist; yPrs1Max = 0.1 + 2.0*PSystem->partDist; yPrs2Min = 0.1 - 2.0*PSystem->partDist; yPrs2Max = 0.1 + 2.0*PSystem->partDist;
	zPrs1min = 0.0 - PSystem->partDist; zPrs2min = 0.1 - PSystem->partDist;
	zPrs1max = 0.0 + PSystem->partDist; zPrs2max = 0.1 + PSystem->partDist;

	// Sensor positions
	// Dam 1610
	// posP1[0] = 1.61; posP1[1] = 0.075; posP1[2] = 0.003;
	// posP2[0] = 1.61; posP2[1] = 0.075; posP2[2] = 0.015;
	// posP3[0] = 1.61; posP3[1] = 0.075; posP3[2] = 0.030;
	// posP4[0] = 1.61; posP4[1] = 0.075; posP4[2] = 0.080;

	// Hydrostatic
	posP1[0] = 0.1; posP1[1] = 0.1; posP1[2] = 0.0;
	posP2[0] = 0.0; posP2[1] = 0.1; posP2[2] = 0.1;

#pragma omp parallel for reduction(+: pndP1, pndP2, pndP3, pndP4, P1wij, P2wij, P3wij, P4wij)
	for(int i=0; i<Particles->numParticles; i++) {
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		
		// // Pressure at a specific particle close to the sensor
		// // Dam 1610
		// if(posXi >= xPrs1Min && posXi <= xPrs1Max && posYi >= yPrs1Min && posYi <= yPrs1Max) {
		// 	// Sensor 1
		// 	if(posZi >= zPrs1min && posZi <= zPrs1max) {
		// 		if(Particles->press[i] > 0.0) {
		// 			double v0 = posP1[0] - posXi;
		// 			double v1 = posP1[1] - posYi;
		// 			double v2 = posP1[2] - posZi;
		// 			double dst2 = v0*v0+v1*v1+v2*v2;
		// 			// Closest fluid particle
		// 			if(dst2 < riP1) {
		// 				P1 = Particles->press[i];
		// 				riP1 = dst2;
		// 			}
		// 		}
		// 	}
		// 	// Sensor 2
		// 	if(posZi >= zPrs2min && posZi <= zPrs2max) {
		// 		if(Particles->press[i] > 0.0) {
		// 			double v0 = posP2[0] - posXi;
		// 			double v1 = posP2[1] - posYi;
		// 			double v2 = posP2[2] - posZi;
		// 			double dst2 = v0*v0+v1*v1+v2*v2;
		// 			// Closest fluid particle
		// 			if(dst2 < riP2) {
		// 				P2 = Particles->press[i];
		// 				riP2 = dst2;
		// 			}
		// 		}
		// 	}
		// 	// Sensor 3
		// 	if(posZi >= zPrs3min && posZi <= zPrs3max) {
		// 		if(Particles->press[i] > 0.0) {
		// 			double v0 = posP3[0] - posXi;
		// 			double v1 = posP3[1] - posYi;
		// 			double v2 = posP3[2] - posZi;
		// 			double dst2 = v0*v0+v1*v1+v2*v2;
		// 			// Closest fluid particle
		// 			if(dst2 < riP3) {
		// 				P3 = Particles->press[i];
		// 				riP3 = dst2;
		// 			}
		// 		}
		// 	}
		// 	// Sensor 4
		// 	if(posZi >= zPrs4min && posZi <= zPrs4max) {
		// 		if(Particles->press[i] > 0.0) {
		// 			double v0 = posP4[0] - posXi;
		// 			double v1 = posP4[1] - posYi;
		// 			double v2 = posP4[2] - posZi;
		// 			double dst2 = v0*v0+v1*v1+v2*v2;
		// 			// Closest fluid particle
		// 			if(dst2 < riP4) {
		// 				P4 = Particles->press[i];
		// 				riP4 = dst2;
		// 			}
		// 		}
		// 	}
		// }

		// // Weighted average pressure
		// double v0 = posP1[0] - posXi;
		// double v1 = posP1[1] - posYi;
		// double v2 = posP1[2] - posZi;
		// double dst2 = v0*v0+v1*v1+v2*v2;
		// if(dst2 < PSystem->reS2) {
		// 	double dst = sqrt(dst2);
		// 	double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
		// 	if(Particles->press[i] > 0.0){
		// 		pndP1 += wS;
		// 		P1wij += Particles->press[i]*wS;
		// 	}
		// }
		
		// v0 = posP2[0] - posXi;
		// v1 = posP2[1] - posYi;
		// v2 = posP2[2] - posZi;
		// dst2 = v0*v0+v1*v1+v2*v2;
		// if(dst2 < PSystem->reS2) {
		// 	double dst = sqrt(dst2);
		// 	double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
		// 	if(Particles->press[i] > 0.0){
		// 		pndP2 += wS;
		// 		P2wij += Particles->press[i]*wS;
		// 	}
		// }
		
		// v0 = posP3[0] - posXi;
		// v1 = posP3[1] - posYi;
		// v2 = posP3[2] - posZi;
		// dst2 = v0*v0+v1*v1+v2*v2;
		// if(dst2 < PSystem->reS2) {
		// 	double dst = sqrt(dst2);
		// 	double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
		// 	if(Particles->press[i] > 0.0){
		// 		pndP3 += wS;
		// 		P3wij += Particles->press[i]*wS;
		// 	}
		// }
		
		// v0 = posP4[0] - posXi;
		// v1 = posP4[1] - posYi;
		// v2 = posP4[2] - posZi;
		// dst2 = v0*v0+v1*v1+v2*v2;
		// if(dst2 < PSystem->reS2) {
		// 	double dst = sqrt(dst2);
		// 	double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
		// 	if(Particles->press[i] > 0.0){
		// 		pndP4 += wS;
		// 		P4wij += Particles->press[i]*wS;
		// 	}
		// }
		
		// Pressure at a specific particle close to the sensor
		// Hydrostatic
		if((posXi >= xPrs1Min && posXi <= xPrs1Max && posYi >= yPrs1Min && posYi <= yPrs1Max) || (posXi >= xPrs2Min && posXi <= xPrs2Max && posYi >= yPrs2Min && posYi <= yPrs2Max)) {
			// Sensor 1
			if(posZi >= zPrs1min && posZi <= zPrs1max) {
				if(Particles->press[i] > 0.0) {
					double v0 = posP1[0] - posXi;
					double v1 = posP1[1] - posYi;
					double v2 = posP1[2] - posZi;
					double dst2 = v0*v0+v1*v1+v2*v2;
					// Closest fluid particle
					if(dst2 < riP1) {
						P1 = Particles->press[i];
						riP1 = dst2;
					}
				}
			}
			// Sensor 2
			if(posZi >= zPrs2min && posZi <= zPrs2max) {
				if(Particles->press[i] > 0.0) {
					double v0 = posP2[0] - posXi;
					double v1 = posP2[1] - posYi;
					double v2 = posP2[2] - posZi;
					double dst2 = v0*v0+v1*v1+v2*v2;
					// Closest fluid particle
					if(dst2 < riP2) {
						P2 = Particles->press[i];
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
		if(dst2 < PSystem->reS2) {
			double dst = sqrt(dst2);
			double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
			if(Particles->press[i] > 0.0){
				pndP1 += wS;
				P1wij += Particles->press[i]*wS;
			}
		}
		
		v0 = posP2[0] - posXi;
		v1 = posP2[1] - posYi;
		v2 = posP2[2] - posZi;
		dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2 < PSystem->reS2) {
			double dst = sqrt(dst2);
			double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
			if(Particles->press[i] > 0.0){
				pndP2 += wS;
				P2wij += Particles->press[i]*wS;
			}
		}

	}

	char *pressTxtFilenameChar = new char[pressTxtFilename.length()+1];
	strcpy (pressTxtFilenameChar, pressTxtFilename.c_str());
	// pressTxtFilenameChar now contains a c-string copy of pressTxtFilename

	// Open the File to write
	pressTxtFile = fopen(pressTxtFilenameChar, "a");
	if(pressTxtFile == NULL) {
		perror ("Error opening press txt file");
		throw std::exception();
	}

	delete[] pressTxtFilenameChar;
	pressTxtFilenameChar = NULL;
	
	// // dam 1610
	// fprintf(pressTxtFile,"\n%lf\t%lf\t%lf\t%lf\t%lf",PSystem->timeCurrent, P1,P2,P3,P4);
	// double Pmean;
	
	// if(pndP1 > 0.0)
	// 	Pmean = P1wij/pndP1;
	// else
	// 	Pmean = 0.0;
	// fprintf(pressTxtFile,"\t%lf",Pmean);
	
	// if(pndP2 > 0.0)
	// 	Pmean = P2wij/pndP2;
	// else
	// 	Pmean = 0.0;
	// fprintf(pressTxtFile,"\t%lf",Pmean);
	
	// if(pndP3 > 0.0)
	// 	Pmean = P3wij/pndP3;
	// else
	// 	Pmean = 0.0;
	// fprintf(pressTxtFile,"\t%lf",Pmean);
	
	// if(pndP4 > 0.0)
	// 	Pmean = P4wij/pndP4;
	// else
	// 	Pmean = 0.0;
	// fprintf(pressTxtFile,"\t%lf",Pmean);
	

	// Hydrostatic
	fprintf(pressTxtFile,"\n%lf\t%lf\t%lf",PSystem->timeCurrent, P1,P2);
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


#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

void MpsInputOutput::writeHeaderTxtFiles(MpsParticleSystem *PSystem, MpsParticle *Particles) {

	if(txtForce == true) {
		char *forceTxtFilenameChar = new char[forceTxtFilename.length()+1];
		strcpy (forceTxtFilenameChar, forceTxtFilename.c_str());
		// forceTxtFilenameChar now contains a c-string copy of forceTxtFilename
		
		// Open the File to write
		forceTxtFile = fopen(forceTxtFilenameChar, "w");
		if(forceTxtFile == NULL) {
			perror ("Error opening force txt file");
			throw std::exception();
		}
		fprintf(forceTxtFile,"Time(s)\tForce(N) X\tY\tZ");
		// Close force file
		fclose(forceTxtFile);

		delete[] forceTxtFilenameChar;
		forceTxtFilenameChar = NULL;
	}

	if(txtPress == true) {
		char *pressTxtFilenameChar = new char[pressTxtFilename.length()+1];
		strcpy (pressTxtFilenameChar, pressTxtFilename.c_str());
		// pressTxtFilenameChar now contains a c-string copy of pressTxtFilename

		// Open the File to write
		pressTxtFile = fopen(pressTxtFilenameChar, "w");
		if(pressTxtFile == NULL) {
			perror ("Error opening press txt file");
			throw std::exception();
		}

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

// Write vtk file with initial bucktes
void MpsInputOutput::writeBuckets(MpsParticleSystem *PSystem, MpsParticle *Particles) {

	char output_filename[256];
	char *output_folder_char = new char[vtuOutputFoldername.length()+1];
	strcpy (output_folder_char, vtuOutputFoldername.c_str());
	// output_folder_char now contains a c-string copy of vtuOutputFoldername

	sprintf(output_filename, "%s/buckets.vtk",output_folder_char);

	delete[] output_folder_char;
	output_folder_char = NULL;

	// ASCII FILE
	fp = fopen(output_filename, "w");

	int nDimX = PSystem->numBucketsX + 1;
	int nDimY = PSystem->numBucketsY + 1;
	int nDimZ = PSystem->numBucketsZ;
	if(PSystem->dim==3) {
		nDimZ += 1;
	}
	
	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"Initial Buckets\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET STRUCTURED_POINTS\n");
	fprintf(fp,"DIMENSIONS %d %d %d\n", nDimX, nDimY, nDimZ);
	fprintf(fp,"ORIGIN %lf %lf %lf\n", PSystem->domainMinX, PSystem->domainMinY, PSystem->domainMinZ);
	fprintf(fp,"SPACING %lf %lf %lf\n", PSystem->bucketSide, PSystem->bucketSide, PSystem->bucketSide);
	fprintf(fp,"CELL_DATA %d\n", PSystem->numBucketsXYZ);
	fprintf(fp,"SCALARS density int 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");
	for(int i=0; i<PSystem->numBucketsXYZ; i++) {
		fprintf(fp,"1 ");
	}

	fclose(fp);
}

//
void MpsInputOutput::stringToChar(char *out_char) {
	out_char = new char[vtuOutputFoldername.length()+1];
	strcpy(out_char, vtuOutputFoldername.c_str());
}

// Delete all files inside the simulation folder
void MpsInputOutput::deleteDirectoryFiles() {
	for (const auto& entry : std::experimental::filesystem::directory_iterator(vtuOutputFoldername)) 
		std::experimental::filesystem::remove_all(entry.path());
}