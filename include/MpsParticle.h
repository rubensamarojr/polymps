// Copyright (c) 2021 Rubens AMARO
// Distributed under the MIT License.

#pragma once

enum boundaryWallType {
	PARTICLE = 0,
	POLYGON = 1
};

enum calcPressType {
	EXPLICIT = 0,
	WEAKLY = 1
};

enum calcPNDType {
	SUM_WIJ = 0,
	MEAN_SUM_WIJ = 1,
	DIFFUSIVE = 2
};

enum slipBC {
	NO_SLIP = 0,
	FREE_SLIP = 1
};

enum repForceType {
	HARADA = 0,
	MITSUME = 1,
	LENNARD_JONES = 2,
	MONAGHAN_KAJTAR = 3
};

enum meshType{
	FIXED = 0,
	DEFORMABLE = 1,
	FORCED = 2
};

enum viscType{
	NEWTONIAN = 0,
	NON_NEWTONIAN = 1
};


class MpsParticle {
public:

	MpsParticle();
	virtual ~MpsParticle();
	// Return time
	double getTime();
	// Display simulation informations at each nIter iterations
	void displayInfo(const int intervalIter);
	// Read input data from file .json to class MpsParticle
	void readInputFile();
	// Read data from file .grid to class MpsParticle
	void readMpsParticleFile(const std::string& grid_file);
	// Call functions to write output files
	void writeOutputFiles();
	// Write data. Format .prof
	void writeProfAscii();
	template <typename T> void SwapEnd(T& var);
	// Write data. Format .vtu (Paraview)
	void writeVtuBinary();
	// Write data. Format .vtu (Paraview) - Only free-surface particles
	void writeVtuBinaryFreeSurface();
	// Write data. Format .vtu (Paraview)
	void writeVtuAscii();
	// Write data. Format .vtu (Paraview) - Only free-surface particles
	void writeVtuAsciiFreeSurface();
	// Write header for vtu files
	void writePvd();
	// Weight function
	double weight(const double dst, const double re, const int wijType);
	// Weight function for gradient
	double weightGradient(const double dst, const double re, const int wijType);
	// Derivate of weight function
	double delWeight(const double dst, const double re, const int wijType);
	// Verify if particle is out of domain
	void checkParticleOutDomain(const int i);
	// Allocation of buckets
	// Murotani et al., 2015. Performance improvements of differential operators code for MPS method on GPU.
	void allocateBuckets();
	// Set parameters
	void setParameters();
	// Set initial PND and number of neighbors
	void setInitialPndNumberOfNeigh();
	// Update particle ID's in buckets
	void updateBuckets();
	// Acceleration due to Laplacian of viscosity and gravity
	void calcViscosityGravity();
	// Prediction of pressure gradient
	void predictionPressGradient();
	// Prediction of pressure gradient (Polygon wall)
	void predictionWallPressGradient();
	// Update velocity and position
	void updateVelocityPosition1st();
	// Check collisions between particles
	void checkCollisions();
	// Set force on wall to zero
	//void particle::WallZeroForce_omp(int nNodes, int nSolids, solid_fem * &solid);
	void setWallForceZero(const int nNodes, double *nodeforceX, double *nodeforceY, double *nodeforceZ);
	// Free-surface particles. NPCD (Polygon wall)
	void calcWallNPCD();
	// Compute PND
	void calcPnd();
	// Diffusive term of density/PND (pndType = calcPNDType:DIFFUSIVE)
	// An enhanced weakly-compressible MPS method for free-surface flows
	// https://doi.org/10.1016/j.cma.2019.112771
	void calcPndDiffusiveTerm();
	// Diffusive term of density/PND (Polygon wall) - Free-slip (pndType = calcPNDType:DIFFUSIVE)
	void calcWallSlipPndDiffusiveTerm();
	// Diffusive term of density/PND (Polygon wall) - No-slip (pndType = calcPNDType:DIFFUSIVE)
	void calcWallNoSlipPndDiffusiveTerm();
	// Update PND (pndType = calcPNDType:DIFFUSIVE)
	void updatePnd();
	// Mean PND at wall, dummy or free-surface particles (pndType = calcPNDType:DIFFUSIVE)
	void meanPndParticlesWallDummySurface();
	// Mean PND (Polygon wall) (pndType = calcPNDType:MEAN_SUM_WIJ)
	void meanWallPnd();
	// Mean PND (pndType = calcPNDType:MEAN_SUM_WIJ)
	void meanPnd();
	// Compute pressure EMPS (mpsType = calcPressType::EXPLICIT) and type of particle
	void calcPressEMPSandParticleBC();
	// Compute pressure WCMPS (mpsType = calcPressType::WEAKLY) and type of particle
	void calcPressWCMPSandParticleBC();
	// Extrapolate pressure to wall and dummy particles
	void extrapolatePressParticlesWallDummy();
	// Extrapolate pressure to inner particles near polygon walls
	void extrapolatePressParticlesNearPolygonWall();
	// Determinant of matrix
	double detMatrix (double M11, double M12, double M13, double M21, double M22, double M23, double M31, double M32, double M33);
	// Inverse of matrix
	int inverseMatrix (int dim, double &M11, double &M12, double &M13, double &M21, double &M22, double &M23, double &M31, double &M32, double &M33);
	// Correction matrix
	void correctionMatrix();
	// Acceleration due to pressure gradient
	void calcPressGradient();
	// Acceleration due to pressure gradient (Polygon wall)
	void calcWallPressGradient();
	// Calculation of the volume of fraction if phase II in the mixture
	void calcVolumeFraction();
	// Viscosity interaction values for "real" fluid particles
	void calcViscosityInteractionVal();
	// Free-slip condition. Viscosity interaction values
	void calcWallSlipViscosityInteractionVal();
	// No-Slip condition. Viscosity interaction values
	void calcWallNoSlipViscosityInteractionVal();
	// Free-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
	void calcWallSlipViscosity();
	// No-Slip condition. Add acceleration due laplacian of viscosity on wall (Polygon wall)
	void calcWallNoSlipViscosity();
	// 3D triangle to xy plane
	// https://math.stackexchange.com/questions/856666/how-can-i-transform-a-3d-triangle-to-xy-plane
	void transformMatrix(double *V1, double *V2, double *V3, double *RM);
	// 3D -> 2D
	void transform3Dto2D(double *P1, double *RM);
	// 2D -> 3D
	void transform2Dto3D(double *P1, double *RM);
	// Force on wall due fluid particles - FSI
//	void forceParticlesToWall(mesh mesh, solid_fem * &solid);
	// Update velocity and positions
	void updateVelocityPosition2nd();
	// Shifting technique
	// Improvements for accuracy and stability in a weakly-compressible particle method
	// https://www.sciencedirect.com/science/article/pii/S0045793016302250
	void calcShifting();
	// Normal vector on the fluid
	// An accurate and stable multiphase moving particle semi-implicit method based on a corrective matrix for all particle interaction models
	// https://onlinelibrary.wiley.com/doi/full/10.1002/nme.5844
	void calcNormalParticles();
	// Shifting technique (Polygon wall)
	void calcWallShifting();
	// Concentration and Gradient of concentration
	void calcConcAndConcGradient();
	// Concentration and Gradient of concentration (Polygon wall)
	void calcWallConcAndConcGradient();
	// Normal vector on the fluid
	void calcNormalConcentration();
	// Update velocity at wall and dummy particles
	void updateVelocityParticlesWallDummy();
	// Pressure sensors
	void writePressSensors();
	// Force
	//void writeForce();
	void writeHeaderTxtFiles();

	// Files
// 	// MpsParticle data json file
// 	FILE* js;
// 	// MpsParticle data vtu files
// 	FILE* fp;
// 	// Pressure txt file
// 	FILE* pressTxt;
	// Force txt file
	FILE* forceTxtFile;

// 	std::string gridFilename;
 	std::string meshRigidFilename;
 	std::string meshDeformableFilename;
 	std::string meshForcedFilename;
 	std::string vtuOutputFoldername;
 	std::string forceTxtFilename;
 	std::string pressTxtFilename;

// 	bool outputPnd;
// 	bool outputNeigh;
// 	bool outputDeviation;
// 	bool outputConcentration;
// 	bool outputAuxiliar;
	
// 	// Flags of type
 	bool femOn;
 	bool forcedOn;
 	bool freeSurfWall;
 	int wallType;
 	int vtuType;
 	bool txtPress;
 	bool txtForce;

// 	///////////// INPUTS /////////////
// 	// Geometry dimensions limits
// 	double domainMinX;		// Minimum value in the x direction of the analysis domain (m)
// 	double domainMinY;		// Minimum value in the y direction of the analysis domain (m)
// 	double domainMinZ;		// Minimum value in the z direction of the analysis domain (m)
// 	double domainMaxX;		// Minimum value in the x direction of the analysis domain (m)
// 	double domainMaxY;		// Minimum value in the y direction of the analysis domain (m)
// 	double domainMaxZ;		// Minimum value in the z direction of the analysis domain (m)

// 	// Physical parameters
// 	double densityFluid;	// Fluid particle density (kg/m3)
// 	double densityWall;		// Wall particle density (kg/m3)
// 	double gravityX;		// Gravity x (m/s2)
// 	double gravityY;		// Gravity y (m/s2)
// 	double gravityZ;		// Gravity z (m/s2)
// 	// Rheological parameters
// 	double KNM_VS1;			// Kinematic viscosity phase 1 (m2/s)
// 	double KNM_VS2;			// Kinematic viscosity phase 2 (m2/s)
// 	double DNS_FL1;			// Fluid particle density phase 1 (kg/m3)
// 	double DNS_FL2;			// Fluid particle density phase 2 (kg/m3)
// 	double DNS_SDT;			// Sediment density (kg/m3)
 	int fluidType;		// Newtonian:0 , Non Newtonian:1 
// 	double N;				// flow behaviour (power law) index
// 	double MEU0;			// consistency index
// 	double PHI;				// friction angle (RAD) lower limit
// 	double PHI_WAL;			// friction angle (RAD)
// 	double PHI_BED;			// friction angle (RAD)
// 	double PHI_2;			// second friction angle Values are based on  Minatti & Paris (2015) upper limit
// 	double cohes;			// cohesiveness coefficient
// 	int Fraction_method;	// Method of calculation of volume of fraction. 1: Linear dist across the interface, 2: smoothed value
// 	//int visc_max;			// maximum viscosity uses to avoid singularity
// 	double DG;				// grain size
// 	double I0;				// I0 value in Meu9I0 rheology     Values are based on  Minatti &  Paris (2015)
// 	double mm;				// Regularization parameter. Control the viscosity grows
// 	int stress_calc_method;	// Method 1; viscosity is directly used in momentum equation. Method 2: first the stress tensor is calculated then it is used in momentum equation
// 	int visc_itr_num;
// 	double visc_error;
// 	double visc_ave;
// 	double Cd;				// Drag coefficient
// 	double VF_min;			// Minimum volume fraction
// 	double VF_max;			// Maximum volume fraction

 	// Numerical
 	double dim;				// Dimension
 	double partDist;		// Average particle distance (m)
 	double timeStep;		// Time step (s)
 	double timeSimulation;	// Time of simulation (s)
 	int iterOutput;			// Number of iterations to determine the output interval
// 	double cflNumber;		// Courant (CFL) condition number
 	int numOfIterations;	// Number of iterations
 	int fileNumber;			// File number
 	double timeCurrent, velMax, CFLcurrent, timerStart, timerEnd;	// Variables to show in screen
 	double reS;				// Influence radius small
 	double reL;				// Influence radius large
 	double reS2;			// Influence radius small to square
 	double reL2;			// Influence radius large to square
// 	double eps_reS;
 	int mpsType;			// Explicit MPS = 0;  Weakly compressible MPS = 1
 	int weightType;			// Weight function re/rij - 1 = 0; re/rij + rij/re - 2 = 1; re/rij - rij/re = 2; pow(1-rij/re,3.0) = 3;
 	int slipCondition;		// No-slip = 0 ; Free-slip = 1 
// 	int gradientType;		// Pressure gradient: Pj - Pmin = 0; Pj + Pi = 1; Pj + Pi - 2*Pmin = 2; ni*Pj/nj + nj*Pi/ni = 3
 	bool gradientCorrection;// Corrected pressure gradient. No = 0; Yes = 1
 	double relaxPress;		// Relaxation factor for pressure correction (<= 1)
// 	double soundSpeed;		// Sound speed (m/s) (10)
// 	double gamma;			// Gamma weakly compressible MPS
 	double shiftingType;	// Adjusted velocity: No = 0; DW*Uij = 1; GradCi = 2
// 	double dri;				// Adjusted velocity paramater (DW*Uij) dri <= 0.01
// 	double coefA;			// Dimensionless number 1-6 (GradCi) coefA = 2.0 provide a good compromise
// 	double machNumber;		// Mach number 0.1 (GradCi)
// 	double VEL_A;			// Adjusted velocity paramater (3) a = 0.9 (NOT IMPLEMENTED !!!)
 	int pndType;			// PND = Soma(wij) = 0; PND = soma(PNDj)*wij/soma(wij) = 1; PND = Diffusive term = 2
// 	double diffusiveCoef;	// Diffusive term coefficient
// 	int repulsiveForceType;				// Wall repulsive force: Harada = 0, Mitsume = 1; Lennard-Jones = 2; Monaghan-Kajtar = 3
// 	double reRepulsiveForce;			// Influence radius for repulsive force
// 	double repForceCoefMitsume;			// Wall coefficent repulsive force (10000/100000) (500000) (Mitusme)
// 	double repForceCoefLennardJones;	// Wall coefficent repulsive force (1-10) (Lennard-Jones)
// 	double repForceCoefMonaghanKajtar;	// Wall coefficent repulsive force (1-10) (Monaghan-Kajtar)
// 	double EPS_RE;			// 
// 	double pndThreshold;	// Surface threshold PND (2D - 0.97/ 3D - 0.93)
// 	double neighThreshold;	// Surface threshold Neighboors
// 	double npcdThreshold;	// NPCD threshold
// 	double collisionRatio;	// Collision ratio
// 	double distLimitRatio;	// Coefficient of distance which does not allow any further access between particles (0.9)
 	int ghost;				// Ghost particle ID
 	int fluid;				// Fluid particle ID
// 	int wall;				// Wal particle ID
// 	int surface;			// Free-surface particle BC
// 	int inner;				// Inner particle BC
// 	int other;				// Other particle BC
// 	int numPartTypes;		// Number of particle types
 	int numOfMeshs;			// Number of meshs
 	int numOfRigidMesh;		// Number of fixed rigid meshs
 	int numOfDeformableMesh;// Number of deformable meshs
 	int numOfForcedMesh;	// Number of forced meshs

// 	///////////// INPUTS /////////////

// 	// Buckets
// 	double bucketSide;			// Length of one bucket side
// 	double bucketSide2;			// Length of one bucket side to square
// 	double invBucketSide;		// Inverse of length of one bucket side
// 	int numBucketsX;			// Number of buckets in the x direction in the analysis domain
// 	int numBucketsY;			// Number of buckets in the y direction in the analysis domain
// 	int numBucketsZ;			// Number of buckets in the z direction in the analysis domain
// 	int numBucketsXY;
// 	int numBucketsXYZ;			// Number of buckets in analysis area
// 	int *firstParticleInBucket;	// First particle number stored in the bucket
// 	int *lastParticleInBucket;	// Last particle number stored in the bucket
// 	int *nextParticleInSameBucket;	// Next particle number in the same bucket

// 	// Model
// 	double pndSmallZero;		// Initial particle number density (small)
// 	double pndLargeZero;		// Initial particle number density (large)
// 	double pndGradientZero;		// Initial particle number density (gradient operator)
// 	double lambdaZero;			// Coefficient λ of Laplacian model
// 	double numNeighZero;		// Initial number of neighboors
// 	double coeffViscosity;		// Coefficient used to calculate viscosity term Multiphase
// 	double coeffPressEMPS;		// Coefficient used to calculate pressure E-MPS
// 	double coeffPressGrad;		// Coefficient used to calculate pressure gradient term
// 	double coeffPressWCMPS;		// Coefficient used to calculate pressure WC-MPS
// 	double coeffShifting1;		// Coefficient used to adjust velocity type 1
// 	double coeffShifting2;		// Coefficient used to adjust velocity type 2
// 	double distCollisionLimit;	// A distance that does not allow further access between particles
// 	double distCollisionLimit2;	// distCollisionLimit to square
// 	double restitutionCollision;
	
// 	// Scalars
 	int numParticles;		// Number of particles
 	int *particleType;		// Particle type
// 	int *particleBC;		// BC particle type
// 	int *numNeigh;			// Number of neighboors

	double *press;			// Particle pressure
	double *pressAverage;	// Time averaged particle pressure
// 	double *pndi;			// PND
// 	double *pndSmall;		// PND small = sum(wij)
// 	double *npcdDeviation2;	// NPCD deviation modulus
// 	double *concentration;	// Concentration
// 	double *velDivergence;	// Divergence of velocity
// 	double *diffusiveTerm;	// Diffusive term
// 	double *Dns, *invDns;	// Density and its inverse

// 	// Vectors
// 	double *acc;				// Particle acceleration
// 	double *accStar;			// Particle acceleration due gravity and viscosity
 	double *pos;				// Particle position
// 	double *vel;				// Particle velocity
// 	double *npcdDeviation;		// NPCD deviation
// 	double *gradConcentration;	// Gradient of concentration
// 	double *correcMatrixRow1;	// Correction matrix - Row 1
// 	double *correcMatrixRow2;	// Correction matrix - Row 2
// 	double *correcMatrixRow3;	// Correction matrix - Row 3
// 	double *normal;				// Particle normal

// 	// Polygons
// 	// Scalars
	int *nearMeshType;				// Particle near deformable mesh
 	bool *particleNearWall;			// Particle near wall (true or false)
 	int *numNeighWallContribution;	// Number of neighbors due wall

 	double *pndWallContribution;			// PND wall
// 	double *deviationDotPolygonNormal;		// Deviation vector X polygonal wall
// 	double *numNeighborsSurfaceParticles;	// Number of free-surface particle neighbors
 	double *distParticleWall2;				// Squared distance of particle to triangle mesh
// 	// Vectors
 	double *particleAtWallPos;	// Particle at wall coordinate
 	double *mirrorParticlePos;	// Mirrored particle coordinate
// 	double *wallParticleForce1;	// Wall-Particle force
// 	double *wallParticleForce2;	// Wall-Particle force
	double *polygonNormal;		// Polygon normal
	
// //	double *Posk;			// Particle coordinates
// //	double *Velk;			// Particle velocity
// //	double *Acv;			// Part

// 	// Non-Newtonian
// 	// Scalars
// 	int *PTYPE;				// Type of fluid
	
// 	double coeffViscMultiphase, NEU;
// 	double *Cv;				// Concentration
// 	double *II;				// Invariant
// 	double *MEU;			// Dynamic viscosity
// 	double *MEU_Y;			// Dynamic viscosity ??
// 	double *Inertia;		//
// 	double *pnew;			// New pressure
// 	double *p_rheo_new;		//
// 	double *RHO;			// Fluid density
// 	double *p_smooth;		//
// 	double *VF;				//
// 	double *S12;			//
// 	double *S13;			//
// 	double *S23;			//
// 	double *S11;			//
// 	double *S22;			//
// 	double *S33;			//

// 	// FSI
// 	// Scalars
	int *elementID;			// Element ID
// 	// Vectors
	double *forceWall;		// Force on wall

// 	double *nodeX, *nodeY, *nodeZ;		// Nodal positions
// 	double *nodeDX, *nodeDY, *nodeDZ;	// Nodal displacements
// 	double *nodeUX, *nodeUY, *nodeUZ;	// Nodal velocities
// 	double *nodeFx, *nodeFy, *nodeFz;	// Nodal equivalent forces
// 	int *elemNode1id, *elemNode2id, *elemNode3id; // Global node ID

 	// Forced rigid wall
	double velVWall[3]; // Uniform wall velocity 

private:
	// Files
	// MpsParticle data json file
	FILE* js;
	// MpsParticle data vtu files
	FILE* fp;
	// Pressure txt file
	FILE* pressTxtFile;

	std::string gridFilename;
	// std::string meshRigidFilename;
	// std::string meshDeformableFilename;
	// std::string meshForcedFilename;
	// std::string vtuOutputFoldername;
	// std::string forceTxtFilename;
	// std::string pressTxtFilename;

	bool outputPnd;
	bool outputNeigh;
	bool outputDeviation;
	bool outputConcentration;
	bool outputAuxiliar;
	
	// Flags of type
	// bool femOn;
	// bool forcedOn;
	// bool freeSurfWall;
	// int wallType;
	// int vtuType;

	///////////// INPUTS /////////////
	// Geometry dimensions limits
	double domainMinX;		// Minimum value in the x direction of the analysis domain (m)
	double domainMinY;		// Minimum value in the y direction of the analysis domain (m)
	double domainMinZ;		// Minimum value in the z direction of the analysis domain (m)
	double domainMaxX;		// Minimum value in the x direction of the analysis domain (m)
	double domainMaxY;		// Minimum value in the y direction of the analysis domain (m)
	double domainMaxZ;		// Minimum value in the z direction of the analysis domain (m)

	// Physical parameters
	double densityFluid;	// Fluid particle density (kg/m3)
	double densityWall;		// Wall particle density (kg/m3)
	double gravityX;		// Gravity x (m/s2)
	double gravityY;		// Gravity y (m/s2)
	double gravityZ;		// Gravity z (m/s2)
	// Rheological parameters
	double KNM_VS1;			// Kinematic viscosity phase 1 (m2/s)
	double KNM_VS2;			// Kinematic viscosity phase 2 (m2/s)
	double DNS_FL1;			// Fluid particle density phase 1 (kg/m3)
	double DNS_FL2;			// Fluid particle density phase 2 (kg/m3)
	double DNS_SDT;			// Sediment density (kg/m3)
	//int fluidType;		// Newtonian:0 , Non Newtonian:1
	double N;				// flow behaviour (power law) index
	double MEU0;			// consistency index
	double PHI;				// friction angle (RAD) lower limit
	double PHI_WAL;			// friction angle (RAD)
	double PHI_BED;			// friction angle (RAD)
	double PHI_2;			// second friction angle Values are based on  Minatti & Paris (2015) upper limit
	double cohes;			// cohesiveness coefficient
	int Fraction_method;	// Method of calculation of volume of fraction. 1: Linear dist across the interface, 2: smoothed value
	//int visc_max;			// maximum viscosity uses to avoid singularity
	double DG;				// grain size
	double I0;				// I0 value in Meu9I0 rheology     Values are based on  Minatti &  Paris (2015)
	double mm;				// Regularization parameter. Control the viscosity grows
	int stress_calc_method;	// Method 1; viscosity is directly used in momentum equation. Method 2: first the stress tensor is calculated then it is used in momentum equation
	int visc_itr_num;
	double visc_error;
	double visc_ave;
	double Cd;				// Drag coefficient
	double VF_min;			// Minimum volume fraction
	double VF_max;			// Maximum volume fraction

	// Numerical
	// double dim;				// Dimension
	// double partDist;		// Average particle distance (m)
	// double timeStep;		// Time step (s)
	// double timeSimulation;	// Time of simulation (s)
	// int iterOutput;			// Number of iterations to determine the output interval
	double cflNumber;		// Courant (CFL) condition number
	// int numOfIterations;	// Number of iterations
	// int fileNumber;			// File number
	// double timeCurrent, velMax, CFLcurrent, timerStart, timerEnd;
	// double reS;				// Influence radius small
	// double reL;				// Influence radius large
	// double reS2;			// Influence radius small to square
	// double reL2;			// Influence radius large to square
	double eps_reS;
	// int mpsType;			// Explicit MPS = 0;  Weakly compressible MPS = 1
	// int weightType;			// Weight function re/rij - 1 = 0; re/rij + rij/re - 2 = 1; re/rij - rij/re = 2; pow(1-rij/re,3.0) = 3;
	// int slipCondition;		// No-slip = 0 ; Free-slip = 1
	int gradientType;		// Pressure gradient: Pj - Pmin = 0; Pj + Pi = 1; Pj + Pi - 2*Pmin = 2; ni*Pj/nj + nj*Pi/ni = 3
	// bool gradientCorrection;// Corrected pressure gradient. No = 0; Yes = 1
	//double relaxPress;		// Relaxation factor for pressure correction (<= 1)
	double soundSpeed;		// Sound speed (m/s) (10)
	double gamma;			// Gamma weakly compressible MPS
	// double shiftingType;	// Adjusted velocity: No = 0; DW*Uij = 1; GradCi = 2
	double dri;				// Adjusted velocity paramater (DW*Uij) dri <= 0.01
	double coefA;			// Dimensionless number 1-6 (GradCi) coefA = 2.0 provide a good compromise
	double machNumber;		// Mach number 0.1 (GradCi)
	double VEL_A;			// Adjusted velocity paramater (3) a = 0.9 (NOT IMPLEMENTED !!!)
	// int pndType;			// PND = Soma(wij) = 0; PND = soma(PNDj)*wij/soma(wij) = 1; PND = Diffusive term = 2
	double diffusiveCoef;	// Diffusive term coefficient
	int repulsiveForceType;				// Wall repulsive force: Harada = 0, Mitsume = 1; Lennard-Jones = 2; Monaghan-Kajtar = 3
	double reRepulsiveForce;			// Influence radius for repulsive force
	double repForceCoefMitsume;			// Wall coefficent repulsive force (10000/100000) (500000) (Mitusme)
	double repForceCoefLennardJones;	// Wall coefficent repulsive force (1-10) (Lennard-Jones)
	double repForceCoefMonaghanKajtar;	// Wall coefficent repulsive force (1-10) (Monaghan-Kajtar)
	double EPS_RE;			// 
	double pndThreshold;	// Surface threshold PND (2D - 0.97/ 3D - 0.93)
	double neighThreshold;	// Surface threshold Neighboors
	double npcdThreshold;	// NPCD threshold
	double collisionRatio;	// Collision ratio
	double distLimitRatio;	// Coefficient of distance which does not allow any further access between particles (0.9)
	// int ghost;				// Ghost particle ID
	// int fluid;				// Fluid particle ID
	int wall;				// Wal particle ID
	int surface;			// Free-surface particle BC
	int inner;				// Inner particle BC
	int other;				// Other particle BC
	int numPartTypes;		// Number of particle types
	// 	int numPartTypes;		// Number of particle types
 	// int numOfMeshs;			// Number of meshs
 	// int numOfRigidMesh;		// Number of fixed rigid meshs
 	// int numOfDeformableMesh;// Number of deformable meshs
 	// int numOfForcedMesh;	// Number of forced meshs

	///////////// INPUTS /////////////

	// Buckets
	double bucketSide;			// Length of one bucket side
	double bucketSide2;			// Length of one bucket side to square
	double invBucketSide;		// Inverse of length of one bucket side
	int numBucketsX;			// Number of buckets in the x direction in the analysis domain
	int numBucketsY;			// Number of buckets in the y direction in the analysis domain
	int numBucketsZ;			// Number of buckets in the z direction in the analysis domain
	int numBucketsXY;
	int numBucketsXYZ;			// Number of buckets in analysis area
	int *firstParticleInBucket;	// First particle number stored in the bucket
	int *lastParticleInBucket;	// Last particle number stored in the bucket
	int *nextParticleInSameBucket;	// Next particle number in the same bucket

	// Model
	double pndSmallZero;		// Initial particle number density (small)
	double pndLargeZero;		// Initial particle number density (large)
	double pndGradientZero;		// Initial particle number density (gradient operator)
	double lambdaZero;			// Coefficient λ of Laplacian model
	double numNeighZero;		// Initial number of neighboors
	double coeffViscosity;		// Coefficient used to calculate viscosity term Multiphase
	double coeffPressEMPS;		// Coefficient used to calculate pressure E-MPS
	double coeffPressGrad;		// Coefficient used to calculate pressure gradient term
	double coeffPressWCMPS;		// Coefficient used to calculate pressure WC-MPS
	double coeffShifting1;		// Coefficient used to adjust velocity type 1
	double coeffShifting2;		// Coefficient used to adjust velocity type 2
	double distCollisionLimit;	// A distance that does not allow further access between particles
	double distCollisionLimit2;	// distCollisionLimit to square
	double restitutionCollision;
	
	// Scalars
	// int numParticles;		// Number of particles
	// int *particleType;		// Particle type
	int *particleBC;		// BC particle type
	int *numNeigh;			// Number of neighboors

	// double *press;			// Particle pressure
	// double *pressAverage;	// Time averaged particle pressure
	double *pndi;			// PND
	double *pndSmall;		// PND small = sum(wij)
	double *npcdDeviation2;	// NPCD deviation modulus
	double *concentration;	// Concentration
	double *velDivergence;	// Divergence of velocity
	double *diffusiveTerm;	// Diffusive term
	double *Dns, *invDns;	// Density and its inverse

	// Vectors
	double *acc;				// Particle acceleration
	double *accStar;			// Particle acceleration due gravity and viscosity
	// double *pos;				// Particle position
	double *vel;				// Particle velocity
	double *npcdDeviation;		// NPCD deviation
	double *gradConcentration;	// Gradient of concentration
	double *correcMatrixRow1;	// Correction matrix - Row 1
	double *correcMatrixRow2;	// Correction matrix - Row 2
	double *correcMatrixRow3;	// Correction matrix - Row 3
	double *normal;				// Particle normal

	// Polygons
	// Scalars
	// int *nearMeshType;				// Particle near deformable mesh
	// bool *particleNearWall;			// Particle near wall (true or false)
	// int *numNeighWallContribution;	// Number of neighbors due wall

	// double *pndWallContribution;			// PND wall
	double *deviationDotPolygonNormal;		// Deviation vector X polygonal wall
	double *numNeighborsSurfaceParticles;	// Number of free-surface particle neighbors
	// double *distParticleWall2;				// Squared distance of particle to triangle mesh
	// Vectors
	// double *particleAtWallPos;	// Particle at wall coordinate
	// double *mirrorParticlePos;	// Mirrored particle coordinate
	double *wallParticleForce1;	// Wall-Particle force
	double *wallParticleForce2;	// Wall-Particle force
	// double *polygonNormal;		// Polygon normal
	
//	double *Posk;			// Particle coordinates
//	double *Velk;			// Particle velocity
//	double *Acv;			// Part

	// Non-Newtonian
	// Scalars
	int *PTYPE;				// Type of fluid
	
	double coeffViscMultiphase, NEU;
	double *Cv;				// Concentration
	double *II;				// Invariant
	double *MEU;			// Dynamic viscosity
	double *MEU_Y;			// Dynamic viscosity ??
	double *Inertia;		//
	double *pnew;			// New pressure
	double *p_rheo_new;		//
	double *RHO;			// Fluid density
	double *p_smooth;		//
	double *VF;				//
	double *S12;			//
	double *S13;			//
	double *S23;			//
	double *S11;			//
	double *S22;			//
	double *S33;			//

	// FSI
	// Scalars
	// int *elementID;			// Element ID
	// Vectors
	// double *forceWall;		// Force on wall

	double *nodeX, *nodeY, *nodeZ;		// Nodal positions
	double *nodeDX, *nodeDY, *nodeDZ;	// Nodal displacements
	double *nodeUX, *nodeUY, *nodeUZ;	// Nodal velocities
	double *nodeFx, *nodeFy, *nodeFz;	// Nodal equivalent forces
	int *elemNode1id, *elemNode2id, *elemNode3id; // Global node ID

	// Forced rigid wall
	// double velVWall[3]; // Uniform wall velocity 
};