/**
 * @defgroup   MPSPARTICLESYSTEM Particle System
 *
 * @brief      This file implements Particle System.
 * @details    Properties shared by all the particles.
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef MPS_INCLUDE_PARTICLESYSTEM_H_
#define MPS_INCLUDE_PARTICLESYSTEM_H_

//#define SHOW_FUNCT_NAME_PART	///print the function name from any location inside a C++ function (useful for investigating programs)

#pragma once

/**
 * @brief      Defines the solid boundary wall type
 */
enum boundaryWallType {
	PARTICLE = 0,	///< Layers of wall and dummy particles
	POLYGON = 1 	///< Polygon mesh wall
};

/**
 * @brief      Defines the particle type
 */
enum partType{
	FLUID = 0,		///< Fluid particle
	WALL = 1,		///< Wall particle
	DUMMY_WALL = 2 	///< Dummy (wall) particle
};

/**
 * @brief      Defines the pressure calculation
 */
enum calcPressType {
	EXPLICIT = 0,			///< Explicit MPS
	WEAKLY = 1,				///< Weakly-compressible MPS
	IMPLICIT_PND = 2,		///< Incompressible MPS with pnd source term
	IMPLICIT_PND_DIVU = 3	///< Incompressible MPS with pnd + divergence of velocity source term
};

/**
 * @brief      Defines the Pressure solver using for Incompressible MPS
 */
enum solvPressType{
	CG = 0,			///< Conjugate gradient
	BICGSTAB = 1	///< Biconjugate gradient stabilized
};

/**
 * @brief      Defines the Particle number density (PND) formulation
 */
enum calcPNDType {
	SUM_WIJ = 0,		///< Sum of the weight funtion
	MEAN_SUM_WIJ = 1,	///< Average Sum of the weight funtion
	DIFFUSIVE = 2		///< Diffusive term used only for explicit of weakly-compressible MPS
};

/**
 * @brief      Defines the free-surface detection
 */
enum calcBCType {
	PND_NEIGH = 0,	///< Based on the particle number density and number of neighbors
	PND_NPCD = 1,	///< Based on the particle number density and centroid deviation
	PND_ARC = 2		///< Sum of the arc method
};

/**
 * @brief      Defines the slip condition at wall
 */
enum slipBC {
	NO_SLIP = 0,	///< No-slip
	FREE_SLIP = 1	///< Free-slip
};

/**
 * @brief      Defines the Polygon wall repulsive force
 */
enum repForceType {
	HARADA = 0,			///< Harada et al. 2008. Improvement of Wall boundary calculation model for MPS method.
	MITSUME = 1,		///< Mitsume et al. 2015. Explicitly represented polygon wall boundary model for the explicit MPS method.
	LENNARD_JONES = 2,	///< Monaghan JJ, 1994. Simulating free surface flows with SPH.
	MONAGHAN_KAJTAR = 3	///< Monaghan et al. 2009. SPH particle boundary forces for arbitrary boundaries.
};

/**
 * @brief      Defines the Polygon wall motion
 */
enum meshType{
	FIXED = 0,		///< Fixed wall
	DEFORMABLE = 1,	///< Deformable wall
	FORCED = 2		///< Forced wall
};

/**
 * @brief      Defines the fluid viscosity behavior
 */
enum viscType{
	NEWTONIAN = 0,		///< Newtonian fluid
	NON_NEWTONIAN = 1	///< Non-Newtonian fluid
};

/**
 * @brief      Defines the particle collision method
 */
enum colType {
	PC = 0,		///< Particle Collision
	DPC = 1		///< Dynamic Particle Collision
};

/**
 * @brief      Domain boundary condition
 */
enum domainBC {
	NONE_BC = 0,	///< None
	PERIODIC = 1	///< Periodic
};

/**
 * @brief      This class describes the particle system engine.
 * @details    Properties shared by all the particles.
 */
class MpsParticleSystem {
public:
	// Constructor declaration
	MpsParticleSystem();
	// Destructor declaration
	virtual ~MpsParticleSystem();

	///////////// JSON FILE INPUTS /////////////
	int wallType;	///< Type of solid wall - Polygon or Layers of particles (Wall + Dummy)
	bool femOn;		///< Finite element method - Off:0, On:1
	bool forcedOn;	///< Forced solid flag - Off:0, On:1
	
	// Geometry dimensions limits
	double domainMinX;		///< Minimum value in the x direction of the analysis domain (m)
	double domainMinY;		///< Minimum value in the y direction of the analysis domain (m)
	double domainMinZ;		///< Minimum value in the z direction of the analysis domain (m)
	double domainMaxX;		///< Maximum value in the x direction of the analysis domain (m)
	double domainMaxY;		///< Maximum value in the y direction of the analysis domain (m)
	double domainMaxZ;		///< Maximum value in the z direction of the analysis domain (m)

	// Domain Boundary Condition
	int domainTypeBC;		///< Type of Domain Boundary Condition
	int numBC;				///< Number of boundary conditions
	int limitTypeBC;		///< Type of domain limit to be used in the Boundary Condition
	bool periodicDirectionX;		///< Input periodic direction in domain X
	bool periodicDirectionY;		///< Input periodic direction in domain Y
	bool periodicDirectionZ;		///< Input periodic direction in domain Z

	// Physical parameters
	double densityFluid;	///< Fluid particle density (kg/m3)
	double densityWall;		///< Wall particle density (kg/m3)
	double KNM_VS1;			///< Kinematic viscosity phase 1 (m2/s)
	int fluidType;			///< Newtonian:0 , Non Newtonian:1
	double gravityX;		///< Gravity x (m/s2)
	double gravityY;		///< Gravity y (m/s2)
	double gravityZ;		///< Gravity z (m/s2)
	
	// Rheological parameters
	double KNM_VS2;			///< Kinematic viscosity phase 2 (m2/s)
	double DNS_FL1;			///< Fluid particle density phase 1 (kg/m3)
	double DNS_FL2;			///< Fluid particle density phase 2 (kg/m3)
	double DNS_SDT;			///< Sediment density (kg/m3)
	double N;				///< flow behaviour (power law) index
	double MEU0;			///< consistency index
	double PHI_1;			///< friction angle (RAD) lower limit
	double PHI_WAL;			///< friction angle (RAD)
	double PHI_BED;			///< friction angle (RAD)
	double PHI_2;			///< second friction angle Values are based on  Minatti & Paris (2015) upper limit
	double cohes;			///< cohesiveness coefficient
	int Fraction_method;	///< Method of calculation of volume of fraction. 1: Linear dist across the interface, 2: smoothed value
	double DG;				///< grain size
	double I0;				///< I0 value in Meu9I0 rheology     Values are based on  Minatti &  Paris (2015)
	double mm;				///< Regularization parameter. Control the viscosity grows
	int stress_calc_method;	///< Method 1; viscosity is directly used in momentum equation. Method 2: first the stress tensor is calculated then it is used in momentum equation
	int visc_itr_num;
	double visc_error;
	double visc_ave;
	double Cd;				///< Drag coefficient
	double VF_min;			///< Minimum volume fraction
	double VF_max;			///< Maximum volume fraction

	// Numerical
	double dim;				///< Dimension
	double partDist;		///< Average particle distance (m)
	double timeStep;		///< Time step (s)
	double timeSimulation;	///< Time of simulation (s)
	int iterOutput;			///< Number of iterations to determine the output interval
	double cflNumber;		///< Courant (CFL) condition number
	int weightType;			///< Weight function re/rij - 1 = 0; re/rij + rij/re - 2 = 1; re/rij - rij/re = 2; pow(1-rij/re,3.0) = 3;
	int slipCondition;		///< No-slip = 0 ; Free-slip = 1
	double reS;				///< Influence radius small
	double reL;				///< Influence radius large
	int gradientType;		///< Pressure gradient: Pj - Pmin = 0; Pj + Pi = 1; Pj + Pi - 2*Pmin = 2; ni*Pj/nj + nj*Pi/ni = 3
	bool gradientCorrection;// Corrected pressure gradient. No = 0; Yes = 1
	double relaxPress;		///< Relaxation factor for pressure correction (<= 1)
	int mpsType;			///< Explicit MPS = 0;  Weakly compressible MPS = 1
	double soundSpeed;		///< Sound speed (m/s) (10)
	double gamma;			///< Gamma weakly compressible MPS
	int solverType;					///< Solver type
	double alphaCompressibility;	///< Compressibility leads to a diagonal dominant matrix
	double relaxPND;		///< Relaxation coefficent for pressure source term PND
	double shiftingType;	///< Adjusted velocity: No = 0; DW*Uij = 1; GradCi = 2
	double dri;				///< Adjusted velocity paramater (DW*Uij) dri <= 0.01
	double coefA;			///< Dimensionless number 1-6 (GradCi) coefA = 2.0 provide a good compromise
	double machNumber;		///< Mach number 0.1 (GradCi)
	double VEL_A;			///< Adjusted velocity paramater (3) a = 0.9 (NOT IMPLEMENTED !!!)
	int pndType;			///< PND = Soma(wij) = 0; PND = soma(PNDj)*wij/soma(wij) = 1; PND = Diffusive term = 2
	double diffusiveCoef;	///< Diffusive term coefficient
	int repulsiveForceType;				///< Wall repulsive force: Harada = 0, Mitsume = 1; Lennard-Jones = 2; Monaghan-Kajtar = 3
	double reRepulsiveForce;			///< Influence radius for repulsive force
	double expectMaxVelocity;			///< Expected maximum velocity
	double repForceCoefMitsume;			///< Wall coefficent repulsive force (10000/100000) (500000) (Mitusme)
	double repForceCoefLennardJones;	///< Wall coefficent repulsive force (1-10) (Lennard-Jones)
	double repForceCoefMonaghanKajtar;	///< Wall coefficent repulsive force (1-10) (Monaghan-Kajtar)
	double EPS_RE;			///< 
	int freeSurfType;		///< Free surface condition
	double pndThreshold;	///< Surface threshold PND (2D - 0.97/ 3D - 0.93)
	double neighThreshold;	///< Surface threshold Neighbors
	double npcdThreshold;	///< Surface threshold NPCD
	double thetaThreshold;	///< Surface threshold ARC
	double normThreshold;	///< Surface threshold Normal
	int collisionType;		///< Particle collision type (PC or DPC)
	double collisionRatio;	///< Collision ratio
	double distLimitRatio;	///< Coefficient of distance which does not allow any further access between particles (0.9)
	double lambdaCollision;		///< Non-dimensional coefficient used to adjust background pressure
	int ghost;				///< Ghost particle ID
	int fluid;				///< Fluid particle ID
	int wall;				///< Wall particle ID
	int dummyWall;			///< Dummy wall particle ID
	int surface;			///< Free-surface particle BC
	int inner;				///< Inner particle BC
	int other;				///< Other particle BC
	int numPartTypes;		///< Number of particle types

	///////////// PARTICLE SYSTEM VARIABLES /////////////

	const double epsilonZero = 10.0*std::numeric_limits<double>::epsilon();	///< Very small number near to zero

	int numOfRigidMesh;		///< Number of fixed rigid meshs
	int numOfDeformableMesh;// Number of deformable meshs
	int numOfForcedMesh;	///< Number of forced meshs
	int numOfMeshs;			///< Number of meshs
	
	// Physical domain limits
	double physDomMinX;		///< Minimum value in the x direction of the physical domain (m)
	double physDomMinY;		///< Minimum value in the y direction of the physical domain (m)
	double physDomMinZ;		///< Minimum value in the z direction of the physical domain (m)
	double physDomMaxX;		///< Maximum value in the x direction of the physical domain (m)
	double physDomMaxY;		///< Maximum value in the y direction of the physical domain (m)
	double physDomMaxZ;		///< Maximum value in the z direction of the physical domain (m)

	// Buckets
	int numBucketsX;			///< Number of buckets in the x direction in the analysis domain
	int numBucketsY;			///< Number of buckets in the y direction in the analysis domain
	int numBucketsZ;			///< Number of buckets in the z direction in the analysis domain
	int numBucketsXY;			///< Number of buckets in the x and y directions in the analysis domain
	int numBucketsXYZ;			///< Number of (all) buckets in analysis area
	double bucketSide;			///< Length of one bucket side
	double bucketSide2;			///< Length of one bucket side to square
	double invBucketSide;		///< Inverse of length of one bucket side
	double periodicLength[3];	///< Periodic length in x, y and z direction
	
	// Model
	double NEU;
	double coeffViscosity;		///< Coefficient used to calculate viscosity term Multiphase
	double coeffViscMultiphase;	///< Coefficient used to calculate viscosity term Multiphase
	double coeffPressGrad;		///< Coefficient used to calculate pressure gradient term
	double coeffPressEMPS;		///< Coefficient used to calculate pressure E-MPS
	double coeffPressWCMPS;		///< Coefficient used to calculate pressure WC-MPS
	double coeffShifting1;		///< Coefficient used to adjust velocity type 1
	double coeffShifting2;		///< Coefficient used to adjust velocity type 2
	double coeffPPE;			///< Coefficient used to PPE
	double coeffPPESource;		///< Coefficient used to PPE source term
	double pndGradientZero;		///< Initial particle number density (gradient operator)
	double pndSmallZero;		///< Initial particle number density (small)
	double pndLargeZero;		///< Initial particle number density (large)
	double lambdaZero;			///< Coefficient λ of Laplacian model
	double numNeighZero;		///< Initial number of neighbors
	double distCollisionLimit;	///< A distance that does not allow further access between particles
	double distCollisionLimit2;	///< distCollisionLimit to square
	double restitutionCollision;///< Restitution related to Kinetic energy variation of the particles

	// Numerical
	double eps_reS;
	double betaPnd;			///< Surface cte PND
	double betaNeigh;		///< Surface cte Neighbors
	double delta2;			///< Surface cte NPCD
	double thetaArc;		///< Surface cte theta ARC
	double hThreshold2;		///< Surface cte radius ARC
	double dstThreshold2;	///< Surface cte radius ARC
	double normThreshold2;	///< Surface cte Normal
	int numOfIterations;	///< Number of iterations
	int fileNumber;			///< File number
	double timeCurrent;		///< Current simulation time
	double velMax;			///< Current Maximum flow velocity
	double CFLcurrent;		///< Current Courant number
	double CFLvisc;			///< Current CFL number related to the viscosity
	double timerStart;		///< Initial time
	double timerEnd;		///< Current time
	double invPartDist;		///< Inverse of average particle distance (m)
	double reS2;			///< Influence radius small to square
	double reL2;			///< Influence radius large to square

	// Forced rigid wall motion
	double uniformVelWall[3]; ///< Uniform wall velocity

protected:
	
private:
	

};

#endif // MPS_INCLUDE_PARTICLESYSTEM_H_