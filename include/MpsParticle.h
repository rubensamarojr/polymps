/**
 * @defgroup   MPSPARTICLE Particle Data and Basic Functions
 *
 * @brief      This file implements Particle data and basic functions.
 * @details    Properties shared by all the particles and some basic functions such as: distance between particles,
 * weight functions, domain lengths used for periodic BC and check particle out of domain.
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef EMPS_INCLUDE_MPSPARTICLE_H_
#define EMPS_INCLUDE_MPSPARTICLE_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "MpsParticleSystem.h"

#pragma once

/**
 * @brief      This class describes the particle engine.
 * @details    Properties shared by all the particles and some basic functions such as: distance between particles,
 * weight functions, domain lengths used for periodic BC and check particle out of domain.
 */
class MpsParticle {
public:
	// Constructor declaration
	MpsParticle();
	// Destructor declaration
	virtual ~MpsParticle();
	
	/**
	 * @brief      Computes the square distance between two particles "i" and "j"
	 * @param[in]  j     Neighboring particle index
	 * @param[in]  rxi   Position X of particle "i"
	 * @param[in]  ryi   Position Y of particle "i"
	 * @param[in]  rzi   Position Z of particle "i"
	 * @param      rx    Component X of vector distance between particles "i" and "j"
	 * @param      ry    Component Y of vector distance between particles "i" and "j"
	 * @param      rz    Component Z of vector distance between particles "i" and "j"
	 * @param      rij2  Particle square distance
	 */
	void sqrDistBetweenParticles(const int j, const double rxi, const double ryi, 
		const double rzi, double &rx, double &ry, double &rz, double &rij2);
	
	/**
	 * @brief      Computes the square distance between two particles "i" and "j"
	 * @details    Distance considering Periodic boundary condition.
	 * @param[in]  j     Neighboring particle index
	 * @param[in]  rxi   Position X of particle "i"
	 * @param[in]  ryi   Position Y of particle "i"
	 * @param[in]  rzi   Position Z of particle "i"
	 * @param      rx    Component X of vector distance between particles "i" and "j"
	 * @param      ry    Component Y of vector distance between particles "i" and "j"
	 * @param      rz    Component Z of vector distance between particles "i" and "j"
	 * @param      rij2  Particle square distance
	 * @param[in]  plx   The domain length in X
	 * @param[in]  ply   The domain length in Y
	 * @param[in]  plz   The domain length in Z
	 */
	void sqrDistBetweenParticles(const int j, const double rxi, const double ryi, 
		const double rzi, double &rx, double &ry, double &rz, double &rij2,
		const double plx, const double ply, const double plz);
	
	/**
	 * @brief      Gets the periodic lengths.
	 * @param[in]  jb       Neighboring bucket index
	 * @param      perlx    The domain length in X (zero if no periodic BC is adopted)
	 * @param      perly    The domain length in Y (zero if no periodic BC is adopted)
	 * @param      perlz    The domain length in Z (zero if no periodic BC is adopted)
	 * @param      PSystem  The physical system
	 */
	void getPeriodicLengths(const int jb, double &perlx, double &perly, double &perlz, MpsParticleSystem *PSystem);
	
	/**
	 * @brief      Weight function
	 * @param[in]  dst      The distance between two particles "i" and "j"
	 * @param[in]  re       Effective radius
	 * @param[in]  wijType  The weight function type
	 * @return     The value of the weight function
	 */
	double weight(const double dst, const double re, const int wijType);
	
	/**
	 * @brief      Weight function for gradient
	 * @param[in]  dst      The distance between two particles "i" and "j"
	 * @param[in]  re       Effective radius
	 * @param[in]  wijType  The weight function type
	 * @return     The value of the weight function
	 */
	double weightGradient(const double dst, const double re, const int wijType);
	
	/**
	 * @brief      Derivate of weight function
	 * @param[in]  dst      The distance between two particles "i" and "j"
	 * @param[in]  re       Effective radius
	 * @param[in]  wijType  The weight function type
	 * @return     The value of the derivative of weight function
	 */
	double delWeight(const double dst, const double re, const int wijType);

	/**
	 * @brief      Add gravity to particle acceleration
	 * @param      PSystem    The physical system
	 * @param      Particles  The particles data
	 */
	void addGravity(MpsParticleSystem *PSystem, MpsParticle *Particles);
	
	/**
	 * @brief      Verify if particle is out of domain
	 * @param      PSystem    The physical system
	 */
	void checkParticleOutDomain(MpsParticleSystem *PSystem);
	
	/**
	 * @brief      Sets the parameters.
	 * @param      PSystem    The physical system
	 */
	void setParameters(MpsParticleSystem *PSystem);
	
	/**
	 * @brief      Set force on wall to zero.
	 * @param[in]  nNodes      The nodes of the mesh
	 * @param      nodeforceX  The componenet X of the force node 
	 * @param      nodeforceY  The componenet Y of the force node
	 * @param      nodeforceZ  The componenet Z of the force node
	 * @param      PSystem     The physical system
	 */
	void setWallForceZero(const int nNodes, double *nodeforceX, double *nodeforceY, double *nodeforceZ, MpsParticleSystem *PSystem);
	

	//void particle::WallZeroForce_omp(int nNodes, int nSolids, solid_fem * &solid);
	// Force on wall due fluid particles - FSI
	// void forceParticlesToWall(mesh mesh, solid_fem * &solid);
	
	// Buckets
	int *firstParticleInBucket;		///< First particle number stored in the bucket
	int *lastParticleInBucket;		///< Last particle number stored in the bucket
	int *nextParticleInSameBucket;	///< Next particle number in the same bucket
	int *bucketPeriodicBC;			///< Periodic Boundary Condition of the bucket

	///////////// INPUTS /////////////
	// Scalars
	int numParticlesZero;	///< Number of particles at the initial instant of simulation
	int numParticles;		///< Number of particles during simulation
	
	///////////// INPUTS - ARRAYS /////////////
	// Scalars
	int    *particleType;	///< Particle type
	double *press;			///< Particle pressure
	double *pressAverage;	///< Time averaged particle pressure
	// Vectors
	double *pos;				///< Particle position
	double *vel;				///< Particle velocity
	double *acc;				///< Particle acceleration
	double *accStar;			///< Particle acceleration due gravity and viscosity
	double *npcdDeviation;		///< NPCD deviation
	double *gradConcentration;	///< Gradient of concentration
	double *correcMatrixRow1;	///< Correction matrix - Row 1
	double *correcMatrixRow2;	///< Correction matrix - Row 2
	double *correcMatrixRow3;	///< Correction matrix - Row 3
	double *normal;				///< Particle normal
	double *dvelCollision;		///< Variation of velocity due collision

	// Polygons
	// Scalars
	int    *nearMeshType;				///< Particle near deformable mesh
	bool   *particleNearWall;			///< Particle near wall (true or false)
	int    *numNeighWallContribution;	///< Number of neighbors due wall
	double *pndWallContribution;		///< PND wall
	double *distParticleWall2;			///< Squared distance of particle to triangle mesh
	// Vectors
	double *particleAtWallPos;	///< Particle at wall coordinate
	double *mirrorParticlePos;	///< Mirrored particle coordinate
	double *polygonNormal;		///< Polygon normal

	// Non-Newtonian
	// Scalars
	int    *PTYPE;			///< Type of fluid
	double *Cv;				///< Concentration
	double *II;				///< Invariant
	double *MEU;			///< Dynamic viscosity
	double *MEU_Y;			///< Dynamic viscosity ??
	double *Inertia;		///< Dimensionless Inertia number
	double *pnew;			///< New pressure
	double *p_rheo_new;		///< New rheology pressure
	double *RHO;			///< Fluid density
	double *p_smooth;		///< Smoothed pressure
	double *VF;				///< Volume fraction
	double *S12;			///< Strain rate tensor. Element 12
	double *S13;			///< Strain rate tensor. Element 13
	double *S23;			///< Strain rate tensor. Element 23
	double *S11;			///< Strain rate tensor. Element 11
	double *S22;			///< Strain rate tensor. Element 22
	double *S33;			///< Strain rate tensor. Element 33

	// FSI
	int    *elementID;		///< Scalar Element ID
	double *forceWall;		///< Vector Force on wall
	
	// Linear system used in Incompressible MPS
	Eigen::VectorXd sourceTerm;		///< Solver source term
	Eigen::VectorXd pressurePPE;	///< Solver pressure

	// Model
	int  *bucketTypeBC;			///< Type of Domain Boundary Condition in Bucket
	bool *periodicDirection;	///< Periodic direction in domain (x, y, or z)
	// double *periodicLength;		///< Periodic length in x, y and z direction
	
	// Scalars
	int *particleBC;		///< BC particle type
	int *numNeigh;			///< Number of neighbors
	double *pndi;			///< Particle number density at particle i
	double *pndki;			///< Particle number density at particle i, step k
	double *pndski;			///< Mean neighbor fluid Particle number density, step k
	double *pndSmall;		///< Particle number density small = sum(wij)
	double *npcdDeviation2;	///< NPCD deviation modulus
	double *concentration;	///< Concentration
	double *velDivergence;	///< Divergence of velocity
	double *diffusiveTerm;	///< Diffusive term
	double *Dns, *invDns;	///< Density and its inverse
	double *pndAuxVar1;		///< Auxiliar variable to be used in PND calculations
	double *pndAuxVar2;		///< Auxiliar variable to be used in PND calculations

	// Polygons
	// Scalars
	double *deviationDotPolygonNormal;		///< Deviation vector X polygonal wall
	double *numNeighborsSurfaceParticles;	///< Number of free-surface particle neighbors
	// Vectors
	double *wallParticleForce1;	///< Wall-Particle force
	double *wallParticleForce2;	///< Wall-Particle force

	// FSI
	// Vectors
	double *nodeX, *nodeY, *nodeZ;		///< Nodal positions
	double *nodeDX, *nodeDY, *nodeDZ;	///< Nodal displacements
	double *nodeUX, *nodeUY, *nodeUZ;	///< Nodal velocities
	double *nodeFx, *nodeFy, *nodeFz;	///< Nodal equivalent forces
	int *elemNode1id, *elemNode2id, *elemNode3id; // Global node ID

protected:
	
private:

};

#endif // EMPS_INCLUDE_MPSPARTICLE_H_