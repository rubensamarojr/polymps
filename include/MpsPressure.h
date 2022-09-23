/**
 * @defgroup   MPSPRESSURE Pressure
 *
 * @brief      This file implements calculations related to the particle Pressure models.
 * @details    Computes particle pressure using Equation of State (EoS) for explicit or weakly-compressible MPS, and solving
 * the Pressure Poisson Equation (PPE) incompressible MPS. Computes divergence of velocity to be used in the source term of PPE. 
 * Set negative pressures to zero. Extrapolates pressure at particle walls. Compute accelerations due to gradient of pressure.
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef MPS_INCLUDE_PRESSURE_H_
#define MPS_INCLUDE_PRESSURE_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "MpsParticleSystem.h"
#include "MpsParticle.h"
#include "MpsBucket.h"

/**
 * @brief      This class describes the mps particle pressure engine.
 * @details    Computes particle pressure using Equation of State (EoS) for explicit or weakly-compressible MPS, and solving
 * the Pressure Poisson Equation (PPE) incompressible MPS. Computes divergence of velocity to be used in the source term of PPE. 
 * Set negative pressures to zero. Extrapolates pressure at particle walls. Compute accelerations due to gradient of pressure.
 */
class MpsPressure {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsPressure();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsPressure();
	
	/**
	 * @brief      Selects the particle pressure approach.
	 * @details    Adjustment of velocity. Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        Eq. (21) in <a href="https://doi.org/10.1016/j.compfluid.2016.07.014" target="_blank">Improvements for accuracy and stability in a weakly-compressible particle method</a>
	 */
	void calcPress(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Allocates memory for pressure and source term Vector
	 *
	 * @param      Particles  The particles data
	 */
	void allocateMemory(MpsParticle *Particles);

	/**
	 * @brief      Calculates the pressure using Explicit MPS approach.
	 *
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void calcPressEMPS(MpsParticleSystem *PSystem, MpsParticle *Particles);
	
	/**
	 * @brief      Calculates the pressure using Weakly-compressible MPS approach.
	 *
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void calcPressWCMPS(MpsParticleSystem *PSystem, MpsParticle *Particles);
	
	/**
	 * @brief      Calculates the pressure using Incompressible MPS with pnd source term.
	 *
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void solvePressurePoissonPnd(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates the pressure using Incompressible MPS with pnd + divergence of velocity source term.
	 *
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void solvePressurePoissonPndDivU(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates the divergence of the velocity.
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcVelDivergence(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Calculates the pressure using explicit MPS approach.
	 * @details    Contribution of the closest polygon wall. Free-slip boundary condition.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcWallSlipVelDivergence(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Calculates the pressure using explicit MPS approach.
	 * @details    Contribution of the closest polygon wall. No-slip boundary condition.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcWallNoSlipVelDivergence(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Set negative pressure to zero.
	 *
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void setZeroOnNegativePressure(MpsParticle *Particles);

	/**
	 * @brief      Extrapolates pressure to wall and dummy particles.
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void extrapolatePressParticlesWallDummy(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Extrapolates pressure to inner particles near polygon walls.
	 *
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void extrapolatePressParticlesNearPolygonWall(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Calculates predicted acceleration due to the gradient of pressure.
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void predictionPressGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates predicted acceleration due to the gradient of pressure.
	 * @details    Contribution of the closest polygon wall. Add repulsive force if normaliwSqrt < reRepulsiveForce.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void predictionWallPressGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates acceleration due to the gradient of pressure.
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcPressGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates acceleration due to the gradient of pressure.
	 * @details    Contribution of the closest polygon wall. Add repulsive force if normaliwSqrt < reRepulsiveForce.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcWallPressGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Calculates a Repulsive force perpendicular to the polygon wall.
	 * @details    It prevents penetrations of the free-surface particles into polygon walls or inner fluid particles 
	 * knuckled edges of polygon. Force proportional to 1∕Δt^2, that is, noticeably sensitive to any change on Δt.
	 * @param      force       The repulsive force vector
	 * @param[in]  normal      The normal vector
	 * @param[in]  normalSqrt  The normal vector magnitude
	 * @param[in]  i           Particle index
	 * @param      Particles   The particles data
	 * @see        Eq. (20) in <a href="https://doi.org/10.1007/s40571-019-00269-6" target="_blank">Parallel analysis system for free-surface flow using MPS method with explicitly represented polygon wall boundary model</a>
	 */
	void repulsiveForceHarada(double *force, const double *normal, const double normalSqrt, const int i, MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Calculates a Repulsive force perpendicular to the polygon wall.
	 * @details    It prevents penetrations of the free-surface particles into polygon walls or inner fluid particles 
	 * knuckled edges of polygon. Adopts a specific coefficient, which should be tuned for each simulation.
	 * @param      force       The repulsive force vector
	 * @param[in]  normal      The normal vector
	 * @param[in]  normalSqrt  The normal vector magnitude
	 * @param[in]  i           Particle index
	 * @param      Particles   The particles data
	 * @see        Eq. (19) in <a href="https://doi.org/10.1007/s40571-019-00269-6" target="_blank">Parallel analysis system for free-surface flow using MPS method with explicitly represented polygon wall boundary model</a>
	 */
	void repulsiveForceMitsume(double *force, const double *normal, const double normalSqrt, const int i, MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Calculates a Repulsive force perpendicular to the polygon wall.
	 * @details    It prevents penetrations of the free-surface particles into polygon walls or inner fluid particles 
	 * knuckled edges of polygon. Since Drep is proportional to the instantaneous flow field, this formulation provides a
	 * dynamic adjustment of the repulsive force, that is, a little or no tuning of the coefficient Crep is required, thereby
	 * significantly improving the numerical stability for a wide range of simulations.
	 * @param      force       The repulsive force vector
	 * @param[in]  normal      The normal vector
	 * @param[in]  normalSqrt  The normal vector magnitude
	 * @param[in]  vMax2       The squared maximum velocity
	 * @param[in]  i           Particle index
	 * @param      Particles   The particles data
	 * @see        Eq. (38) in <a href="https://doi.org/10.1002/fld.5083" target="_blank">Three-dimensional weakly compressible moving particle simulation coupled with geometrically nonlinear shell for hydro-elastic free-surface flows</a>
	 */
	void repulsiveForceLennardJones(double *force, const double *normal, const double normalSqrt, const double vMax2, const int i, MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Calculates a Repulsive force perpendicular to the polygon wall.
	 * @details    It prevents penetrations of the free-surface particles into polygon walls or inner fluid particles 
	 * knuckled edges of polygon.
	 * @param      force       The repulsive force vector
	 * @param[in]  normal      The normal vector
	 * @param[in]  normalSqrt  The normal vector magnitude
	 * @param[in]  vMax2       The squared maximum velocity
	 * @param[in]  i           Particle index
	 * @param      Particles   The particles data
	 * @see        <a href="https://doi.org/10.1016/j.cpc.2009.05.008" target="_blank">SPH particle boundary forces for arbitrary boundaries</a>
	 */
	void repulsiveForceMonaghanKajtar(double *force, const double *normal, const double normalSqrt, const double vMax2, const int i, MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Public function to get value of solver iterator.
	 *
	 * @return     The solver iterator.
	 */
	int getSolverIter() {
		return solverIter;
	}

	/**
	 * @brief      Public function to get the solver error.
	 *
	 * @return     The solver error.
	 */
	double getSolverError() {
		return solverError;
	}

protected:

	/**
	 * @brief		Solve linear system using Conjugate Gradient.
	 * @details		Diagonal, symmetrical and sparse matrix.
	 * @param[in]	p_mat      The sparse matrix of the PPE
	 * @param  		Particles  The particles data
	 */
	void solveConjugateGradient(Eigen::SparseMatrix<double> p_mat, MpsParticle *Particles);

	/**
	 * @brief 		Solve linear system using Biconjugate Gradient Stabilized.
	 * @details 	Diagonal, symmetrical and sparse matrix.
	 * @param[in]	p_mat      The sparse matrix of the PPE
	 * @param  		Particles  The particles data
	 */
	void solveBiConjugateGradientStabilized(Eigen::SparseMatrix<double> p_mat, MpsParticle *Particles);
	
private:

	int solverIter;					///< Solver number of iterations
	double solverError;				///< Solver estimated error
	
};

#endif // MPS_INCLUDE_PRESSURE_H_