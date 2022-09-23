/**
 * @defgroup   MPSBOUNDARYCONDITION Boundary Conditions
 *
 * @brief      This file implements the Boundary Conditions.
 * @details    Updates type of particle (free-surface, inner), set Bucket as Periodic or None and 
 * manages Inflow/Outflow.
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef MPS_INCLUDE_BOUNDARYCONDITION_H_
#define MPS_INCLUDE_BOUNDARYCONDITION_H_

#include "MpsParticleSystem.h"
#include "MpsParticle.h"
#include "MpsBucket.h"

/**
 * @brief      This class describes the boundary conditions engine.
 * @details    Updates type of particle (free-surface, inner), set Bucket as Periodic or None and 
 * manages Inflow/Outflow.
 */
class MpsBoundaryCondition {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsBoundaryCondition();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsBoundaryCondition();
	
	/**
	 * @brief      Set Boundary Condition of the bucket
	 * @details    None or Periodic.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void setBucketBC(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Updates type of particle
	 * @details    Free-surface or inner particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @see        Eqs. (13) and (18) in <a href="https://doi.org/10.1016/j.cma.2010.12.001" target="_blank">Step-by-step improvement of MPS method in simulating violent free-surface motions and impact-loads</a>
	 * @see        Eqs. (17) in <a href="https://doi.org/10.1002/fld.4213" target="_blank">Fluid interface detection technique based on neighborhood particles centroid deviation (NPCD) for particle methods</a>
	 */
	void updateParticleBC(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	
protected:
	
private:
	
};

#endif // MPS_INCLUDE_BOUNDARYCONDITION_H_