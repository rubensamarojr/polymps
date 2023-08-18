/**
 * @defgroup   MPSPARTICLEVELPOS Particle Velocity Position
 *
 * @brief      This file implements calculations related to the particle velocity and position.
 * @details    Updates particles velocities and positions. 
 *
 * @author     Rubens Amaro
 * @date       2022
 */


#ifndef MPS_INCLUDE_PARTICLEVELPOS_H_
#define MPS_INCLUDE_PARTICLEVELPOS_H_

#include "MpsParticleSystem.h"
#include "MpsParticle.h"
#include "MpsBucket.h"

/**
 * @brief      This class describes the particle collision engine.
 * @details    Updates particles velocities and positions.
 */
class MpsParticleVelPos {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsParticleVelPos();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsParticleVelPos();

	/**
	 * @brief      Prediction of particle velocity and position
	 * @param      PSystem    The physical system
	 * @param      Particles  The particles data
	 */
	void predictionVelocityPosition(MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Correction of particle velocity and position
	 * @details    Besides velocity and position for each particle, the maximum velocity (PSystem->velMax) 
	 * of the fluid domain is calculated herein.
	 * @param      PSystem    The physical system
	 * @param      Particles  The particles data
	 */
	void correctionVelocityPosition(MpsParticleSystem *PSystem, MpsParticle *Particles);

	// 
	/**
	 * @brief      Update velocity at wall and dummy particles
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The physical system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        <a href="https://doi.org/10.1016/j.cma.2010.12.001" target="_blank">Step-by-step improvement of MPS method in simulating violent free-surface motions and impact-loads</a>
	 */
	void updateVelocityParticlesWallDummy(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	
protected:
	
private:
	
};

#endif // MPS_INCLUDE_PARTICLEVELPOS_H_