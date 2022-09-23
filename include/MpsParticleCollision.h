/**
 * @defgroup   MPSPARTICLECOLLISION Particle Collision
 *
 * @brief      This file implements calculations related to the particle Collision models.
 * @details    Adjusts the particles positions, then avoiding numerical instabilities related to particle clustering.
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef MPS_INCLUDE_PARTICLECOLLISION_H_
#define MPS_INCLUDE_PARTICLECOLLISION_H_

#include "MpsParticleSystem.h"
#include "MpsParticle.h"
#include "MpsBucket.h"

/**
 * @brief      This class describes the particle collision engine.
 * @details    Adjusts the particles positions, then avoiding numerical instabilities related to particle clustering.
 */
class MpsParticleCollision {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsParticleCollision();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsParticleCollision();
	
	/**
	 * @brief      Check collisions between particles
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The physical system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        <a href="https://doi.org/10.1016/j.cma.2010.12.001" target="_blank">Step-by-step improvement of MPS method in simulating violent free-surface motions and impact-loads</a>
	 */
	void checkParticleCollisions(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Check collisions between particles using Dynamic Particle Collision
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The physical system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        <a href="https://doi.org/10.1016/j.jcp.2021.110202" target="_blank">Enhanced weakly-compressible MPS method for violent free-surface flows: Role of particle regularization techniques</a>
	 */
	void checkDynamicParticleCollisions(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
protected:
	
private:
	
};

#endif // MPS_INCLUDE_PARTICLECOLLISION_H_