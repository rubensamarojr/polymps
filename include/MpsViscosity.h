/**
 * @defgroup   MPSVISCOSITY Viscosity
 *
 * @brief      This file implements calculations related to the fluid viscosity.
 * @details    Compute accelerations due to Laplacian of velocity (+ gravity acceleration).
 * For non-Newtonian fluid, computes volume of fraction and viscosity interaction values.
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef MPS_INCLUDE_VISCOSITY_H_
#define MPS_INCLUDE_VISCOSITY_H_

#include "MpsParticleSystem.h"
#include "MpsParticle.h"
#include "MpsBucket.h"

/**
 * @brief      This class describes the mps viscosity engine.
 * @details    Compute accelerations due to Laplacian of velocity (+ gravity acceleration).
 * For non-Newtonian fluid, computes volume of fraction and viscosity interaction values.
 */
class MpsViscosity {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsViscosity();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsViscosity();
	
	/**
	 * @brief      Calculates accelerations due to Laplacian of velocity
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcViscosity(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates the volume of fraction if second fluid (phase II) is in the mixture.
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcVolumeFraction(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Calculates viscosity interaction values for "real" fluid particles
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        <a href="https://doi.org/10.1016/j.compfluid.2018.06.023" target="_blank">Meshfree particle numerical modelling of sub-aerial and submerged landslides</a>
	 */
	void calcViscosityInteractionVal(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates the viscosity interaction values.
	 * @details    Contribution of the closest polygon wall. Free-slip boundary condition.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcWallSlipViscosityInteractionVal(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	// No-Slip condition. Viscosity interaction values
	
	/**
	 * @brief      Calculates the viscosity interaction values.
	 * @details    Contribution of the closest polygon wall. No-slip boundary condition.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcWallNoSlipViscosityInteractionVal(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates acceleration due to the laplacian of velocity.
	 * @details    Contribution of the closest polygon wall. Free-slip boundary condition.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcWallSlipViscosity(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	// 
	
	/**
	 * @brief      Calculates acceleration due to the laplacian of velocity.
	 * @details    Contribution of the closest polygon wall. No-slip boundary condition.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcWallNoSlipViscosity(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
protected:
	
private:
	
};

#endif // MPS_INCLUDE_VISCOSITY_H_