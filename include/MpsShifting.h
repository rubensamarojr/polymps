/**
 * @defgroup   MPSSHIFTING Particle Shifting
 *
 * @brief      This file implements calculations related to the particle Shifting models.
 * @details    Adjusts the particles positions or velocities to improve the uniformity of the particle distribution.
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef MPS_INCLUDE_SHIFTING_H_
#define MPS_INCLUDE_SHIFTING_H_

#include "MpsParticleSystem.h"
#include "MpsParticle.h"
#include "MpsBucket.h"

/**
 * @brief      This class describes the particle shifting engine.
 * @details    Adjusts the particles positions or velocities to improve the uniformity of the particle distribution.
 */
class MpsShifting {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsShifting();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsShifting();
	
	/**
	 * @brief      Adjustment of particle velocity.
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        Eq. (21) in <a href="https://doi.org/10.1016/j.compfluid.2016.07.014" target="_blank">Improvements for accuracy and stability in a weakly-compressible particle method</a>
	 */
	void calcShifting(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Calculates normal vector on the fluid particle.
	 * @details    The surface normal is a weighted summation of direction vectors. Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        Eq. (39) in <a href="https://onlinelibrary.wiley.com/doi/full/10.1002/nme.5844" target="_blank">An accurate and stable multiphase moving particle semi-implicit method based on a corrective matrix for all particle interaction models</a>
	 */
	void calcNormalParticles(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Adjustment of particle velocity.
	 * @details    Contribution of the closest polygon wall.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        We adpated to polygon wall the Eq. (21) from <a href="https://doi.org/10.1016/j.compfluid.2016.07.014" target="_blank">Improvements for accuracy and stability in a weakly-compressible particle method</a>
	 */
	void calcWallShifting(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Calculates normal vector on the fluid particle.
	 * @details    The surface normal is a weighted summation of direction vectors. Contribution of the closest polygon wall.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        We adpated to polygon wall the Eq. (39) from <a href="https://onlinelibrary.wiley.com/doi/full/10.1002/nme.5844" target="_blank">An accurate and stable multiphase moving particle semi-implicit method based on a corrective matrix for all particle interaction models</a>
	 */
	void calcWallNormalParticles(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Calculates the Concentration and Gradient of concentration on the fluid particle, and adjusts its position.
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        Eq. (15) and (16) in <a href="https://doi.org/10.1002/fld.5083" target="_blank">Three-dimensional weakly compressible moving particle simulation coupled with geometrically nonlinear shell for hydro-elastic free-surface flows</a>
	 */
	void calcConcAndConcGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates the Concentration and Gradient of concentration on the fluid particle, and adjusts its position.
	 * @details    Contribution of the closest polygon wall.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        Eq. (33) in <a href="https://doi.org/10.1002/fld.5083" target="_blank">Three-dimensional weakly compressible moving particle simulation coupled with geometrically nonlinear shell for hydro-elastic free-surface flows</a>
	 */
	void calcWallConcAndConcGradient(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates normal vector on the fluid particle.
	 * @details    Normalized gradient of concentration.
	 * @param      Particles  The particles data
	 */
	void calcNormalConcentration(MpsParticle *Particles);
	

protected:
	
private:
	
};

#endif // MPS_INCLUDE_SHIFTING_H_