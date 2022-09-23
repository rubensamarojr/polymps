/**
 * @defgroup   MPSPNDNEIGH PND, Neighbors and Deviation
 *
 * @brief      This file implements calculations related to the PND, number of neighbors and weighted deviation vector.
 * @details    Computes Particle Number Density (PND), number of neighbors and Neighborhood Particles Centroid Deviation (NPCD).
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef MPS_INCLUDE_PNDNEIGH_H_
#define MPS_INCLUDE_PNDNEIGH_H_

#include "MpsParticleSystem.h"
#include "MpsParticle.h"
#include "MpsBucket.h"

/**
 * @brief      This class describes the PND, Number of neighbors and NPCD engines.
 * @details    Computes Particle Number Density (PND), number of neighbors and Neighborhood Particles Centroid Deviation (NPCD).
 */
class MpsPndNeigh {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsPndNeigh();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsPndNeigh();
	
	/**
	 * @brief      Set initial PND, number of neighbors and weighted deviation.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void setInitialPndNumberOfNeighNPCD(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Calculates weighted deviation.
	 * @details    Contribution of the closest polygon wall.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcWallNPCD(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Calculates PND, number of neighbors and weighted deviation.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void calcPndnNeighNPCD(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates PND based on continuity equation including a diffusive term.
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        Eq. (19) in <a href="https://doi.org/10.1002/fld.5083" target="_blank">Three-dimensional weakly compressible moving particle simulation coupled with geometrically nonlinear shell for hydro-elastic free-surface flows</a>
	 */
	void calcPndDiffusiveTerm(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates PND based on continuity equation including a diffusive term.
	 * @details    Contribution of the closest polygon wall. Free-slip boundary condition.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        Eq. (35) and (42) in <a href="https://doi.org/10.1002/fld.5083" target="_blank">Three-dimensional weakly compressible moving particle simulation coupled with geometrically nonlinear shell for hydro-elastic free-surface flows</a>
	 */
	void calcWallSlipPndDiffusiveTerm(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Calculates PND based on continuity equation including a diffusive term.
	 * @details    Contribution of the closest polygon wall. No-slip boundary condition.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        Eq. (35) and (45) in <a href="https://doi.org/10.1002/fld.5083" target="_blank">Three-dimensional weakly compressible moving particle simulation coupled with geometrically nonlinear shell for hydro-elastic free-surface flows</a>
	 */
	void calcWallNoSlipPndDiffusiveTerm(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Updates diffusive PND.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void updatePnd(MpsParticle *Particles);
	
	/**
	 * @brief      Mean diffusive PND at wall, dummy or free-surface particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void meanPndParticlesWallDummySurface(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Mean PND computed as the sum of the weight function.
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void meanPnd(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

	/**
	 * @brief      Mean PND computed as the sum of the weight function.
	 * @details    Contribution of the closest polygon wall.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void meanWallPnd(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Mean PND computed as the sum of the weight function.
	 * @details    Contribution of the fluid neighboring particles.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 */
	void meanNeighFluidPnd(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);

protected:
	
private:
	
};

#endif // MPS_INCLUDE_PNDNEIGH_H_