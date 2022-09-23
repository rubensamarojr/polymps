/**
 * @defgroup   MPSBUCKET Neighboring Buckets
 *
 * @brief      This file implements Neighboring Buckets.
 * @details    Updates the neighboring particle candidates list.
 *
 * @author     Rubens Amaro
 * @date       2022
 */


#ifndef MPS_INCLUDE_BUCKET_H_
#define MPS_INCLUDE_BUCKET_H_

#include "MpsParticleSystem.h"
#include "MpsParticle.h"

/**
 * @brief      This class describes the Bucket engines.
 * @details    Updates the neighboring particle candidates list.
 */
class MpsBucket {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsBucket();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsBucket();
	
	/**
	 * @brief      Allocation of memory for buckets data.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @see        <a href="https://doi.org/10.1007/s40571-015-0059-2" target="_blank">Performance improvements of differential operators code for MPS method on GPU</a>
	 */
	void allocateMemory(MpsParticleSystem *PSystem, MpsParticle *Particles);
	
	/**
	 * @brief      Calculates the physical domain limits and adjusts the buckets values.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void calcDomainLimits(MpsParticleSystem *PSystem, MpsParticle *Particles);
	
	/**
	 * @brief      Update particle ID's in buckets.
	 * @details    Updates the neighboring particle candidates list.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @see        <a href="https://doi.org/10.1007/s40571-015-0059-2" target="_blank">Performance improvements of differential operators code for MPS method on GPU</a>
	 */
	void updateParticlesID(MpsParticleSystem *PSystem, MpsParticle *Particles);
	
	/**
	 * @brief      Return the bucket coordinates for particle "i"
	 * @param      bx       Bucket coordinate X
	 * @param      by       Bucket coordinate Y
	 * @param      bz       Bucket coordinate Z
	 * @param[in]  rxi      Particle position X
	 * @param[in]  ryi      Particle position Y
	 * @param[in]  rzi      Particle position Z
	 * @param      PSystem    The particle system
	 * @see        Eq. (13) in <a href="https://doi.org/10.1007/s40571-015-0059-2" target="_blank">Performance improvements of differential operators code for MPS method on GPU</a>
	 */
	void bucketCoordinates(int &bx, int &by, int &bz,
		const double rxi, const double ryi, const double rzi, MpsParticleSystem *PSystem);

protected:
	
private:

	// double bucketSide;			///< Length of one bucket side
	// double bucketSide2;			///< Length of one bucket side to square
	// double invBucketSide;		///< Inverse of length of one bucket side
	
	// Copy data from periodic buckets to border buckets
	void copyDataBetweenBuckets(const int b, MpsParticleSystem *PSystem, MpsParticle *Particles);
	
};

#endif // MPS_INCLUDE_BUCKET_H_