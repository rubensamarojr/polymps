/**
 * @defgroup   MPSINFLOWOUTFLOW Mps Inflow Outflow
 *
 * @brief      This file maneges Inflow/Outflow boundary conditions.
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef MPS_INCLUDE_INFLOWOUTFLOW_H_
#define MPS_INCLUDE_INFLOWOUTFLOW_H_

#include "MpsParticleSystem.h"
#include "MpsParticle.h"
//#include "MpsBucket.h"

/**
 * @brief      This class describes the Inflow/Outflow engine.
 * @details    Creates or removes particles, updates array of particles, and 
 * manages Inflow/Outflow.
 */
class MpsInflowOutflow {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsInflowOutflow();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsInflowOutflow();

	/**
	 * @brief      This class describes a plane.
	 * @details    The normalized plane equation is written component wise as ax + by + cz − d = 0, 
	 * where normal = (a, b, c) and d = dot(normal,p) for a given point p on the plane
	 * @se
	 */
	class Plane {
	public:
		double a, b, c, d;	///< Normalized plane equation components ax + by + cz − d = 0
		double normal[3];	///< Normalized normal vector
	};
	/**
	 * Plane of the iNFLOW/oUTFLOW plane (interface)
	 */
	Plane Pio;

	void setPlaneFromPointNormal(MpsParticleSystem *PSystem, const int ii);

	double calcSignedDistance(const double *P1);

	/**
	 * @brief      Computes the correction matrix
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The physical system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        <a href="https://doi.org/10.1016/j.cma.2010.12.001" target="_blank">Step-by-step improvement of MPS method in simulating violent free-surface motions and impact-loads</a>
	 */
	//void correctionMatrix(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	
protected:
	
private:
	
};

#endif // MPS_INCLUDE_INFLOWOUTFLOW_H_