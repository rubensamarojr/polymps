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
		int ID;				///< Plane ID
		double press;		///< Plane inOuflow cte press
		double velocity[3];	///< Plane inOuflow cte velocity vector
	};
	/**
	 * Plane of the iNFLOW/oUTFLOW plane (interface)
	 */
	Plane Pio;

	/**
	 * @brief      Constructs a plane given a point and normal
	 * @details    Constructs a plane in 3D space specifying a single point on the plane and 
	 * the surface normal pointing toward the fluid domain
	 * @param      PSystem  The p system
	 * @param[in]  ii       The new value
	 */
	void setPlaneFromPointNormal(MpsParticleSystem *PSystem, const int ii);

	/**
	 * @brief      Calculates the euclidean signed distance between point P = (Px, Py, Pz) and the plane.
	 * @param[in]  Px    X component of Point P
	 * @param[in]  Py    Y component of Point P
	 * @param[in]  Pz    Z component of Point P
	 *
	 * @return     The signed distance.
	 */
	double calcSignedDistance(const double Px, const double Py, const double Pz);


	// Impose motion to particles
	void imposeMotionParticles(MpsParticleSystem *PSystem, MpsParticle *Particles);

	// Check particles in the Inflow/Outflow region
	void checkInOutflowParticles(MpsParticleSystem *PSystem, MpsParticle *Particles);

	// Clear variables only of the inOutflow (IO) particles
	void clearVariablesInOutflowParticles(MpsParticleSystem *PSystem, MpsParticle *Particles);

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