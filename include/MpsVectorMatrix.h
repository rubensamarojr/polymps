/**
 * @defgroup   MPSVECTORMATRIX Vector Matrix
 *
 * @brief      This file implements vector and matrix algebra operations.
 * @details    Inverse of matrix, determinant of matrix and spatial transformations.
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef MPS_INCLUDE_VECTORMATRIX_H_
#define MPS_INCLUDE_VECTORMATRIX_H_

#include "MpsParticleSystem.h"
#include "MpsParticle.h"
#include "MpsBucket.h"

/**
 * @brief      This class describes the vector/matrix engine.
 * @details    Inverse of matrix, determinant of matrix and spatial transformations.
 */
class MpsVectorMatrix {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsVectorMatrix();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsVectorMatrix();

	/**
	 * @brief      Computes the correction matrix
	 * @details    Contribution of the neighboring particles.
	 * @param      PSystem    The physical system
	 * @param      Particles  The particles data
	 * @param      Buckets    The buckets data
	 * @see        <a href="https://doi.org/10.1016/j.cma.2010.12.001" target="_blank">Step-by-step improvement of MPS method in simulating violent free-surface motions and impact-loads</a>
	 */
	void correctionMatrix(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets);
	
	/**
	 * @brief      Computes the determinant of a matrix of dimension PSystem->dim x PSystem->dim
	 * @param[in]  M11   Element 11
	 * @param[in]  M12   Element 12
	 * @param[in]  M13   Element 13
	 * @param[in]  M21   Element 21
	 * @param[in]  M22   Element 22
	 * @param[in]  M23   Element 23
	 * @param[in]  M31   Element 31
	 * @param[in]  M32   Element 32
	 * @param[in]  M33   Element 33
	 * @return     The value of the determinant of matrix
	 */
	double detMatrix(const double M11, const double M12, const double M13, const double M21, const double M22, const double M23, 
		const double M31, const double M32, const double M33);
	
	/**
	 * @brief      Computes the inverse of a matrix of dimension PSystem->dim x PSystem->dim
	 * @param      M11      Element 11
	 * @param      M12      Element 12
	 * @param      M13      Element 13
	 * @param      M21      Element 21
	 * @param      M22      Element 22
	 * @param      M23      Element 23
	 * @param      M31      Element 31
	 * @param      M32      Element 32
	 * @param      M33      Element 33
	 * @param      PSystem    The physical system
	 * @return     The value of the inverse of matrix
	 */
	int inverseMatrix(double &M11, double &M12, double &M13, double &M21, double &M22, double &M23, 
		double &M31, double &M32, double &M33, const int dimension, const double epsZero);
	
	/**
	 * @brief      Transform a 3D triangle to xy plane
	 * @param[in]  V1       Vertex 1
	 * @param[in]  V2       Vertex 2
	 * @param[in]  V3       Vertex 3
	 * @param      RM       The rotation matrix
	 * @param      PSystem    The physical system
	 * @see        See <a href="https://math.stackexchange.com/questions/856666/how-can-i-transform-a-3d-triangle-to-xy-plane" target="_blank">https://math.stackexchange.com/questions/856666/how-can-i-transform-a-3d-triangle-to-xy-plane</a>
	 */
	void transformMatrix(const double *V1, const double *V2, const double *V3, double *RM, const double epsZero);
	
	/**
	 * @brief      Transform 3D -> 2D
	 * @param      P1    The point converted from 3D to 2D
	 * @param[in]  RM    The rotation matrix
	 */
	void transform3Dto2D(double *P1, const double *RM);

	/**
	 * @brief      Transform 2D -> 3D
	 * @param      P1    The point converted from 2D to 3D
	 * @param[in]  RM    The rotation matrix
	 */
	void transform2Dto3D(double *P1, const double *RM);
	
	
protected:
	
private:
	
};

#endif // MPS_INCLUDE_VECTORMATRIX_H_