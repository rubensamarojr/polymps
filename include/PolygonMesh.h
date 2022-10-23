/**
 * @defgroup   POLYGONMESH Polygon Mesh
 *
 * @brief      This file implements the functions related to a Polygonal Mesh.
 * @details    Find closest point on the mesh from a particle and corrects the particle numer density (PND) and number of neighboors.
 * Reads and writes polygon mesh STL files. Interpolates the value of wall weight function and number of neighboors due to to polygon wall.
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef EMPS_INCLUDE_POLYGONMESH_H_
#define EMPS_INCLUDE_POLYGONMESH_H_

#pragma once
#include <math.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <igl/AABB.h>
#include <igl/cotmatrix.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/readSTL.h>
//#include <igl/remove_duplicates.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/writePLY.h>
#include <igl/writeSTL.h>

//#define SHOW_FUNCT_NAME_POLY	///print the function name from any location inside a C++ function (useful for investigating programs)

// Polygon walls class

/**
 * @brief      This class describes a polygon mesh.
 * @details    Find closest point on the mesh from a particle and corrects the particle numer density (PND) and number of neighboors.
 * Reads and writes polygon mesh STL files. Interpolates the value of wall weight function and number of neighboors due to to polygon wall.
 */
class PolygonMesh {
public:
	
	/**
	 * @brief      Constructs a new instance.
	 */
	PolygonMesh();
	//PolygonMesh(const std::string& path, const int nP);
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~PolygonMesh();
	
	/**
	 * @brief      Initializes the elements of a PolygonMesh class.
	 * @param[in]  nP    Numer of particles
	 */
	void initPolygonMesh(const int nP);
	
	/**
	 * @brief      Reads a polygon mesh file.
	 * @details    STL format file
	 * @param[in]  path  The path
	 */
	void readPolygonMeshFile(const std::string& path);

	/**
	 * @brief      Writes a polygon mesh file.
	 * @details    STL format file
	 * @param[in]  mesh_ID  The mesh id
	 * @param[in]  path     The path
	 * @param[in]  iF       The file suffix number
	 */
	void writePolygonMeshFile(const int mesh_ID, const std::string& path, const int iF);
	
	/**
	 * @brief      Creates wall weight (Zij) and number of neighboors functions (numNeighWall)
	 * @param[in]  dim      The dimension
	 * @param[in]  wijType  The weight function type
	 * @param[in]  lo       The particle distance
	 * @param[in]  reL      The large effective radius
	 * @param[in]  reS      The small effective radius
	 */
	void initWijnNeigh(int dim, int wijType, double lo, double reL, double reS);
	
	/**
	 * @brief      Interpolate the value of wall weight function (Zij) from a given data
	 * @details    Returns interpolated value at x from parallel arrays ( xData, numNeighWall )
	 * Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing 
	 * boolean argument extrapolate determines behaviour beyond ends of array (if needed)
	 * @param[in]  re           The effective radius
	 * @param[in]  x            Position at space
	 * @param[in]  extrapolate  True to determines behaviour beyond ends of array (if needed)
	 * @return     Value of wall weight function (Zij)
	 * @see        Eq.(36) in <a href="https://doi.org/10.1002/fld.5083" target="_blank">Three-dimensional weakly compressible moving particle simulation coupled with geometrically nonlinear shell for hydro-elastic free-surface flows</a>
	 */
	double interpolateWij(double re, double x, bool extrapolate);
	
	/**
	 * @brief      Interpolate the value of number of neighboors function from a given data
	 * @details    Returns interpolated value at x from parallel arrays ( xData, numNeighWall )
	 * Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing 
	 * boolean argument extrapolate determines behaviour beyond ends of array (if needed)
	 * @param[in]  re           The effective radius
	 * @param[in]  x            Position at space
	 * @param[in]  extrapolate  True to determines behaviour beyond ends of array (if needed)
	 * @return     Number of neighboors
	 * @see        Eq. (37) in <a href="https://doi.org/10.1002/fld.5083" target="_blank">Three-dimensional weakly compressible moving particle simulation coupled with geometrically nonlinear shell for hydro-elastic free-surface flows</a>
	 */
	int interpolateNumNeighWall(double re, double x, bool extrapolate);
	
	// Find closest point on the mesh from a particle and corrects the PND and number of neighboors
	// Real-time collision detection, Ericson, Chapter 5 - Pg 141 - function ClosestPtPointTriangle
	// void closestPointPNDBoundary();
	
	/**
	 * @brief      Find closest point on the mesh from a particle and corrects the PND and number of neighboors
	 * @details    For Triangle meshes, the AABB tree is used to accelerate point-mesh closest point queries given a 
	 * mesh (V,F) and a query point P (Particle) find the closest point C in the triangle face or vertex
	 * @param[in]  reS2        The squared small effective radius
	 * @param[in]  reL2        The squared large effective radius
	 * @param[in]  nP          Number of particles
	 * @param[in]  wijType     The weight function type
	 * @param      Typ         The array with particle type
	 * @param[in]  fld         The fluid material identifier
	 * @param[in]  msh_id      The mesh identifier
	 * @param[in]  sta_id      The static mesh identifier
	 * @param[in]  fem_id      The deformable mesh identifier
	 * @param[in]  frw_id      The forced mesh identifier
	 * @param[in]  Pos         The array with particles position
	 * @param      wallPos     The array with wall particles position
	 * @param      mirrorPos   The array with mirrored particles position
	 * @param      riw2        The array with squared distance of particle to polygon wall
	 * @param      elementID   The array with elements id
	 * @param      meshID      The array with mesh id
	 * @param      NormalWall  The array with normal at wall
	 */
	void closestPointPNDBoundaryAABB(double reS2, double reL2, int nP, int wijType, int *Typ, int fld, int msh_id, int sta_id,
		int fem_id, int frw_id, const double *Pos, double *wallPos, double *mirrorPos, double *riw2, int *elementID, int *meshID, double *NormalWall);
	//	double *wallPos, double *mirrorPos, double *riw2, double *niw, int *numNeighw, int *elementID, std::vector<int>& particlesNearMesh);
	
	/**
	 * @brief      Update array with ID of particles near the mesh
	 * @param[in]  reS2                    The squared small effective radius
	 * @param[in]  reL2                    The squared large effective radius
	 * @param[in]  nP                      The number of particles
	 * @param[in]  wijType                 The weight function type
	 * @param      Typ                     The array with particle type
	 * @param[in]  fld                     The fluid material identifier
	 * @param[in]  riw2                    The array with squared distance of particle to polygon wall
	 * @param      niw                     The array with value of weight function due to polygon wall
	 * @param      numNeighw               The array with number of neighboors due to polygon wall
	 * @param      meshParticlesNeighbors  The vector with particles near the polygon wall
	 * @param      Nw                      Auxiliar parameter only to show particles near to polygon wall
	 */
	void updateParticlesNearPolygonMesh(double reS2, double reL2, int nP, int wijType, int *Typ, int fld, const double *riw2,
		double *niw, int *numNeighw, std::vector<int>& meshParticlesNeighbors, bool *Nw);
	
	/**
	 * @brief      Correction of velocity due the wall gradient of pressure
	 * @param[in]  timer  The timer
	 * @warning    This function is not implemented.
	 */
	void correctVelocityWallPressure(const double timer);
	
	/**
	 * @brief      Correction of velocity due the wall gradient of pressure
	 * @param[in]  timer  The timer
	 * @warning    This function is not implemented.
	 */
	void correctVelocityWallWithTensorPressure(const double timer);
	
	/**
	 * @brief      Correction of velocity due the laplacian of velocity
	 * @param[in]  timer  The timer
	 * @warning    This function is not implemented.
	 */
	void correctVelocityWallViscositySlip(const double timer);

	/**
	 * @brief      Correction of velocity due the laplacian of velocity
	 * @param[in]  timer  The timer
	 * @warning    This function is not implemented.
	 */
	void correctVelocityWallViscosityNoSlip(const double timer);
	
	/**
	 * @brief      Updates node positions
	 *
	 * @param      nodeX   The node vertex x
	 * @param      nodeY   The node vertex y
	 * @param      nodeZ   The node vertex z
	 * @param[in]  nodeDX  The node displacement dx
	 * @param[in]  nodeDY  The node displacement dy
	 * @param[in]  nodeDZ  The node displacement dz
	 */
	void updatePolygonMesh(double *nodeX, double *nodeY, double *nodeZ, const double *nodeDX, const double *nodeDY, const double *nodeDZ);
	
	/**
	 * @brief      Updates forced rigid wall
	 *
	 * @param      nodeX     The node x
	 * @param      nodeY     The node y
	 * @param      nodeZ     The node z
	 * @param      velVWall  The velocity at wall
	 * @param[in]  dt        The time step
	 * @param[in]  time      The time
	 */
	void updateForcedPolygonMesh(double *nodeX, double *nodeY, double *nodeZ, double *velVWall, const double dt, const double time);

	igl::AABB<Eigen::MatrixXd,3> treeMesh;
	Eigen::MatrixXd meshVertices;	///< Mesh Vertices
	Eigen::MatrixXi meshFaces;		///< Mesh Faces ID
	Eigen::MatrixXd meshNormals;	///< Mesh Normals
	Eigen::MatrixXd NV;				///< New mesh vertices. After remove duplicated vertices
	Eigen::MatrixXi NF;				///< New mesh faces ID. After remove duplicated vertices
	Eigen::MatrixXd NNormals;		///< New normal vertices. After remove duplicated vertices
	//Eigen::MatrixXd FN;			///< New vertices forces. After remove duplicated vertices
	// Describes containers of neighbor particles of the mesh
	//using meshParticlesNeighbors = std::vector<int>;
	//std::vector<int> meshParticlesNeighbors;
	//std::vector<int> meshParticles;
	//std::vector<int> particlesNearMesh;

private:
	Eigen::MatrixXd temporary_position;				///< Temporary particle position vector
	//Eigen::MatrixXd position_transpose;			///< Transpose of particle position vector
	//Eigen::MatrixXd wall_temporary_position;		///< Temporary wall point position
	//Eigen::Matrix3Xd mirror_temporary_position;	///< Temporary mirror particle position
	//Eigen::Matrix3Xd Xim_position;				///< Temporary mirror particle position
	//Eigen::Matrix3Xd Xwall_position;				///< Wall point position
	//Eigen::Matrix3Xd r_iwall;						///< Vector i-wall
	Eigen::MatrixXd xWall;							///< Auxiliar wall point position
	Eigen::VectorXd sqr_distance;					///< unsigned squared distance
	Eigen::VectorXi element_indice;					///< Element indice
	std::vector<double> xDataPND;					///< Particle-wall position
	std::vector<double> xDataNeigh;					///< Particle-wall position
	std::vector<double> Zij;						///< Wall weight function
	std::vector<int> nNeighWall;					///< Number of neighboor
};

#endif // EMPS_INCLUDE_POLYGONMESH_H_