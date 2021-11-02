// Copyright (c) 2019 Rubens AMARO
// Distributed under the MIT License.

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

// Polygon walls class
class PolygonMesh {
public:

	PolygonMesh();
	//PolygonMesh(const std::string& path, const int nP);
	virtual ~PolygonMesh();
	// Initi elements of the class
	void initPolygonMesh(const int nP);
	// Read mesh file (STL)
	void readPolygonMeshFile(const std::string& path);
	// Write mesh file (STL)
	void writePolygonMeshFile(const int mesh_ID, const std::string& path, const int iF);
	// Create wall weight and number of neighboors functions
	void initWijnNeigh(int dim, int wijType, double lo, double reL, double reS);
	// Interpolate the value of wall weight function from a given data
	double interpolateWij(double re, double x, bool extrapolate);
	// Interpolate the value of number of neighboors function from a given data
	int interpolateNumNeighWall(double re, double x, bool extrapolate);
	// Find closest point on the mesh from a particle and corrects the PND and number of neighboors
	// Real-time collision detection, Ericson, Chapter 5 - Pg 141 - function ClosestPtPointTriangle
//	void closestPointPNDBoundary();
	// For Triangle meshes, the AABB tree is used to accelerate point-mesh closest point queries given a 
	// mesh (V,F) and a query point P (Particle) find the closest point C in the triangle face or vertex
	// Libigl
	void closestPointPNDBoundaryAABB(double reS2, double reL2, int nP, int wijType, int *Typ, int fld, int msh_id, int sta_id,
	int fem_id, int frw_id, double *Pos, double *wallPos, double *mirrorPos, double *riw2, int *elementID, int *meshID, double *NormalWall);
	//	double *wallPos, double *mirrorPos, double *riw2, double *niw, int *numNeighw, int *elementID, std::vector<int>& particlesNearMesh);
	// Update vector with ID of particles near the mesh
	void updateParticlesNearPolygonMesh(double reS2, double reL2, int nP, int wijType, int *Typ, int fld, double *riw2,
	double *niw, int *numNeighw, std::vector<int>& meshParticlesNeighbors, bool *Nw);
	// Correction of velocity due the wall gradient of pressure
	void correctVelocityWallPressure(const double timer);
	void correctVelocityWallWithTensorPressure(const double timer);
	// Correction of velocity due the laplacian of velocity
	void correctVelocityWallViscositySlip(const double timer);
	void correctVelocityWallViscosityNoSlip(const double timer);
	// Update node positions
	void updatePolygonMesh(double *nodeX, double *nodeY, double *nodeZ, double *nodeDX, double *nodeDY, double *nodeDZ);
	// Update forced rigid wall
	void updateForcedPolygonMesh(double *nodeX, double *nodeY, double *nodeZ, double *velVWall, const double dt, const double time);

	igl::AABB<Eigen::MatrixXd,3> treeMesh;
	Eigen::MatrixXd meshVertices;	// Mesh Vertices
	Eigen::MatrixXi meshFaces;		// Mesh Faces ID
	Eigen::MatrixXd meshNormals;	// Mesh Normals
	Eigen::MatrixXd NV;				// New mesh vertices. After remove duplicated vertices
	Eigen::MatrixXi NF;				// New mesh faces ID. After remove duplicated vertices
	Eigen::MatrixXd NNormals;		// New normal vertices. After remove duplicated vertices
	//Eigen::MatrixXd FN;			// New vertices forces. After remove duplicated vertices
	// Describes containers of neighbor particles of the mesh
	//using meshParticlesNeighbors = std::vector<int>;
	//std::vector<int> meshParticlesNeighbors;
	//std::vector<int> meshParticles;
	//std::vector<int> particlesNearMesh;

private:
	Eigen::MatrixXd temporary_position;				// Temporary particle position vector
	//Eigen::MatrixXd position_transpose;			// Transpose of particle position vector
	//Eigen::MatrixXd wall_temporary_position;		// Temporary wall point position
	//Eigen::Matrix3Xd mirror_temporary_position;	// Temporary mirror particle position
	//Eigen::Matrix3Xd Xim_position;				// Temporary mirror particle position
	//Eigen::Matrix3Xd Xwall_position;				// Wall point position
	//Eigen::Matrix3Xd r_iwall;						// Vector i-wall
	Eigen::MatrixXd xWall;							// Auxiliar wall point position
	Eigen::VectorXd sqr_distance;					// unsigned squared distance
	Eigen::VectorXi element_indice;					// Element indice
	std::vector<double> xDataPND;					// Particle-wall position
	std::vector<double> xDataNeigh;					// Particle-wall position
	std::vector<double> Zij;						// Wall weight function
	std::vector<int> nNeighWall;					// Number of neighboor
};

#endif // EMPS_INCLUDE_POLYGONMESH_H_