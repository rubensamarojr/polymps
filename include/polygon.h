// Copyright (c) 2019 Rubens AMARO
// Distributed under the MIT License.

#pragma once
#include <math.h>


#include <igl/cotmatrix.h>

#include <igl/readSTL.h>
#include <igl/writePLY.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/AABB.h>
#include <igl/writeSTL.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>

class mesh {
public:

  mesh(const std::string& path, const int nP);
  virtual ~mesh() {};
  // Read mesh file (STL)
  void readMeshFile(const std::string& path);
  // Interpolate the value of wall weight function from a given data
  double interpolateWij(int dim, int wijType, double re, double x, bool extrapolate);
  // Interpolate the value of wall first order derivative weight function from a given data
  double interpolateDwij(int dim, int wijType, double re, double x, bool extrapolate);
  // Find closest point on the mesh from a particle and corrects the PND and number of neighboors
  // Real-time collision detection, Ericson, Chapter 5 - Pg 141 - function ClosestPtPointTriangle
  void closestPointPNDBoundary();
  // For Triangle meshes, the AABB tree is used to accelerate point-mesh closest point queries given a 
  // mesh (V,F) and a query point P (Particle) find the closest point C in the triangle face or vertex
  // Libigl
  void closestPointPNDBoundaryAABB(int dim, double re2, int nP,  int wijType, int *Typ, int fld, double *Pos, 
    double *wallPos, double *mirrorPos, double *niw, int *numNeigh, std::vector<int>& meshParticlesNeighbors);
  // Correction of velocity due the wall gradient of pressure
  void correctVelocityWallPressure(const double timer);
  void correctVelocityWallWithTensorPressure(const double timer);
  // Correction of velocity due the laplacian of velocity
  void correctVelocityWallViscositySlip(const double timer);
  void correctVelocityWallViscosityNoSlip(const double timer);

  igl::AABB<Eigen::MatrixXd,3> treeMesh;
  Eigen::MatrixXd meshVertices;   // Mesh vertices
  Eigen::MatrixXi meshFaces;      // Mesh Faces
  Eigen::MatrixXd meshNormals;    // Mesh Normals

  // Describes containers of neighbor particles of the mesh
  //using meshParticlesNeighbors = std::vector<int>;
  //std::vector<int> meshParticlesNeighbors;
  //std::vector<int> meshParticles;
  //std::vector<int> particlesNearMesh;
  
private:
  Eigen::MatrixXd temporary_position; // Temporary particle position vector
  //Eigen::MatrixXd position_transpose;           // Transpose of particle position vector
  Eigen::MatrixXd wall_temporary_position;      // Temporary wall point position
  Eigen::Matrix3Xd mirror_temporary_position;      // Temporary mirror particle position
  //Eigen::Matrix3Xd Xim_position;                // Temporary mirror particle position
  //Eigen::Matrix3Xd Xwall_position;            // Wall point position
  Eigen::MatrixXd xWall;                        // Auxiliar wall point position
  Eigen::VectorXd sqr_distance;                 // unsigned squared distance
  Eigen::VectorXi element_indice;               // Element indice
  Eigen::Matrix3Xd r_iwall;                     // Vector i-wall
};
