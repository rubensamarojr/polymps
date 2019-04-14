// Copyright (c) 2019 Rubens AMARO
// Distributed under the MIT License.
#include "polygon.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <chrono>


mesh::mesh(const std::string& path, const int nP) {
  //xWall = Eigen::MatrixXd::Zero(getSize(), 3);
  // If closestPointPNDBoundaryAABB is used then comment 2 lines bellow ??
	sqr_distance = Eigen::VectorXd::Zero(nP);
	element_indice = Eigen::VectorXi::Zero(nP);
  /////////////////////////////////////////////////
	mirror_temporary_position = Eigen::Matrix3Xd::Zero(3, nP);
	wall_temporary_position = Eigen::Matrix3Xd::Zero(3, nP);
	r_iwall = Eigen::Matrix3Xd::Ones(3, nP);
  
	temporary_position = Eigen::MatrixXd::Zero(nP, 3);
//	wall_temporary_position = Eigen::VectorXd::Zero(getSize());
//	Xwall_position = Eigen::Matrix3Xd::Zero(3, getSize());
//	Xim_position = Eigen::Matrix3Xd::Zero(3, getSize()); // Tanaka source
//  meshV = Eigen::MatrixXd::Zero(getSize()*3, 3);
//  meshF = Eigen::MatrixXi::Zero(getSize(), 3);
//  meshN = Eigen::MatrixXd::Zero(getSize(), 3);
//  readMeshFile(path);
//  updateParticleNumberDensity();
//  setInitialParticleNumberDensity();
//  setLaplacianLambda();
//  checkSurfaceParticles();
}

void mesh::readMeshFile(const std::string& path) {
  // Load a mesh
  igl::readSTL(path, meshVertices, meshFaces, meshNormals);
  // Find the bounding box
  Eigen::Vector3d m = meshVertices.colwise().minCoeff();
  Eigen::Vector3d M = meshVertices.colwise().maxCoeff();
  std::cout << "m: " << m << std::endl;
  std::cout << "M: " << M << std::endl;
  // Static mesh V, F
  treeMesh.init(meshVertices, meshFaces);
}

// Returns interpolated value at x from parallel arrays ( xData, yData )
// Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
// boolean argument extrapolate determines behaviour beyond ends of array (if needed)
double mesh::interpolate(int dim, int wijType, double re, double x, bool extrapolate)
{
   	// Ratio rij/re
   	std::vector<double> xData = {2.000000e-02,4.000000e-02,6.000000e-02,8.000000e-02,1.000000e-01,1.200000e-01,1.400000e-01,1.600000e-01,1.800000e-01,2.000000e-01,
   								2.200000e-01,2.400000e-01,2.600000e-01,2.800000e-01,3.000000e-01,3.200000e-01,3.400000e-01,3.600000e-01,3.800000e-01,4.000000e-01,
   								4.200000e-01,4.400000e-01,4.600000e-01,4.800000e-01,5.000000e-01,5.200000e-01,5.400000e-01,5.600000e-01,5.800000e-01,6.000000e-01,
   								6.200000e-01,6.400000e-01,6.600000e-01,6.800000e-01,7.000000e-01,7.200000e-01,7.400000e-01,7.600000e-01,7.800000e-01,8.000000e-01,
   								8.200000e-01,8.400000e-01,8.600000e-01,8.800000e-01,9.000000e-01,9.200000e-01,9.400000e-01,9.600000e-01,9.800000e-01,1 };
   	/*
   	// Mesh surface at center of wall particle
   	// Weight function (Zij)
   	std::vector<double> yData;
   	if (condition_.dimension == 2) {
   		yData= {5.324775e+01,2.807619e+01,1.958352e+01,1.526405e+01,1.260999e+01,1.078801e+01,9.441093e+00,8.390885e+00,7.538466e+00,6.824641e+00,
				6.211883e+00,5.675301e+00,5.197765e+00,4.767128e+00,4.374544e+00,4.022795e+00,3.701045e+00,3.402413e+00,3.123826e+00,2.862821e+00,
				2.646657e+00,2.448958e+00,2.262829e+00,2.087143e+00,1.920942e+00,1.763403e+00,1.629744e+00,1.506469e+00,1.389248e+00,1.277630e+00,
				1.171215e+00,1.069645e+00,9.725965e-01,8.797768e-01,7.909196e-01,7.057820e-01,6.241412e-01,5.457926e-01,4.705477e-01,3.982322e-01,
				3.286849e-01,2.617561e-01,1.973067e-01,1.363636e-01,1.111111e-01,8.695652e-02,6.382979e-02,4.166667e-02,2.040816e-02,0};
   	}
   	else {
   		yData= {5.324775e+01,2.807619e+01,1.958352e+01,1.526405e+01,1.260999e+01,1.078801e+01,9.441093e+00,8.390885e+00,7.538466e+00,6.824641e+00,
				6.211883e+00,5.675301e+00,5.197765e+00,4.767128e+00,4.374544e+00,4.022795e+00,3.701045e+00,3.402413e+00,3.123826e+00,2.862821e+00,
				2.646657e+00,2.448958e+00,2.262829e+00,2.087143e+00,1.920942e+00,1.763403e+00,1.629744e+00,1.506469e+00,1.389248e+00,1.277630e+00,
				1.171215e+00,1.069645e+00,9.725965e-01,8.797768e-01,7.909196e-01,7.057820e-01,6.241412e-01,5.457926e-01,4.705477e-01,3.982322e-01,
				3.286849e-01,2.617561e-01,1.973067e-01,1.363636e-01,1.111111e-01,8.695652e-02,6.382979e-02,4.166667e-02,2.040816e-02,0};
	}
	*/
	// Mesh surface tangent of wall particle
   	// Weight function (Zij)
   	std::vector<double> yData;
   	yData.clear();
   	if (dim == 2) {
   		if (wijType == 0){
   			yData={5.241062e+00,4.806382e+00,4.410487e+00,4.054778e+00,3.730637e+00,3.429948e+00,3.149566e+00,2.886979e+00,2.666136e+00,2.467269e+00,
   				2.280089e+00,2.103452e+00,1.936385e+00,1.778054e+00,1.641818e+00,1.517942e+00,1.400164e+00,1.288030e+00,1.181136e+00,1.079119e+00,
   				9.816527e-01,8.884420e-01,7.992183e-01,7.137363e-01,6.317716e-01,5.531178e-01,4.775850e-01,4.049977e-01,3.351934e-01,2.680213e-01,
   				2.033415e-01,1.410236e-01,1.134677e-01,8.921162e-02,6.598985e-02,4.373757e-02,2.239533e-02,1.908397e-03,0.000000e+00,0.000000e+00,
   				0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00};
		}
   		else if (wijType = 1){}
   		else if (wijType = 2){}
   	}
   	else {
   		if (wijType == 0){
   			// wij = re/r - 1
   			yData={8.806738e+00,8.181128e+00,7.606205e+00,7.077553e+00,6.583733e+00,6.113077e+00,5.662959e+00,5.231327e+00,4.868527e+00,4.537055e+00,
   				4.218042e+00,3.910643e+00,3.614143e+00,3.327929e+00,3.065546e+00,2.817415e+00,2.577400e+00,2.345152e+00,2.120345e+00,1.902679e+00,
   				1.691872e+00,1.487658e+00,1.289784e+00,1.098012e+00,9.121101e-01,7.871401e-01,6.699803e-01,5.570121e-01,4.480352e-01,3.428608e-01,
   				2.413113e-01,1.432186e-01,1.134677e-01,8.921162e-02,6.598985e-02,4.373757e-02,2.239533e-02,1.908397e-03,0.000000e+00,0.000000e+00,
   				0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00};
		}
		else if (wijType == 1){
			// wij = re/r + r/re - 2
			yData={4.283453e+00,3.870382e+00,3.499389e+00,3.164078e+00,2.859517e+00,2.582070e+00,2.328922e+00,2.097837e+00,1.886323e+00,1.691247e+00,
				1.511369e+00,1.345706e+00,1.193412e+00,1.053747e+00,9.258529e-01,8.086067e-01,7.014818e-01,6.040244e-01,5.158145e-01,4.364597e-01,
				3.655917e-01,3.028630e-01,2.479442e-01,2.005228e-01,1.603010e-01,1.262202e-01,9.705254e-02,7.255457e-02,5.250003e-02,3.667598e-02,
				2.488189e-02,1.692868e-02,1.156289e-02,7.306856e-03,4.085086e-03,1.832813e-03,4.905643e-04,3.635042e-06,0.000000e+00,0.000000e+00,
				0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00};
		}
		else if (wijType == 2){
			// wij = re/r - r/re
			yData={1.333002e+01,1.249187e+01,1.171302e+01,1.099103e+01,1.030795e+01,9.644084e+00,8.996997e+00,8.364818e+00,7.850731e+00,7.382863e+00,
				6.924716e+00,6.475580e+00,6.034875e+00,5.602110e+00,5.205239e+00,4.826222e+00,4.453319e+00,4.086279e+00,3.724876e+00,3.368899e+00,
				3.018153e+00,2.672453e+00,2.331625e+00,1.995500e+00,1.663919e+00,1.448060e+00,1.242908e+00,1.041470e+00,8.435703e-01,6.490457e-01,
				4.577407e-01,2.695085e-01,2.153724e-01,1.711164e-01,1.278946e-01,8.564234e-02,4.430009e-02,3.813159e-03,0.000000e+00,0.000000e+00,
				0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00};
		}
		else if (wijType == 3){
			// wij = (1 - r/re)^3
			yData={9.067981e-01,8.370517e-01,7.701090e-01,7.060559e-01,6.449566e-01,5.868508e-01,5.317509e-01,4.796413e-01,4.304862e-01,3.843025e-01,
				3.411226e-01,3.009548e-01,2.637830e-01,2.295669e-01,1.982462e-01,1.697668e-01,1.440694e-01,1.210741e-01,1.006813e-01,8.277293e-02,
				6.721277e-02,5.384792e-02,4.250950e-02,3.301353e-02,2.516181e-02,1.875377e-02,1.363023e-02,9.640624e-03,6.628462e-03,4.431749e-03,
				2.883400e-03,1.811604e-03,1.058238e-03,5.494491e-04,2.372314e-04,7.358514e-05,1.051031e-05,6.910701e-09,0.000000e+00,0.000000e+00,
				0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00};
		}
   	}
	
   	int sizeData = xData.size();
   	double rij_re = x/re;

   	int i = 0;                                                                  // find left end of interval for interpolation
   	if ( rij_re >= xData[sizeData - 2] )										// special case: beyond right end
   	{
   		i = sizeData - 2;
   	}
   	else
   	{
		while ( rij_re > xData[i+1] ) i++;
   	}
   	double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
   	if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
   	{
		if ( rij_re < xL ) yR = yL;
		if ( rij_re > xR ) yL = yR;
   	}

   	double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

   	return yL + dydx * ( rij_re - xL );                                         // linear interpolation
}

// Find closest point on the mesh from a particle and corrects the PND and number of neighboors
// Libigl
void mesh::closestPointPNDBoundaryAABB(int dim, double re2, int nP, int *Typ, int fld, double *Pos, 
	double *wallPos, double *mirrorPos, double *niw, int *numNeigh, std::vector<int>& particlesNearMesh) {

  	// MPS -> libigl
#pragma omp parallel for
	for(int i=0;i<nP;i++){
		if(Typ[i] == fld){
//			temporary_position.row(i).x() = Pos[i*3  ];
//			temporary_position.row(i).y() = Pos[i*3+1];
//			temporary_position.row(i).z() = Pos[i*3+2];

			temporary_position(i,0) = Pos[i*3  ];
			temporary_position(i,1) = Pos[i*3+1];
			temporary_position(i,2) = Pos[i*3+2];
		}
	}
  
  //position_transpose = position.transpose(); // Tanaka - position instead of temporary !!!

  // Find closest point Xwall_temporary_position in element "element_indice"
  //igl::point_mesh_squared_distance(temporary_position,meshVertices,meshFaces,sqr_distance,element_indice,xWall);
  // Static mesh
  treeMesh.squared_distance(meshVertices,meshFaces,temporary_position,sqr_distance,element_indice,xWall);
  
  // libigl -> MPS
#pragma omp parallel for
	for(int i=0;i<nP;i++){
		if(Typ[i] == fld){
			// Point on the mesh
			wallPos[i*3  ] = xWall.row(i).x();
			wallPos[i*3+1] = xWall.row(i).y();
			wallPos[i*3+2] = xWall.row(i).z();
			// Mirror particle position Xm = Xi + 2*(Xw - Xi)
			mirrorPos[i*3  ] = Pos[i*3  ] + 2*(wallPos[i*3  ] - Pos[i*3  ]);
			mirrorPos[i*3+1] = Pos[i*3+1] + 2*(wallPos[i*3+1] - Pos[i*3+1]);
			mirrorPos[i*3+2] = Pos[i*3+2] + 2*(wallPos[i*3+2] - Pos[i*3+2]);
		}
	}
  // Point on the mesh
//  wall_temporary_position = xWall.transpose();
  // Mirror particle position Xm = Xi + 2*(Xw - Xi)
//  mirror_temporary_position = temporary_position + 2*(wall_temporary_position - temporary_position);
  // Vector distance between i and wall
//  r_iwall = temporary_position - wall_temporary_position;

	// PND and Neighboorhood correction
  	//voxel_ratio.setZero();
  	particlesNearMesh.clear();
// The following code for example fills std::vectors in parallel and then combines them in the end. 
// As long as your main loop/fill function is the bottleneck this should work well in general and be thread safe.
// https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector/18671256#18671256
#pragma omp parallel
	{
	  	std::vector<int> particlesNearMesh_private;
#pragma omp for nowait //fill vec_private in parallel
		for(int i=0;i<nP;i++){
			niw[i] = 0;
			if(Typ[i] == fld){
				// Compute wall weight function Z(xij)
				//if(i==6){
				//	printf("x:%lf y:%lf z:%lf x:%lf y:%lf z:%lf\n", Pos[i*3],Pos[i*3+1],Pos[i*3+2],mirrorPos[i*3],mirrorPos[i*3+1],mirrorPos[i*3+2]);
				//	double dd = sqrt((pow(0.5*(Pos[i*3]-mirrorPos[i*3]),2)+pow(0.5*(Pos[i*3+1]-mirrorPos[i*3+1]),2)+pow(0.5*(Pos[i*3+2]-mirrorPos[i*3+2]),2)));
				//	printf("dist: %f %f\n",sqrt(sqr_distance(i)), dd);
				//}
				double x2 = sqr_distance(i);
				//printf("x: %f \n",x);
				if (x2 < re2) {

					double x = sqrt(x2);
					double re = sqrt(re2);
					// Add particle ID
					particlesNearMesh_private.push_back(i);
					//particlesNearMesh.push_back(i);

					// PND due wall weight Z(xij)
					niw[i] = interpolate(dim, 0, re, x, true);
					//printf("niw: %f\n",niw[i]);
					// PND due wall weight Z(xij) to Gradient
					//niw2[i] = interpolate(dim, 2, re, x, true);

					/*
					if (dim == 2) {
				      	// Mesh surface tangent of wall particle
				      	if (x <= 0.1403*re)
				          numNeigh[i] += 8;
				        else if (x <= 0.3466*re)
				          numNeigh[i] += 6;
				        else if (x <= 0.6*re)
				          numNeigh[i] += 4;
				        else if (x <= 1.3466*re)
				          numNeigh[i] += 3;
				        else if (x <= 1.6*re)
				          numNeigh[i] += 1;
			  		}
			    	else {
				    	// Mesh surface tangent of wall particle
				      	if (x <= 0.0524*re)
				          numNeigh[i] += 22;
				        else if (x <= 0.1403*re)
				          numNeigh[i] += 18;
				        else if (x <= 0.3466*re)
				          numNeigh[i] += 14;
				        else if (x <= 0.6*re)
				          numNeigh[i] += 10;
				      	else if (x <= 1.0524*re)
				          numNeigh[i] += 9;
				        else if (x <= 1.3466*re)
				          numNeigh[i] += 5;
				        else if (x <= 1.6*re)
				          numNeigh[i] += 1;
			    	}
			    	*/
			    }
			}
		}
#pragma omp critical
    	particlesNearMesh.insert(particlesNearMesh.end(), particlesNearMesh_private.begin(), particlesNearMesh_private.end());
	}
}