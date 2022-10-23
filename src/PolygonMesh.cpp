// Copyright (c) 2019 Rubens AMARO
// Distributed under the MIT License.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include "PolygonMesh.h"

// Constructor declaration
PolygonMesh::PolygonMesh()
{
}
// Destructor declaration
PolygonMesh::~PolygonMesh()
{
}

void PolygonMesh::initPolygonMesh(const int nP) {
	//xWall = Eigen::MatrixXd::Zero(getSize(), 3);
	// If closestPointPNDBoundaryAABB is used then comment 2 lines bellow ??
	sqr_distance = Eigen::VectorXd::Zero(nP);
	element_indice = Eigen::VectorXi::Zero(nP);
/////////////////////////////////////////////////
//mirror_temporary_position = Eigen::Matrix3Xd::Zero(3, nP);
//wall_temporary_position = Eigen::Matrix3Xd::Zero(3, nP);
//r_iwall = Eigen::Matrix3Xd::Ones(3, nP);

	temporary_position = Eigen::MatrixXd::Zero(nP, 3);
//	wall_temporary_position = Eigen::VectorXd::Zero(getSize());
//	Xwall_position = Eigen::Matrix3Xd::Zero(3, getSize());
//	Xim_position = Eigen::Matrix3Xd::Zero(3, getSize()); // Tanaka source
//	meshV = Eigen::MatrixXd::Zero(getSize()*3, 3);
//	meshF = Eigen::MatrixXi::Zero(getSize(), 3);
//	meshN = Eigen::MatrixXd::Zero(getSize(), 3);
//	readMeshFile(path);
//	updateParticleNumberDensity();
//	setInitialParticleNumberDensity();
//	setLaplacianLambda();
//	checkSurfaceParticles();
}

void PolygonMesh::readPolygonMeshFile(const std::string& path) {

	FILE * filename = fopen(path.c_str(),"rb");
	if(NULL==filename)
	{
		fprintf(stderr,"IOError: %s could not be opened...\n", path.c_str());
		return;
	}
	// Load a mesh
	bool success = igl::readSTL(filename, meshVertices, meshFaces, meshNormals);
	if(success == false)
	{
		fprintf(stderr,"IOError: STL mesh %s could not be read\n", path.c_str());
		return;
	}

	// Find the bounding box
	Eigen::Vector3d m = meshVertices.colwise().minCoeff();
	Eigen::Vector3d M = meshVertices.colwise().maxCoeff();
	std::cout << "Minimum mesh limit: " << m.transpose() << std::endl;
	std::cout << "Maximum mesh limit: " << M.transpose() << std::endl;
	// Static mesh V, F
	//treeMesh.init(meshVertices, meshFaces);

	for(int mf=0;mf<meshFaces.rows();mf++) {
		int node0 = meshFaces(mf,0);
		int node1 = meshFaces(mf,1);
		int node2 = meshFaces(mf,2);

		Eigen::Vector3d a0, a1, a2;
		a0 << meshVertices(node0,0) , meshVertices(node0,1) , meshVertices(node0,2);
		a1 << meshVertices(node1,0) , meshVertices(node1,1) , meshVertices(node1,2);
		a2 << meshVertices(node2,0) , meshVertices(node2,1) , meshVertices(node2,2);
		//std::cout << "mf: " << mf << std::endl;
		//std::cout << "a0: " << a0 << std::endl;
		//std::cout << "a1: " << a1 << std::endl;
		//std::cout << "a2: " << a2 << std::endl;
	}

	std::cout << " Mesh file name " << path << std::endl;
	std::cout << " original  mesh containts " << meshVertices.rows() << " vertices and " << meshFaces.rows() << " faces" << std::endl;

	// Remove duplicated vertices upto a uniqueness tolerance (epsilon)
	Eigen::MatrixXi SVI, SVJ;
	//Eigen::VectorXi VI;
	//igl::remove_duplicate_vertices(meshVertices, meshFaces, 0.0, SV, SVI, SVJ, SF);
	//igl::remove_duplicate_vertices(meshVertices, 0.0, SV, SVI, SVJ);
	//igl::remove_duplicates(meshVertices, meshFaces, NV, NF, VI, 2.0e-15);

	// Mesh in (meshVertices, meshFaces)
	//igl::remove_duplicate_vertices(meshVertices, 2.0e-15, NV, SVI, SVJ);
	igl::remove_duplicate_vertices(meshVertices, meshFaces, 2.0e-15, NV, SVI, SVJ, NF);
	
	//Eigen::MatrixXd NormalAux;
	//NNormals = Eigen::MatrixXd::Zero(NF.rows(), 3);
	// Compute face normals via vertex position list, face list
	igl::per_face_normals(NV, NF, NNormals);
	//igl::per_face_normals(NV, NF, NormalAux);
	
	// Forces on nodes
	//FN = Eigen::MatrixXd::Zero(NF.rows(), 3);
	//igl::per_face_normals(NV, NF, FN);

	// Static mesh V, F
	treeMesh.init(NV, NF);
	//std::cout << "nVertices: " << meshVertices.rows() << std::endl;
}

void PolygonMesh::writePolygonMeshFile(const int mesh_ID, const std::string& path, const int iF) {

	char numstr[21], mesh_IDstr[21]; // enough to hold all numbers up to 64-bits
	sprintf(numstr, "%05d", iF);
	sprintf(mesh_IDstr, "%02d", mesh_ID);
	std::string outout_filename = path + "/mesh" + mesh_IDstr + "_" + numstr + ".stl";

//  #include <sstream>

	//std::ostringstream outout_filename;
	//outout_filename << path << "select logged from login where id = " << iF;
	//std::string query(outout_filename.str());

	//char outout_filename[256];

	//sprintf(outout_filename, path"/mesh%05d",iF);

	//sprintf(outout_filename, OUT_FOLDER"/output%05d.vtu",iF);
	//Eigen::Transform<double,3,Eigen::Affine> t = Eigen::Transform<double,3,Eigen::Affine>::Identity();
	//t.rotate(Eigen::AngleAxisd(0.5 * M_PI * timer.getCurrentDeltaTime(), Eigen::Vector3d::UnitZ()));
	//t.translate(Eigen::Vector3d(0, 0, 0)) ;
	//meshVertices = (t * meshVertices.transpose().colwise().homogeneous()).transpose();
	//std::cout << "mVertices: " << meshVertices << " \n";

	//std::string output_index = ;
	//std::string output_index2 = (boost::format(path) % output_index).str();

	//std::string output_index3 = (boost::format(path + "output_%1%.stl") % output_index).str();

	igl::writeSTL(outout_filename,NV,NF,NNormals);
	//igl::writeSTL(outout_filename,NV,NF); 
}

// Update node positions based on Finite Element Method computation - NOT WORKING !!!
void PolygonMesh::updatePolygonMesh(double *nodeX, double *nodeY, double *nodeZ, const double *nodeDX, const double *nodeDY, const double *nodeDZ) {
	// Update node positions
	int nNodes = NV.rows();
#pragma omp parallel for
	for(int nn=0;nn<nNodes;nn++) {
		nodeX[nn] += nodeDX[nn];
		nodeY[nn] += nodeDY[nn];
		nodeZ[nn] += nodeDZ[nn];
		NV(nn,0) = nodeX[nn];
		NV(nn,1) = nodeY[nn];
		NV(nn,2) = nodeZ[nn];
	}

	Eigen::MatrixXd NormalAux;
	// Compute face normals via vertex position list, face list
	igl::per_face_normals(NV, NF, NormalAux);
#pragma omp parallel for
	for(int ff=0;ff<NF.rows();ff++) {
		NNormals(ff,0) = NormalAux(ff,0);
		NNormals(ff,1) = NormalAux(ff,1);
		NNormals(ff,2) = NormalAux(ff,2);
	}
}

// Update forced motion rigid wall
void PolygonMesh::updateForcedPolygonMesh(double *nodeX, double *nodeY, double *nodeZ, double *velVWall, const double dt, const double time) {
	// Update node positions
	// Set the motion here
	// Liao 2015
	double f1 = -300.0*time*time*time + 75.0*time*time;
	double f2 = -300.0*(time+dt)*(time+dt)*(time+dt) + 75.0*(time+dt)*(time+dt);
	double vel = (f2 - f1)/dt;
	double tfim = 0.13;

	// Yilmaz 2021
	//double vel = 2.0;
	//double tfim = 0.15;

	if(time > tfim)
		vel = 0.0;

	velVWall[0] = 0.0;
	velVWall[1] = 0.0;
	velVWall[2] = vel;

	double dx = 0.0;
	double dy = 0.0;
	double dz = vel*dt;

	int nNodes = NV.rows();
#pragma omp parallel for
	for(int nn=0;nn<nNodes;nn++) {
		nodeX[nn] += dx;
		nodeY[nn] += dy;
		nodeZ[nn] += dz;
		NV(nn,0) = nodeX[nn];
		NV(nn,1) = nodeY[nn];
		NV(nn,2) = nodeZ[nn];
	}

	Eigen::MatrixXd NormalAux;
	// Compute face normals via vertex position list, face list
	igl::per_face_normals(NV, NF, NormalAux);
#pragma omp parallel for
	for(int ff=0;ff<NF.rows();ff++) {
		NNormals(ff,0) = NormalAux(ff,0);
		NNormals(ff,1) = NormalAux(ff,1);
		NNormals(ff,2) = NormalAux(ff,2);
	}

#ifdef SHOW_FUNCT_NAME_POLY
	// print the function name (useful for investigating programs)
	std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif
}

// Creation of the wall weight (Zij) and number of neighboors functions (numNeighWall)
void PolygonMesh::initWijnNeigh(int dim, int wijType, double lo, double reL, double reS) {
	double reS2 = reS*reS;
	double reL2 = reL*reL;
	// Number of points of the functions
	int n = 50;

	// Tangent wall surface
	double X0 = 0.0;	double Y0 = 0.0;	double Z0 = 0.0;

	int limiteS = ceil(reS/lo);
	int limiteL = ceil(reL/lo);
	// Position of "i" changing along the verical axis (distance from the plane wall)
	for(int i=1;i<=n;i++) {
		if(dim == 2) {
			double xi, yi, riw, riw_re;
			xi=0.0;	yi=i*reS/n;
			// Distance to wall surface
			riw = sqrt((X0-xi)*(X0-xi) + (Y0-yi)*(Y0-yi));
			// Ratio rij/re
			riw_re = riw/reS;
			xDataPND.push_back(riw_re);
			riw_re = riw/reL;
			xDataNeigh.push_back(riw_re);

			double n0 = 0.0;
			// Neighborhood
			for(int ix= -limiteS;ix<=limiteS;ix++) {
			for(int iy= -limiteS;iy<=0;iy++) {
				double xj = lo * (double)ix;
				double yj = lo * (double)(iy - 0.5);
				double rij2 = (xj-xi)*(xj-xi)+(yj-yi)*(yj-yi);
				if(rij2 < reS2) {
					if(rij2 == 0.0) continue;
					double rij = sqrt(rij2);
					if(wijType == 0)
						n0 += reS/rij - 1.0;
					else if(wijType == 1)
						n0 += reS/rij + rij/reS - 2.0;
					else if(wijType == 2)
						n0 += reS/rij - rij/reS;
					else if(wijType == 3)
						n0 += pow(1.0-rij/reS,3.0);
					else if(wijType == 4)
						n0 += pow(1.0-rij/reS,2.0);
					else
						n0 += reS/rij - 1.0;
				}
			}}
			// Particle number density due wall
			Zij.push_back(n0);

			int nNeigh = 0;
			n0 = 0.0;
			// Neighborhood
			for(int ix= -limiteL;ix<=limiteL;ix++) {
			for(int iy= -limiteL;iy<=0;iy++) {
				double xj = lo * (double)ix;
				double yj = lo * (double)(iy - 0.5);
				double rij2 = (xj-xi)*(xj-xi)+(yj-yi)*(yj-yi);
				if(rij2 < reL2) {
					if(rij2 == 0.0) continue;
					nNeigh += 1;
				}
			}}
			// Number of neighboors due wall
			nNeighWall.push_back(nNeigh);
		}
		if(dim == 3) {
			double xi, yi, zi, riw, riw_re;
			xi=0.0;	yi=i*reS/n;	zi=0.0;
			// Distance to wall surface
			riw = sqrt((X0-xi)*(X0-xi) + (Y0-yi)*(Y0-yi) + (Z0-zi)*(Z0-zi));
			// Ratio rij/re
			riw_re = riw/reS;
			xDataPND.push_back(riw_re);
			riw_re = riw/reL;
			xDataNeigh.push_back(riw_re);
			
			//std::cout << " xData[" << i << "]: " << riw << std::endl;

			double n0 = 0.0;
			// Neighborhood
			for(int ix= -limiteS;ix<=limiteS;ix++) {
			for(int iy= -limiteS;iy<=0      ;iy++) {
			for(int iz= -limiteS;iz<=limiteS;iz++) {
				double xj = lo * (double)ix;
				double yj = lo * (double)(iy - 0.5);
				double zj = lo * (double)iz;
				double rij2 = (xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zj-zi)*(zj-zi);
				if(rij2 < reS2) {
					if(rij2 == 0.0) continue;
					double rij = sqrt(rij2);
					if(wijType == 0)
						n0 += reS/rij - 1.0;
					else if(wijType == 1)
						n0 += reS/rij + rij/reS - 2.0;
					else if(wijType == 2)
						n0 += reS/rij - rij/reS;
					else if(wijType == 3)
						n0 += pow(1.0-rij/reS,3.0);
					else if(wijType == 4)
						n0 += pow(1.0-rij/reS,2.0);
					else
						n0 += reS/rij - 1.0;
				}
			}}}
			// Particle number density due wall
			Zij.push_back(n0);

			int nNeigh = 0;
			n0 = 0.0;
			// Neighborhood
			for(int ix= -limiteL;ix<=limiteL;ix++) {
			for(int iy= -limiteL;iy<=0      ;iy++) {
			for(int iz= -limiteL;iz<=limiteL;iz++) {
				double xj = lo * (double)ix;
				double yj = lo * (double)(iy - 0.5);
				double zj = lo * (double)iz;
				double rij2 = (xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zj-zi)*(zj-zi);
				if(rij2 < reL2) {
					if(rij2 == 0.0) continue;
					nNeigh += 1;
				}
			}}}
			// Number of neighboors due wall
			nNeighWall.push_back(nNeigh);
		}
	}

	// Auxiliary print
	bool printData = false;
	
	if(printData){
		std::cout << " xDataPND" << std::endl;
		for(int i=0;i<xDataPND.size();i++)
			std::cout << xDataPND[i] << ", ";
		std::cout << std::endl;
		std::cout << " Zij" << std::endl;
		for(int i=0;i<Zij.size();i++)
			std::cout << Zij[i] << ", ";
		std::cout << std::endl;
		std::cout << " xDataNeigh" << std::endl;
		for(int i=0;i<xDataNeigh.size();i++)
			std::cout << xDataNeigh[i] << ", ";
		std::cout << std::endl;
		std::cout << " nNeighWall" << std::endl;
		for(int i=0;i<nNeighWall.size();i++)
			std::cout << nNeighWall[i] << ", ";
		std::cout << std::endl;
	}
}

// Weight function (Zij)
// Returns interpolated value at x from parallel arrays ( xData, Zij )
// Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
// boolean argument extrapolate determines behaviour beyond ends of array (if needed)
double PolygonMesh::interpolateWij(double re, double x, bool extrapolate)
{
	int sizeData = xDataPND.size();
	double rij_re = x/re;

	int i = 0;									// find left end of interval for interpolation
	if(rij_re >= xDataPND[sizeData - 2]) {		// special case: beyond right end
		i = sizeData - 2;
	}
	else {
		while(rij_re > xDataPND[i+1]) i++;
	}
	double xL = xDataPND[i], yL = Zij[i], xR = xDataPND[i+1], yR = Zij[i+1];	// points on either side (unless beyond ends)
	if(!extrapolate) {							// if beyond ends of array and not extrapolating
		if (rij_re < xL) yR = yL;
		if (rij_re > xR) yL = yR;
	}

	double dydx = ( yR - yL ) / ( xR - xL );	// gradient

	return yL + dydx * ( rij_re - xL );			// linear interpolation
}

// Number of neighboors due wall (numNeighWall)
// Returns interpolated value at x from parallel arrays ( xData, numNeighWall )
// Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
// boolean argument extrapolate determines behaviour beyond ends of array (if needed)
int PolygonMesh::interpolateNumNeighWall(double re, double x, bool extrapolate)
{
	int sizeData = xDataNeigh.size();
	double rij_re = x/re;

	int i = 0;									// find left end of interval for interpolation
	if(rij_re >= xDataNeigh[sizeData - 2]) {	// special case: beyond right end
		i = sizeData - 2;
	}
	else {
		while(rij_re > xDataNeigh[i+1]) i++;
	}
	double xL = xDataNeigh[i], yL = (double)nNeighWall[i], xR = xDataNeigh[i+1], yR = (double)nNeighWall[i+1];  // points on either side (unless beyond ends)
	if (!extrapolate) {							// if beyond ends of array and not extrapolating
		if (rij_re < xL) yR = yL;
		if (rij_re > xR) yL = yR;
	}

	double dydx = ( yR - yL ) / ( xR - xL );	// gradient

	return (int) (yL + dydx * ( rij_re - xL ));	// linear interpolation
}

// Find closest point on the mesh from a particle and corrects the PND and number of neighboors
// Libigl
void PolygonMesh::closestPointPNDBoundaryAABB(double reS2, double reL2, int nP, int wijType, int *Typ, int fld, int msh_id, int sta_id, 
	int fem_id, int frw_id, const double *Pos, double *wallPos, double *mirrorPos, double *riw2, int *elementID, int *meshID, double *NormalWall) {
//  double *Pos, double *wallPos, double *mirrorPos, double *riw2, double *niw, int *numNeighw, int *elementID, std::vector<int>& particlesNearMesh) {

	// MPS -> libigl
#pragma omp parallel for
	for(int i=0;i<nP;i++) {
//		if(Typ[i] == fld) {
//		temporary_position.row(i).x() = Pos[i*3  ];
//		temporary_position.row(i).y() = Pos[i*3+1];
//		temporary_position.row(i).z() = Pos[i*3+2];

			temporary_position(i,0) = Pos[i*3  ];
			temporary_position(i,1) = Pos[i*3+1];
			temporary_position(i,2) = Pos[i*3+2];
	}

	//position_transpose = position.transpose(); // Tanaka - position instead of temporary !!!

	// Find closest point Xwall_temporary_position in element "element_indice"
	// The search is faster for Static Mesh using threeMesh
	// Dynamic mesh
	if(msh_id == fem_id || msh_id == frw_id)
		igl::point_mesh_squared_distance(temporary_position,NV,NF,sqr_distance,element_indice,xWall);
	// Static mesh
	if(msh_id == sta_id)
		treeMesh.squared_distance(NV,NF,temporary_position,sqr_distance,element_indice,xWall);

	// libigl -> MPS
#pragma omp parallel for
	for(int i=0;i<nP;i++) {
//    if(Typ[i] == fld) {
		// Squared distance of particle to triangle mesh
		if(sqr_distance(i) < riw2[i]) {

			riw2[i] = sqr_distance(i);
			// Point on the mesh
			wallPos[i*3  ] = xWall.row(i).x();
			wallPos[i*3+1] = xWall.row(i).y();
			wallPos[i*3+2] = xWall.row(i).z();
			// Mirror particle position Xm = Xi + 2.0*(Xw - Xi)
			mirrorPos[i*3  ] = Pos[i*3  ] + 2.0*(wallPos[i*3  ] - Pos[i*3  ]);
			mirrorPos[i*3+1] = Pos[i*3+1] + 2.0*(wallPos[i*3+1] - Pos[i*3+1]);
			mirrorPos[i*3+2] = Pos[i*3+2] + 2.0*(wallPos[i*3+2] - Pos[i*3+2]);

			// Set if particle is close to a deformable, forced or fixed mesh
			if(msh_id == fem_id)
				meshID[i] = 1;
			else if (msh_id == frw_id)
				meshID[i] = 2;
			else
				meshID[i] = 0;

			// If is deformable mesh
			// Element ID
			if(msh_id == fem_id)
				elementID[i] = element_indice(i);
			// Polygon normal
			int eID = element_indice(i);
			//NormalWall[i*3  ] = meshNormals(eID,0);
			//NormalWall[i*3+1] = meshNormals(eID,1);
			//NormalWall[i*3+2] = meshNormals(eID,2);

			NormalWall[i*3  ] = NNormals(eID,0);
			NormalWall[i*3+1] = NNormals(eID,1);
			NormalWall[i*3+2] = NNormals(eID,2);

			//printf("i:%d riw2:%lf\n", i,riw2[i]);
			//  double dd = sqrt((pow(0.5*(Pos[i*3]-mirrorPos[i*3]),2)+pow(0.5*(Pos[i*3+1]-mirrorPos[i*3+1]),2)+pow(0.5*(Pos[i*3+2]-mirrorPos[i*3+2]),2)));
			//  printf("dist: %f %f\n",sqrt(sqr_distance(i)), dd);
		}
	}
  // Point on the mesh
//  wall_temporary_position = xWall.transpose();
  // Mirror particle position Xm = Xi + 2.0*(Xw - Xi)
//  mirror_temporary_position = temporary_position + 2.0*(wall_temporary_position - temporary_position);
  // Vector distance between i and wall
//  r_iwall = temporary_position - wall_temporary_position;
/*
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
	for(int i=0;i<nP;i++) {
	  niw[i] = 0;
	  numNeighw[i] = 0;
//      if(Typ[i] == fld) {
	// Compute wall weight function Z(xij)
	//if(i==6) {
	//  printf("x:%lf y:%lf z:%lf x:%lf y:%lf z:%lf\n", Pos[i*3],Pos[i*3+1],Pos[i*3+2],mirrorPos[i*3],mirrorPos[i*3+1],mirrorPos[i*3+2]);
	//  double dd = sqrt((pow(0.5*(Pos[i*3]-mirrorPos[i*3]),2)+pow(0.5*(Pos[i*3+1]-mirrorPos[i*3+1]),2)+pow(0.5*(Pos[i*3+2]-mirrorPos[i*3+2]),2)));
	//  printf("dist: %f %f\n",sqrt(sqr_distance(i)), dd);
	//}
//        double x2 = sqr_distance(i);
	double x2 = riw2[i];
	//printf("x: %f \n",x);
	if (x2 < reL2) {

	  double x = sqrt(x2);
	  double reS = sqrt(reS2);
	  // Add particle ID
	  particlesNearMesh_private.push_back(i);
	  //particlesNearMesh.push_back(i);

	  // PND due wall weight Z(xij)
	  niw[i] = interpolateWij(reS, x, true);
	  //printf("niw: %f\n",niw[i]);
	  // PND due wall weight Z(xij) to Gradient
	  //niw2[i] = interpolateWij(dim, 2, re, x, true);

	  double reL = sqrt(reL2);
	  // Number of neighboors due wall (numNeighWall)
	  numNeighw[i] = interpolateNumNeighWall(reL, x, true);
	  }
	}
#pragma omp critical
	  particlesNearMesh.insert(particlesNearMesh.end(), particlesNearMesh_private.begin(), particlesNearMesh_private.end());
  }
*/

#ifdef SHOW_FUNCT_NAME_POLY
	// print the function name (useful for investigating programs)
	std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif
}

// Update vector with ID of particles near the mesh
void PolygonMesh::updateParticlesNearPolygonMesh(double reS2, double reL2, int nP, int wijType, int *Typ, int fld, const double *riw2,
	double *niw, int *numNeighw, std::vector<int>& particlesNearMesh, bool *Nw) {

	particlesNearMesh.clear();

// The following code for example fills std::vectors in parallel and then combines them in the end. 
// As long as your main loop/fill function is the bottleneck this should work well in general and be thread safe.
// https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector/18671256#18671256
#pragma omp parallel
	{
		std::vector<int> particlesNearMesh_private;
#pragma omp for nowait //fill vec_private in parallel
		for(int i=0;i<nP;i++) {
			niw[i] = 0.0;
			numNeighw[i] = 0;
		    // if(Typ[i] == fld) {
			double x2 = riw2[i];
			if (x2 < reL2) {
				Nw[i]=true; // Only to show particles near polygon
				double x = sqrt(x2);
				double reS = sqrt(reS2);
				// Add particle ID
				particlesNearMesh_private.push_back(i);
				//particlesNearMesh.push_back(i);

				// PND due wall weight Z(xij)
				niw[i] = interpolateWij(reS, x, true);
				//printf("niw: %f\n",niw[i]);
				// PND due wall weight Z(xij) to Gradient
				//niw2[i] = interpolateWij(dim, 2, re, x, true);

				double reL = sqrt(reL2);
				// Number of neighboors due wall (numNeighWall)
				numNeighw[i] = interpolateNumNeighWall(reL, x, true);

				
				// if (dim == 2) {
				// 	// Mesh surface tangent of wall particle
				// 	if (x <= 0.1403*re)
				// 		numNeighw[i] += 8;
				// 	else if (x <= 0.3466*re)
				// 		numNeighw[i] += 6;
				// 	else if (x <= 0.6*re)
				// 		numNeighw[i] += 4;
				// 	else if (x <= 1.3466*re)
				// 		numNeighw[i] += 3;
				// 	else if (x <= 1.6*re)
				// 		numNeighw[i] += 1;
				// }
				// else {
				// 	// Mesh surface tangent of wall particle
				// 	if (x <= 0.0524*re)
				// 		numNeighw[i] += 22;
				// 	else if (x <= 0.1403*re)
				// 		numNeighw[i] += 18;
				// 	else if (x <= 0.3466*re)
				// 		numNeighw[i] += 14;
				// 	else if (x <= 0.6*re)
				// 		numNeighw[i] += 10;
				// 	else if (x <= 1.0524*re)
				// 		numNeighw[i] += 9;
				// 	else if (x <= 1.3466*re)
				// 		numNeighw[i] += 5;
				// 	else if (x <= 1.6*re)
				// 		numNeighw[i] += 1;
				// }
				
			}
		}
#pragma omp critical
	particlesNearMesh.insert(particlesNearMesh.end(), particlesNearMesh_private.begin(), particlesNearMesh_private.end());
	}

#ifdef SHOW_FUNCT_NAME_POLY
	// print the function name (useful for investigating programs)
	std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif
}
