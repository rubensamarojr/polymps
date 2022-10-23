// Copyright (c) 2022 Rubens AMARO
// Distributed under the MIT License.

#include <experimental/filesystem> 	///< numeric_limits
#include <iostream>					///< cout
#include "MpsVectorMatrix.h"

using namespace std;

// Constructor declaration
MpsVectorMatrix::MpsVectorMatrix()
{
}
// Destructor declaration
MpsVectorMatrix::~MpsVectorMatrix()
{
}

// Correction matrix
void MpsVectorMatrix::correctionMatrix(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsBucket *Buckets) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<Particles->numParticles; i++) {
//	if(particleType[i] == PSystem->fluid) {
		Particles->correcMatrixRow1[i*3  ] = 0.0; Particles->correcMatrixRow1[i*3+1] = 0.0; Particles->correcMatrixRow1[i*3+2] = 0.0;
		Particles->correcMatrixRow2[i*3  ] = 0.0; Particles->correcMatrixRow2[i*3+1] = 0.0; Particles->correcMatrixRow2[i*3+2] = 0.0;
		Particles->correcMatrixRow3[i*3  ] = 0.0; Particles->correcMatrixRow3[i*3+1] = 0.0; Particles->correcMatrixRow3[i*3+2] = 0.0;
		double posXi = Particles->pos[i*3  ];	double posYi = Particles->pos[i*3+1];	double posZi = Particles->pos[i*3+2];
		double posMirrorXi = Particles->mirrorParticlePos[i*3  ];	double posMirrorYi = Particles->mirrorParticlePos[i*3+1];	double posMirrorZi = Particles->mirrorParticlePos[i*3+2];
		
		int ix, iy, iz;
		Buckets->bucketCoordinates(ix, iy, iz, posXi, posYi, posZi, PSystem);
		int minZ = (iz-1)*((int)(PSystem->dim-2.0)); int maxZ = (iz+1)*((int)(PSystem->dim-2.0));
		for(int jz=minZ;jz<=maxZ;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*PSystem->numBucketsXY + jy*PSystem->numBucketsX + jx;
			int j = Particles->firstParticleInBucket[jb];
			if(j == -1) continue;
			double plx, ply, plz;
			Particles->getPeriodicLengths(jb, plx, ply, plz, PSystem);
			while(true) {
				double v0ij, v1ij, v2ij, v0imj, v1imj, v2imj, dstij2, dstimj2;
				
				// Particle square distance r_ij^2 = (Xj - Xi_temporary_position)^2
				Particles->sqrDistBetweenParticles(j, posXi, posYi, posZi, v0ij, v1ij, v2ij, dstij2, plx, ply, plz);
				// Mirror particle square distance r_imj^2 = (Xj - Xim_temporary_position)^2
				Particles->sqrDistBetweenParticles(j, posMirrorXi, posMirrorYi, posMirrorZi, v0imj, v1imj, v2imj, dstimj2, plx, ply, plz);

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < PSystem->reS2 && (dstij2 < dstimj2 || PSystem->wallType == boundaryWallType::PARTICLE)) {
				if(j != i) {
					double dst = sqrt(dstij2);
					double wS = Particles->weight(dst, PSystem->reS, PSystem->weightType);
					double invDstij2 = 1.0/dstij2;
					Particles->correcMatrixRow1[i*3  ] += wS*v0ij*v0ij*invDstij2;	Particles->correcMatrixRow1[i*3+1] += wS*v0ij*v1ij*invDstij2;	Particles->correcMatrixRow1[i*3+2] += wS*v0ij*v2ij*invDstij2;
					Particles->correcMatrixRow2[i*3  ] += wS*v1ij*v0ij*invDstij2;	Particles->correcMatrixRow2[i*3+1] += wS*v1ij*v1ij*invDstij2;	Particles->correcMatrixRow2[i*3+2] += wS*v1ij*v2ij*invDstij2;
					Particles->correcMatrixRow3[i*3  ] += wS*v2ij*v0ij*invDstij2;	Particles->correcMatrixRow3[i*3+1] += wS*v2ij*v1ij*invDstij2;	Particles->correcMatrixRow3[i*3+2] += wS*v2ij*v2ij*invDstij2;
				}}
				j = Particles->nextParticleInSameBucket[j];
				if(j == -1) break;
			}
		}}}

		// coeffPressGrad is a negative cte (-dim/noGrad)
		Particles->correcMatrixRow1[i*3  ] *= -PSystem->coeffPressGrad;	Particles->correcMatrixRow1[i*3+1] *= -PSystem->coeffPressGrad;	Particles->correcMatrixRow1[i*3+2] *= -PSystem->coeffPressGrad;
		Particles->correcMatrixRow2[i*3  ] *= -PSystem->coeffPressGrad;	Particles->correcMatrixRow2[i*3+1] *= -PSystem->coeffPressGrad;	Particles->correcMatrixRow2[i*3+2] *= -PSystem->coeffPressGrad;
		Particles->correcMatrixRow3[i*3  ] *= -PSystem->coeffPressGrad;	Particles->correcMatrixRow3[i*3+1] *= -PSystem->coeffPressGrad;	Particles->correcMatrixRow3[i*3+2] *= -PSystem->coeffPressGrad;

		// Inverse of the matrix
		int rcv = inverseMatrix(Particles->correcMatrixRow1[i*3],Particles->correcMatrixRow1[i*3+1],Particles->correcMatrixRow1[i*3+2],
								Particles->correcMatrixRow2[i*3],Particles->correcMatrixRow2[i*3+1],Particles->correcMatrixRow2[i*3+2],
								Particles->correcMatrixRow3[i*3],Particles->correcMatrixRow3[i*3+1],Particles->correcMatrixRow3[i*3+2],
								(int)PSystem->dim, PSystem->epsilonZero);

		if(Particles->numNeigh[i] < 4) {
			Particles->correcMatrixRow1[i*3  ] = 1.0;	Particles->correcMatrixRow1[i*3+1] = 0.0;	Particles->correcMatrixRow1[i*3+2] = 0.0;
			Particles->correcMatrixRow2[i*3  ] = 0.0;	Particles->correcMatrixRow2[i*3+1] = 1.0;	Particles->correcMatrixRow2[i*3+2] = 0.0;
			Particles->correcMatrixRow3[i*3  ] = 0.0;	Particles->correcMatrixRow3[i*3+1] = 0.0;	Particles->correcMatrixRow3[i*3+2] = 1.0;
		}

//		if(i == 200) {
//			printf("\n X %e %e %e ", Particles->correcMatrixRow1[i*3  ], Particles->correcMatrixRow1[i*3+1], Particles->correcMatrixRow1[i*3+2]);
//			printf("\n Y %e %e %e ", Particles->correcMatrixRow2[i*3  ], Particles->correcMatrixRow2[i*3+1], Particles->correcMatrixRow2[i*3+2]);
//			printf("\n Z %e %e %e \n", Particles->correcMatrixRow3[i*3  ], Particles->correcMatrixRow3[i*3+1], Particles->correcMatrixRow3[i*3+2]);
//		}
	}

#ifdef SHOW_FUNCT_NAME_PART
	// print the function name (useful for investigating programs)
	cout << __PRETTY_FUNCTION__ << endl;
#endif
}

// Determinant of matrix
double MpsVectorMatrix::detMatrix(const double M11, const double M12, const double M13, 
	const double M21, const double M22, const double M23, const double M31, const double M32, const double M33) {
	
	return (M11*M22*M33 + M12*M23*M31 + M13*M21*M32)
			- (M13*M22*M31 + M12*M21*M33 + M11*M23*M32);
}

// Inverse of matrix
int MpsVectorMatrix::inverseMatrix(double &M11, double &M12, double &M13, double &M21, double &M22, double &M23, 
	double &M31, double &M32, double &M33, const int dimension, const double epsZero) {
	
	double M[3][3], Maux[3][3];

	Maux[0][0] = M11;	Maux[0][1] = M12;	Maux[0][2] = M13;
	Maux[1][0] = M21;	Maux[1][1] = M22;	Maux[1][2] = M23;
	Maux[2][0] = M31;	Maux[2][1] = M32;	Maux[2][2] = M33;

	if(dimension == 2) {
		Maux[0][2] = Maux[1][2] = Maux[2][0] = Maux[2][1] = 0.0;
		Maux[2][2] = 1.0;
	}

	// Convert matrix to identity
	for(int i = 0; i < 3; i++)
	for(int j = 0; j < 3; j++) {
		if(i == j) M[i][j] = 1.0;
		else M[i][j] = 0.0;
	}

	double detM = detMatrix(Maux[0][0], Maux[0][1], Maux[0][2],
							Maux[1][0], Maux[1][1], Maux[1][2],
							Maux[2][0], Maux[2][1], Maux[2][2]);
	if(detM <= epsZero) {
		M11 = 1.0;	M12 = 0.0;	M13 = 0.0;
		M21 = 0.0;	M22 = 1.0;	M23 = 0.0;
		M31 = 0.0;	M32 = 0.0;	M33 = 1.0;
		return 0;
	}

	for(int k = 0; k < dimension; k++) {

		/*if(fabs(Maux[k][k]) <= epsZero) {
			M11 = 1.0;	M12 = 0.0;	M13 = 0.0;
			M21 = 0.0;	M22 = 1.0;	M23 = 0.0;
			M31 = 0.0;	M32 = 0.0;	M33 = 1.0;
			return 0;
		}*/

		for(int i = 0; i < dimension; i++) {
			if(i == k) continue;

			double m = Maux[i][k]/Maux[k][k];

			for(int j = 0; j < dimension; j++) {
				Maux[i][j] -= m*Maux[k][j];
				M[i][j] -= m*M[k][j];
			}
		}
	}

	for(int i = 0; i < dimension; i++)
		for(int j = 0; j < dimension; j++)
			M[i][j] /= Maux[i][i];

	M11 = M[0][0];	M12 = M[0][1];	M13 = M[0][2];
	M21 = M[1][0];	M22 = M[1][1];	M23 = M[1][2];
	M31 = M[2][0];	M32 = M[2][1];	M33 = M[2][2];

	return 1;
}

// 3D triangle to xy plane
// https://math.stackexchange.com/questions/856666/how-can-i-transform-a-3d-triangle-to-xy-plane
void MpsVectorMatrix::transformMatrix(const double *V1, const double *V2, const double *V3, double *RM, const double epsZero) {
	// double A[3];
	double B[3],C[3],U[3],V[3],W[3];
	// Translate the vertex "V1" to origin 0,0,0
	for(unsigned int i = 0; i < 3; i++) {
		//A[i] = V1[i] - V1[i];
		B[i] = V2[i] - V1[i];
		C[i] = V3[i] - V1[i];
	}
	// Define the vector U pointing from A to B, and normalise it.
	double normB = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
	if(normB > epsZero) {
		for(unsigned int i = 0; i < 3; i++)
			U[i] = B[i]/normB;
	}
	else {
		for(unsigned int i = 0; i < 3; i++)
			U[i] = 0.0;
	}
	// Take the cross product which is at right angles to the triangle.
	W[0] = U[1] * C[2] - U[2] * C[1];
	W[1] = U[2] * C[0] - U[0] * C[2];
	W[2] = U[0] * C[1] - U[1] * C[0];
	// Normalise it to give the unit vector
	double normW = sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);
	if(normW > epsZero) {
		for(unsigned int i = 0; i < 3; i++)
			W[i] = W[i]/normW;
	}
	else {
		for(unsigned int i = 0; i < 3; i++)
			W[i] = 0.0;
	}
	// Find the cross-product automatically a unit vector.
	V[0] = U[1] * W[2] - U[2] * W[1];
	V[1] = U[2] * W[0] - U[0] * W[2];
	V[2] = U[0] * W[1] - U[1] * W[0];
	// In coordinates corresponding to the basis {U, V, W }, the triangle lies. The rotation matrix that carries the usual
	// (1,0,0) to U, the usual (0,1,0) to V, and the usual (,0,0,1) to W has U, V and W as its columns.
	// You want the other direction, so take the inverse - which for a rotation matrix is just the transpose.
	for(unsigned int i = 0; i < 3; i++) {
		RM[i  ] = U[i];
		RM[i+3] = V[i];
		RM[i+6] = W[i];
	}

	/*
	// Projection of P represented in triangle coordinate system 3D (X,Y,Z) -> 2D (U,V,W).
	proj_pn = RM*proj_p';
	// Square approximation support in triangle coordinate system 3D (X,Y,Z) -> 2D (U,V,W). Here, we assume the
	// value in W as 0.
	ls = sqrt(re^2-p_pb^2);
	squareSupport = [[-ls+proj_pn(1);ls+proj_pn(1);ls+proj_pn(1);-ls+proj_pn(1)],...
	[-ls+proj_pn(2);-ls+proj_pn(2);ls+proj_pn(2);ls+proj_pn(2)]];
	// Triangle represented in triangle coordinate system 3D
	// (X,Y,Z) -> 2D (U,V,W). Here we assume the value in W as 0
	triangle2D = [[an(1);bn(1);cn(1)],[an(2);bn(2);cn(2)]];
	// Sutherland Hodgman clipping - 2D (U,V,W). Here, the value in W is assumed 0 for all points.
	clippedPolygon2D = sutherlandHodgman(triangle2D,squareSupport);
	*/
}

// 3D -> 2D
void MpsVectorMatrix::transform3Dto2D(double *P1, const double *RM) {
	double A[3];
	// Point represented in triangle coordinate system (U,V,W). The value in W is the same for all points (3D -> 2D).
	// Trasnpose of RM
	for(unsigned int i = 0; i < 3; i++) {
		A[i] = RM[3*i]*P1[0]+RM[3*i+1]*P1[1]+RM[3*i+2]*P1[2];
	}
	for(unsigned int i = 0; i < 3; i++) {
		P1[i] = A[i];
	}
}

// 2D -> 3D
void MpsVectorMatrix::transform2Dto3D(double *P1, const double *RM) {
	double A[3];
	// Point represented in coordinate system (X,Y,Z) (2D -> 3D).
	for(unsigned int i = 0; i < 3; i++) {
		A[i] = RM[i]*P1[0]+RM[i+3]*P1[1]+RM[i+6]*P1[2];
	}
	for(unsigned int i = 0; i < 3; i++) {
		P1[i] = A[i];
	}
}