#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include "polygon.h"

//Comment below to simulate only particles (fluid and wall particles)
#define POLYGON_ON

//Uncomment below to simulate coupled FEM
//#define FEM_ON

//Unomment below to use forced mesh wall (Work only if deformable mesh is the second mesh)
//#define FORCED_ON

//Uncomment below to print only free-surface and wall particles
//#define ONLY_FREE_WALL

//ASCII or BINARY output files - BINARY NOT WORKING !!!
//#define BINARY

// Mesh file
// Number of fixed meshes
#ifdef POLYGON_ON
#define NUM_POLY 1
#else
#define NUM_POLY 0
#endif

// Number of deformable meshes (Finite Element Method) - NOT WORKING !!!
#ifdef FEM_ON
#define NUM_FEM 1
#else
#define NUM_FEM 0
#endif

// Number of foced meshes
#ifdef FORCED_ON
#define NUM_FORC 1	
#else
#define NUM_FORC 0
#endif


// Number of meshes
enum {NUM_MESHS = (NUM_POLY + NUM_FEM + NUM_FORC),};

#define STA_ID 0	// Static mesh file ID
#define FEM_ID 1	// Deformable mesh file ID used for Finite Element Method - NOT WORKING !!!
#define FRW_ID 2	// Forced rigid wall file ID

// Grid file
//#define IN_FILE "input/dambreak_fluid.prof"
//#define IN_FILE "input/brumadinho_fluid_lo10p00.prof"
//#define IN_FILE "input/brumadinho_fluid_lo05p00_desl.prof"
#define IN_FILE "input/dam1610_3D_fluid_lo0p010_mps.grid"
//#define IN_FILE "input/hydro_3D_lo0p010_fluid_mps.grid"

// Static rigid wall
#define IN_MESH_0 "input/dam1610.stl"
//#define IN_MESH "input/BRUMADINHO_space10_model_lucas_top_0p50.stl"
//#define IN_MESH "input/damErosion3D.stl"
//#define IN_MESH "input/hydro_3D_box.stl"

// Deformable wall (FEM) - NOT WORKING !!!
//#define IN_MESH_1 "input/elastic_beam_16.stl" // Liao

// Forced rigid wall
//#define IN_MESH_2 "input/tank_liao_forced_gate.stl"

// Output folder
//#define OUT_FOLDER "BRUMADINHO_lo10p00_rho1500_v04p00e-01"
//#define OUT_FOLDER "BRUMADINHO_lo05p00_rho1500_v04p00e-01"
//#define OUT_FOLDER "dam1610_NONEW_02"
#define OUT_FOLDER "dambreak01"
//#define OUT_FOLDER "hydro_3D_lo0p010_03"

// Outputs (Paraview)
#define OUTPND 1
#define OUTNEIGH 0
#define OUTDEV 0
#define OUTAUX 0
#define OUTCONC 0

// Geometry
#define PCL_DST 0.01				// Average particle distance (m) (10/5) (0.01)
#define MIN_X  (0.0 - PCL_DST*3)	// Minimum value in the x direction of the analysis domain (m) (0)
#define MIN_Y  (0.0 - PCL_DST*3)	// Minimum value in the y direction of the analysis domain (m) (0)
#define MIN_Z  (0.0 - PCL_DST*3)	// Minimum value in the z direction of the analysis domain (m) (700) (0)
#define MAX_X  (1.7 + PCL_DST*3)	// Maximum value in the x direction of the analysis domain (m) (6000) (1.65) (1.7)
#define MAX_Y  (0.2 + PCL_DST*3)	// Maximum value in the y direction of the analysis domain (m) (5300) (0.15) (0.2)
#define MAX_Z  (0.7 + PCL_DST*30)	// Maximum value in the z direction of the analysis domain (m) (1300) (0.70) (0.6)

// Model
#define DIM 3				// Dimension
#define GST -1				// Ghost particle ID
#define FLD 0				// Fluid particle ID
#define WLL 1				// Wal particle ID
#define SRF 1				// Free-surface particle BC
#define INR 0				// Inner particle BC
#define OTR -1				// Other particle BC
#define NUM_TYP 2			// Number of particle types	
// Physical
#define DNS_FLD 1000		// Fluid particle density (kg/m3)
#define DNS_WLL 1000		// Wall particle density (kg/m3)
#define KNM_VS1 0.000001			// Kinematic viscosity phase 1 (m2/s)
#define KNM_VS2 0.000001		// Kinematic viscosity phase 2 (m2/s)
#define DNS_FL1 1000.0		// Fluid particle density phase 1 (kg/m3)
#define DNS_FL2 1540.0		// bulk particle density phase 2 (kg/m3)
#define DNS_SDT	1540.0		// sediment density (kg/m3)
#define G_X 0.0				// Gravity acceleration x component (m/s2)
#if DIM==2
#define G_Y -9.81			// Gravity acceleration y component (m/s2)
#define G_Z  0.00			// Gravity acceleration z component (m/s2)
#else
#define G_Y  0.00			// Gravity acceleration y component (m/s2)
#define G_Z -9.81			// Gravity acceleration z component (m/s2)
#endif
// Rheological parameters
#define Fluid2_type 0		// Newtonian:0 , Non Newtonian:1 
#define N 1.2				// flow behaviour (power law) index
#define MEU0 0.03			// consistency index
#define PHI 0.541			// friction angle (RAD) lower limit
#define PHI_WAL 0.541		// friction angle (RAD)
#define PHI_BED 0.541		// friction angle (RAD)
#define PHI_2 0.6			// second friction angle Values are based on  Minatti & Paris (2015) upper limit
#define cohes 0.0				// cohesiveness coefficient
#define Fraction_method 2   // Method of calculation of volume of fraction. 1: Linear dist across the interface, 2: smoothed value
//#define visc_max 20			// maximum viscosity uses to avoid singularity
#define DG 0.0035			// grain size
#define I0 0.75				// I0 value in Meu9I0 rheology     Values are based on  Minatti &  Paris (2015)
#define mm 100.0			// Regularization parameter. Control the viscosity grows
#define stress_cal_method 1	// Method 1; viscosity is directly used in momentum equation. Method 2: first the stress tensor is calculated then it is used in momentum equation
#define visc_itr_num 1
#define visc_error 0.0
#define visc_ave 0.0
#define Cd 0.47			// Drag coefficient
#define VF_min 0.25		// Minimum volume fraction
#define VF_max 0.65		// Maximum volume fraction

// Numerical
#define DT 0.00025				// Time step (s) (0.02/0.01) (0.00025)
#define FIN_TIM 1.00		// Time of simulation (s)
#define OPT_FQC 100			// Number of iterations to determine the output interval
#define reS2D	3.1			// Influence radius small 2D
#define reL2D	3.1			// Influence radius large 2D
#define reS3D	2.1			// Influence radius small 3D
#define reL3D	2.1			// Influence radius large 3D
#define MPS_TYP	1			// Explicit MPS = 0;  Weakly compressible MPS = 1
#define WGT_TYP	0			// Weight function re/rij - 1 = 0; re/rij + rij/re - 2 = 1; re/rij - rij/re = 2; pow(1-rij/re,3.0) = 3;
#define GRD_TYP	3			// Pressure gradient: Pj - Pmin = 0; Pj + Pi = 1; Pj + Pi - 2*Pmin = 2; ni*Pj/nj + nj*Pi/ni = 3
#define GRD_COR	0			// Corrected pressure gradient. No = 0; Yes = 1
#define RLX_PRS 1.0			// Relaxation factor for pressure correction (<= 1)
#define SND 15.00			// Sound speed (m/s) (10)
#define GAM 7.0				// Gamma weakly compressible MPS
#define ADJ_VEL 2			// Adjusted velocity: No = 0; DW*Uij = 1; GradCi = 2
#define DRI 0.01			// Adjusted velocity paramater (DW*Uij) DRI <= 0.01
#define COEF_A 1.0			// Dimensionless number 1-6 (GradCi) COEF_A = 2.0 provides a good compromise
#define MACH 0.1			// Mach number 0.1 (GradCi)
#define VEL_A 0.9			// Adjusted velocity paramater (3) a = 0.9 (NOT IMPLEMENTED !!!)
#define PND_CAL	2			// PND = Soma(wij) = 0; PND = soma(PNDj)*wij/soma(wij) = 1; PND = Diffusion term = 2
#define DFS_MPS 0.35		// Diffusion term
#define EPS_RE 0.01			// 
#define REP_FOR 2			// Wall repulsive force: Harada = 0, Mitsume = 1; Lennard-Jones = 2; Monaghan-Kajtar = 3
#define REP_RE	0.5			// Influence radius for repulsive force
#define MIT_REP 40000000.0		// Wall coefficent repulsive force (10000/100000) (500000) (Mitusme)
#define LNJ_REP 2.0			// Wall coefficent repulsive force (1-10) (Lennard-Jones)
#define KJT_REP 1.0			// Wall coefficent repulsive force (1-10) (Monaghan-Kajtar)
#define SLP 0				// No-slip = 0 ; Free-slip = 1 
#define CRT_NUM 0.2			// Courant (CFL) condition number
#define PND_TRS 0.98		// Surface threshold PND (2D - 0.97/ 3D - 0.93)
#define NGH_TRS 0.85		// Surface threshold Neighboors
#define DLT_TRS 0.2			// NPCD threshold
#define COL_RAT 0.2			// Collision ratio
#define DST_LMT_RAT 0.85	// Coefficient of distance which does not allow any further access between particles (0.9)

#if WGT_TYP==0
#define WEI(dst, re) ((re/dst) - 1.0)	// Weight function
#define WEI_GRAD(dst, re) ((re/dst) - 1.0)	// Weight function Gradient
#define DEL_WEI(dst, re) (-re/(dst*dst))	// Derivate of weight function
#elif WGT_TYP==1
#define WEI(dst, re) (re/dst + dst/re - 2)	// Weight function
#define WEI_GRAD(dst, re) (re/dst + dst/re - 2)	// Weight function Gradient
#define DEL_WEI(dst, re) (-re/(dst*dst) + 1/re)	// Derivate of weight function
#elif WGT_TYP==2
#define WEI(dst, re) (re/dst - dst/re)	// Weight function
#define WEI_GRAD(dst, re) (re/dst - dst/re)	// Weight function Gradient
#define DEL_WEI(dst, re) (-re/(dst*dst) - 1/re)	// Derivate of weight function
#elif WGT_TYP==3
#define WEI(dst, re) ((1-dst/re)*(1-dst/re)*(1-dst/re))	// Weight function
#define WEI_GRAD(dst, re) ((1-dst/re)*(1-dst/re)*(1-dst/re))	// Weight function Gradient
#define DEL_WEI(dst, re) (-3/re*(1-(dst/re))*(1-(dst/re)))	// Derivate of weight function
#endif
//#define WEI(dst, re) ((re/dst+dst/re) - 2.0)	// Weight function
//#define WEI(dst, re) (re/dst-dst/re)	// Weight function
//#define WEI_GRAD(dst, re) (re/dst-dst/re)	// Weight function Gradient
//#define WEI_WEND(dst, re) ((1-dst/re)*(1-dst/re)*(1-dst/re)*(1-dst/re)*(1+4*dst/re))	// Wendland weight function
//#define WEI_WEND(dst, re) ((1-dst*dst/(re*re)))
//#define WEI_WEND(dst, re) ((1-dst/re)*(1-dst/re)*(1-dst/re))	// Weight function


FILE* fp;
char filename[256];
int iLP,iF;
double TIM, vMax, Crt, timer_sta, timer_end;
int nP;
double *Acc,*Pos,*Vel,*Prs,*pav,*pndi,*dev,*pndS,*devMod2,*Ci,*gradCi,*CMr1,*CMr2,*CMr3,*DIVi,*DIi;
int *Typ,*numNeigh,*Bc;
double reS,reS2,reL,reL2,eps_reS,re_rep;
double DB,DB2,DBinv;
int nBx,nBy,nBz,nBxy,nBxyz;
int *bfst,*blst,*nxt;
double n0S,n0L,lmd,A1,A2,A3,A4,A5,A6,rlim,rlim2,COL,nNeigh0,n0Grad;
double Dns[NUM_TYP],invDns[NUM_TYP];
// Polygon
double *wallPos,*mirrorPos,*niw,*riw2;
int *Nw,*numNeighw,*meshID;
double *F1,*F2,*Normal,*AccStar,*NormalWall,*devXnormal,*numNeighFree;
int wijType;
//double *Posk, *Velk, *Acv;

// Non-Newtonian
double A1_M, NEU;
double *Cv, *II, *MEU, *MEU_Y, *Inertia, *pnew, *p_rheo_new, *RHO, *p_smooth, *VF;
double *S12, *S13, *S23, *S11, *S22, *S33;
int *PTYPE;

// Vector with ID of particles near to mesh
std::vector<int> partNearMesh;

// Auxiliary mesh
int *elementID; // Element indice

// Forced rigid wall
double *nodeFRWX, *nodeFRWY, *nodeFRWZ; // Nodal positions
double velVWall[3]; // Uniform wall velocity 

// Return time
double get_dtime(void) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 0.000001);
}

// Verify if particle is out of domain
void ChkPcl(int i) {
	if(	Pos[i*3  ]>MAX_X || Pos[i*3  ]<MIN_X ||
		Pos[i*3+1]>MAX_Y || Pos[i*3+1]<MIN_Y ||
		Pos[i*3+2]>MAX_Z || Pos[i*3+2]<MIN_Z) {
		Typ[i] = GST; Bc[i] = OTR;

		Acc[i*3]=Acc[i*3+1]=Acc[i*3+2]=0.0;
		AccStar[i*3]=AccStar[i*3+1]=AccStar[i*3+2]=0.0;
		Vel[i*3]=Vel[i*3+1]=Vel[i*3+2]=0.0;
		dev[i*3]=dev[i*3+1]=dev[i*3+2]=0.0;
		Prs[i]=pav[i]=pndi[i]=pndS[i]=devMod2[i]=Ci[i]= 0.0;
		gradCi[i*3]=gradCi[i*3+1]=gradCi[i*3+2]=0.0;
		CMr1[i*3]=CMr1[i*3+1]=CMr1[i*3+2]=0.0;
		CMr2[i*3]=CMr2[i*3+1]=CMr2[i*3+2]=0.0;
		CMr3[i*3]=CMr3[i*3+1]=CMr3[i*3+2]=0.0;
		numNeigh[i]=Nw[i]=numNeighw[i]=meshID[i]=0;numNeighFree[i]=0.0;
		niw[i]=DIVi[i]=DIi[i]=0.0;
		F1[i*3]=F1[i*3+1]=F1[i*3+2]=0.0;
		F2[i*3]=F2[i*3+1]=F2[i*3+2]=0.0;
	}
}

// Read input data
void RdDat(void) {
	fp = fopen(IN_FILE, "r");
	if(fp == NULL) perror ("Error opening grid file");

	int zeroZero;
	fscanf(fp,"%d",&zeroZero);
	fscanf(fp,"%d",&nP);
	printf("nP: %d\n",nP);
	Acc = (double*)malloc(sizeof(double)*nP*3);			// Particle acceleration
	Pos = (double*)malloc(sizeof(double)*nP*3);			// Particle coordinates
	Vel = (double*)malloc(sizeof(double)*nP*3);			// Particle velocity
	Prs = (double*)malloc(sizeof(double)*nP);			// Particle pressure
	pav = (double*)malloc(sizeof(double)*nP);			// Time averaged particle pressure
	Typ = (int*)malloc(sizeof(int)*nP);					// Particle type
	Normal = (double*)malloc(sizeof(double)*nP*3);		// Particle normal
	dev = (double*)malloc(sizeof(double)*nP*3);			// NPCD deviation
	devMod2 = (double*)malloc(sizeof(double)*nP);		// NPCD deviation modulus
	Ci = (double*)malloc(sizeof(double)*nP);			// Concentration
	gradCi = (double*)malloc(sizeof(double)*nP*3);		// Gradient of concentration
	CMr1 = (double*)malloc(sizeof(double)*nP*3);		// Correction matrix - Row 1
	CMr2 = (double*)malloc(sizeof(double)*nP*3);		// Correction matrix - Row 2
	CMr3 = (double*)malloc(sizeof(double)*nP*3);		// Correction matrix - Row 3
	numNeigh = (int*)malloc(sizeof(int)*nP);			// Number of neighboors
	pndi = (double*)malloc(sizeof(double)*nP);			// PND
	pndS = (double*)malloc(sizeof(double)*nP);			// PND small = sum(wij)
	Bc = (int*)malloc(sizeof(int)*nP);					// BC particle type
	AccStar = (double*)malloc(sizeof(double)*nP*3);		// Particle acceleration due gravity and viscosity

	DIVi = (double*)malloc(sizeof(double)*nP);			// Divergence of velocity
	DIi = (double*)malloc(sizeof(double)*nP);			// Diffusive term

	// Polygons
	wallPos = (double*)malloc(sizeof(double)*nP*3);		// Particle at wall coordinate
	mirrorPos = (double*)malloc(sizeof(double)*nP*3);	// Mirrored particle coordinate
	niw = (double*)malloc(sizeof(double)*nP);			// PND wall
	Nw = (int*)malloc(sizeof(int)*nP);					// Particle near wall
	numNeighw = (int*)malloc(sizeof(int)*nP);			// Number of neighbors due wall
	F1 = (double*)malloc(sizeof(double)*nP*3);			// Wall-Particle force
	F2 = (double*)malloc(sizeof(double)*nP*3);			// Wall-Particle force
	riw2 = (double*)malloc(sizeof(double)*nP);			// Squared distance of particle to triangle mesh
	meshID = (int*)malloc(sizeof(int)*nP);				// Particle near deformable mesh
	NormalWall = (double*)malloc(sizeof(double)*nP*3);	// Polygon normal
	devXnormal = (double*)malloc(sizeof(double)*nP);	// Deviation vector X polygonal wall
	numNeighFree = (double*)malloc(sizeof(double)*nP);	// Number of free-surface particle neighbors

//	Posk = (double*)malloc(sizeof(double)*nP*3);	// Particle coordinates
//	Velk = (double*)malloc(sizeof(double)*nP*3);	// Particle velocity
//	Acv = (double*)malloc(sizeof(double)*nP*3);	// Part

	// Non-Newtonian
	Cv = (double*)malloc(sizeof(double)*nP);			// Concentration
	II = (double*)malloc(sizeof(double)*nP);			// Invariant
	PTYPE = (int*)malloc(sizeof(int)*nP);				// Type of fluid
	MEU = (double*)malloc(sizeof(double)*nP);			// Dynamic viscosity
	MEU_Y = (double*)malloc(sizeof(double)*nP);			// Dynamic viscosity ??
	Inertia = (double*)malloc(sizeof(double)*nP);		//
	pnew = (double*)malloc(sizeof(double)*nP);			// New pressure
	p_rheo_new = (double*)malloc(sizeof(double)*nP);	//
	RHO = (double*)malloc(sizeof(double)*nP);			// Fluid density
	p_smooth = (double*)malloc(sizeof(double)*nP);		//
	VF = (double*)malloc(sizeof(double)*nP);			//
	S12 = (double*)malloc(sizeof(double)*nP);			//
	S13 = (double*)malloc(sizeof(double)*nP);			//
	S23 = (double*)malloc(sizeof(double)*nP);			//
	S11 = (double*)malloc(sizeof(double)*nP);			//
	S22 = (double*)malloc(sizeof(double)*nP);			//
	S33 = (double*)malloc(sizeof(double)*nP);			//


	elementID = (int*)malloc(sizeof(int)*nP);			// Element ID

	for(int i=0; i<nP; i++) {
		int a[2];
		double b[8];

		// Uncomment here to read .prof file
		//fscanf(fp," %d %d %lf %lf %lf %lf %lf %lf %lf %lf",&a[0],&a[1],&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7]);
		// Uncomment here to read .grid file
		a[0] = 0;
		fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf",&a[1],&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7]);
		Typ[i]=a[1];
		Pos[i*3]=b[0];	Pos[i*3+1]=b[1];	Pos[i*3+2]=b[2];
		Vel[i*3]=b[3];	Vel[i*3+1]=b[4];	Vel[i*3+2]=b[5];
		Prs[i]=b[6];	pav[i]=b[7];

//		printf("X: %d %lf %lf %lf %lf %lf %lf %lf %lf\n",Typ[i],Pos[3*i],Pos[3*i+1],Pos[3*i+2],Vel[3*i],Vel[3*i+1],Vel[3*i+2],Prs[i],pav[i]);
//		Posk[i*3]=b[0];	Posk[i*3+1]=b[1];	Posk[i*3+2]=b[2];
//		Velk[i*3]=b[3];	Velk[i*3+1]=b[4];	Velk[i*3+2]=b[5];
	}
	fclose(fp);
	for(int i=0; i<nP; i++) {ChkPcl(i);}
	for(int i=0;i<nP*3;i++) {
		Acc[i]=0.0;Normal[i]=0.0;AccStar[i]=0.0;NormalWall[i]=0.0;CMr1[i]=0.0;CMr2[i]=0.0;CMr3[i]=0.0;/*Acv[i]=0.0;*/
		wallPos[i]=0.0;mirrorPos[i]=0.0;F1[i]=0.0;F2[i]=0.0;dev[i]=0.0;gradCi[i]=0.0;
	}
	for(int i=0; i<nP; i++) {
		niw[i]=0.0;numNeigh[i]=0;pndi[i]=0.0;Bc[i]=0;Nw[i]=0;numNeighw[i]=0;pndS[i]=0.0;devMod2[i]=0.0;Ci[i]=0.0;Cv[i]=0.0;
		II[i]=0.0;PTYPE[i]=0;MEU[i]=0.0;MEU_Y[i]=0.0;Inertia[i]=0.0;pnew[i]=0.0;p_rheo_new[i]=0.0;p_smooth[i]=0.0;VF[i]=0.0;
		S12[i]=0.0;S13[i]=0.0;S23[i]=0.0;S11[i]=0.0;S22[i]=0.0;S33[i]=0.0;elementID[i]=0;DIVi[i]=0.0;DIi[i]=0.0;
		meshID[i]=0;devXnormal[i]=0.0;numNeighFree[i]=0.0;
		riw2[i]=10e8*PCL_DST;

		/*
		// Assign type and density
		if(Pos[i*3+2] <= 0.3) {
			PTYPE[i]=2;
			RHO[i] = DNS_FL2;
			// CHANGED Only for the first time step
			MEU[i] = KNM_VS2 * DNS_FL2;
		}
		else {
			PTYPE[i]=1;
			RHO[i] = DNS_FL1;
			// CHANGED Only for the first time step
			MEU[i] = KNM_VS1 * DNS_FL1;
		}
		*/
		// Assign type and density
		if(Typ[i]==1) {
			Typ[i]=0;
			PTYPE[i]=2;
			RHO[i] = DNS_FL2;
			// CHANGED Only for the first time step
			MEU[i] = KNM_VS2 * DNS_FL2;
		}
		else {
			Typ[i]=0;
			PTYPE[i]=1;
			RHO[i] = DNS_FL1;
			// CHANGED Only for the first time step
			MEU[i] = KNM_VS1 * DNS_FL1;
		}
	}
}
/*
void RdMes(void) {
	fp = fopen(IN_FILE, "r");
	fscanf(fp,"%d",&nP);
	mesh mesh(IN_FILE, nP);
	mesh.readMeshFile(IN_MESH_0);
}
*/

// Write data. Format .prof
void WrtDat(void) {
	char outout_filename[256];
	sprintf(outout_filename, "output%05d.prof",iF);
	fp = fopen(outout_filename, "w");
	fprintf(fp,"%d\n",nP);
	for(int i=0; i<nP; i++) {
		int a[2];
		double b[9];
		a[0]=i;	a[1]=Typ[i];
		b[0]=Pos[i*3];	b[1]=Pos[i*3+1];	b[2]=Pos[i*3+2];
		b[3]=Vel[i*3];	b[4]=Vel[i*3+1];	b[5]=Vel[i*3+2];
		b[6]=Prs[i];		b[7]=pav[i]/OPT_FQC;
		b[8]=pndi[i];
		fprintf(fp," %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",a[0],a[1],b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[8]);
		pav[i]=0.0;
	}
	fclose(fp);
	iF++;
}

// https://stackoverflow.com/questions/105252
// https://stackoverflow.com/questions/10913666/error-writing-binary-vtk-files
// https://stackoverflow.com/questions/55829282/write-vtk-file-in-binary-format
template <typename T>
void SwapEnd(T& var)
{
	char* varArray = reinterpret_cast<char*>(&var);
	for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
		std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}

// Write data. Format .vtu (Paraview) - NOT WORKING !!!
void WrtVtuBinary(void) {

	char outout_filename[256];
	sprintf(outout_filename, OUT_FOLDER"/output%05d.vtu",iF);

	int nParticles = 0;
#ifdef ONLY_FREE_WALL
	for(int i=0; i<nP; i++){
		if(Bc[i] == SRF || Typ[i] == WLL)
			nParticles++;
	}
#else
	nParticles = nP;
#endif

	// BINARY FILE
	std::ofstream file;
	file.open(outout_filename, std::ios::out | std::ios::binary);

	// Header
	//file << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
	file << "<VTKFile type='UnstructuredGrid' version='1.0' byte_order='LittleEndian' header_type='UInt64'>" << std::endl;
	file << "  <UnstructuredGrid>" << std::endl;
	file << "    <Piece NumberOfPoints='" << nParticles << "' NumberOfCells='" << nParticles << "'>" << std::endl;
	
	// Point data
	file << "      <PointData>" << std::endl;
	file << "        <DataArray type='Float32' Name='Velocity' NumberOfComponents='3' format='binary'>" << std::endl;
	for(size_t i=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		{
			float ptx = (float)Vel[i*3  ];
			float pty = (float)Vel[i*3+1];
			float ptz = (float)Vel[i*3+2];

			SwapEnd(ptx);
			SwapEnd(pty);
			SwapEnd(ptz);

			file.write(reinterpret_cast<char*>(&ptx), sizeof(float));
			file.write(reinterpret_cast<char*>(&pty), sizeof(float));
			file.write(reinterpret_cast<char*>(&ptz), sizeof(float));
		}
	}
	file << std::endl << "        </DataArray>" << std::endl;
	
	file << "        <DataArray type='Float32' Name='Pressure' format='binary'>" << std::endl;
	for(size_t i=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		{
			float ptv = (float)Prs[i];
			SwapEnd(ptv);
			file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
		}
	}
	file << std::endl << "        </DataArray>" << std::endl;
	
	file << "<        DataArray type='Int32' Name='BC' format='binary'>" << std::endl;
	for(size_t i=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		{
			int32_t type = Bc[i];
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
		}
	}
	file << std::endl << "        </DataArray>" << std::endl;

	if(OUTPND) {
		file << "        <DataArray type='Float32' Name='pnd' format='binary'>" << std::endl;
		for(size_t i=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
			if(Bc[i] == SRF || Typ[i] == WLL)
#endif
			{
				float ptv = pndi[i];
				SwapEnd(ptv);
				file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
			}
		}
		file << std::endl << "        </DataArray>" << std::endl;
		
		file << "        <DataArray type='Float32' Name='pndS' format='binary'>" << std::endl;
		for(size_t i=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
			if(Bc[i] == SRF || Typ[i] == WLL)
#endif
			{
				float ptv = pndS[i];
				SwapEnd(ptv);
				file.write(reinterpret_cast<char*>(&ptv), sizeof(float));
			}
		}
		file << std::endl << "        </DataArray>" << std::endl;
		
	}

	if(OUTNEIGH) {
		file << "        <DataArray type='Int32' Name='nNeigh' format='binary'>" << std::endl;
		for(size_t i=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
			if(Bc[i] == SRF || Typ[i] == WLL)
#endif
			{
				int32_t type = numNeigh[i];
				int type_i = static_cast<int>(type);
				SwapEnd(type_i);
				file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
			}
		}
		file << std::endl << "        </DataArray>" << std::endl;
		
	}

	file << "      </PointData>" << std::endl;

	// Cell data
	file << "      <CellData>" << std::endl;
	file << "      </CellData>" << std::endl;

	//Points
	file << "      <Points>" << std::endl;
	file << "        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='binary'>" << std::endl;
	for(size_t i=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		{
			float ptx = (float)Pos[i*3  ];
			float pty = (float)Pos[i*3+1];
			float ptz = (float)Pos[i*3+2];

			SwapEnd(ptx);
			SwapEnd(pty);
			SwapEnd(ptz);

			file.write(reinterpret_cast<char*>(&ptx), sizeof(float));
			file.write(reinterpret_cast<char*>(&pty), sizeof(float));
			file.write(reinterpret_cast<char*>(&ptz), sizeof(float));
		}
	}
	file << std::endl <<"        </DataArray>" << std::endl;
	file << "      </Points>" << std::endl;

	// Cells
	file << "      <Cells>" << std::endl;
	file << "        <DataArray type='Int64' Name='connectivity' format='binary'>" << std::endl;
	for(size_t i=0, ii=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		{
			int64_t type = ii;
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
			ii++;
		}
	}
	file << std::endl << "        </DataArray>" << std::endl;
	
	file << "        <DataArray type='Int64' Name='offsets' format='binary'>" << std::endl;
	for(size_t i=0, ii=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		{
			int64_t type = ii+1;
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
			ii++;
		}
	}
	file << std::endl << "        </DataArray>" << std::endl;

	file << "        <DataArray type='UInt8' Name='types' format='binary'>" << std::endl;
	for(size_t i=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		{
			uint8_t type = 1;
			int type_i = static_cast<int>(type);
			SwapEnd(type_i);
			file.write(reinterpret_cast<char*>(&type_i), sizeof(int));
		}
	}
	file << std::endl << "        </DataArray>" << std::endl;
	file << "      </Cells>" << std::endl;
	file << "    </Piece>" << std::endl;
	file << "  </UnstructuredGrid>" << std::endl;
	file << "</VTKFile>" << std::endl;

	file.close();

	iF++;
}


// Write data. Format .vtu (Paraview)
void WrtVtuAscii(void) {

	char outout_filename[256];
	sprintf(outout_filename, OUT_FOLDER"/output%05d.vtu",iF);

	int nParticles = 0;
#ifdef ONLY_FREE_WALL
	for(int i=0; i<nP; i++){
		if(Bc[i] == SRF || Typ[i] == WLL)
			nParticles++;
	}
#else
	nParticles = nP;
#endif

	// ASCII FILE
	fp = fopen(outout_filename, "w");

	// Header
	fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
	fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
	fprintf(fp,"  <UnstructuredGrid>\n");
 	fprintf(fp,"    <Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",nParticles,nParticles);

 	// Points
 	fprintf(fp,"      <Points>\n");
 	fprintf(fp,"        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='ascii'>\n          ");
	for(int i=0; i<nP; i++){
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		{
			fprintf(fp,"%lf %lf %lf ",Pos[i*3],Pos[i*3+1],Pos[i*3+2]);
		}
	}
 	fprintf(fp,"\n        </DataArray>\n");
  	fprintf(fp,"      </Points>\n");

  	// Point data
  	fprintf(fp,"      <PointData>\n");
 	fprintf(fp,"        <DataArray type='Float32' Name='Velocity' NumberOfComponents='3' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		//double val=sqrt(Vel[i*3]*Vel[i*3]+Vel[i*3+1]*Vel[i*3+1]+Vel[i*3+2]*Vel[i*3+2]);
		fprintf(fp,"%f %f %f ",(float)Vel[i*3],(float)Vel[i*3+1],(float)Vel[i*3+2]);
	}
 	fprintf(fp,"\n        </DataArray>\n");
 	//fprintf(fp,"        <DataArray type='Float32' Name='pressave' format='ascii'>\n");
	//for(int i=0; i<nP; i++) {	fprintf(fp,"%f ",(float)(pav[i]/OPT_FQC));}
	//fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='pressure' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		fprintf(fp,"%f ",(float)Prs[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='BC' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		fprintf(fp,"%d ",Bc[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	if(OUTPND) {
 		fprintf(fp,"        <DataArray type='Float32' Name='pnd' format='ascii'>\n");
		for(int i=0; i<nP; i++) {	
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
			fprintf(fp,"%f ",(float)pndi[i]);
		}
 		fprintf(fp,"\n        </DataArray>\n");
 		fprintf(fp,"        <DataArray type='Float32' Name='pndS' format='ascii'>\n");
		for(int i=0; i<nP; i++) {	
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
			fprintf(fp,"%f ",(float)pndS[i]);
		}
 		fprintf(fp,"\n        </DataArray>\n");
 	}
 	if(OUTNEIGH) {
 		fprintf(fp,"        <DataArray type='Int32' Name='nNeigh' format='ascii'>\n");
		for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
			fprintf(fp,"%d ",numNeigh[i]);
		}
 		fprintf(fp,"\n        </DataArray>\n");
 	}
 	if(OUTDEV) {
//		fprintf(fp,"        <DataArray type='Float32' Name='devSquare' format='ascii'>\n");
//		for(int i=0; i<nP; i++) {	fprintf(fp,"%f ",(float)devMod2[i]);}
// 		fprintf(fp,"\n        </DataArray>\n");
 		fprintf(fp,"        <DataArray type='Float32' Name='dev' NumberOfComponents='3' format='ascii'>\n");
 		for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
			fprintf(fp,"%f %f %f ",(float)dev[i*3],(float)dev[i*3+1],(float)dev[i*3+2]);
		}
 		fprintf(fp,"\n        </DataArray>\n");
 		fprintf(fp,"        <DataArray type='Float32' Name='NormalWall' NumberOfComponents='3' format='ascii'>\n");
 		for(int i=0; i<nP; i++) {
		fprintf(fp,"%f %f %f ",(float)NormalWall[i*3],(float)NormalWall[i*3+1],(float)NormalWall[i*3+2]);
		}
 		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='devXnormal' format='ascii'>\n");
//		for(int i=0; i<nP; i++) {	fprintf(fp,"%f ",(float)devXnormal[i]);}
// 		fprintf(fp,"\n        </DataArray>\n");
 	}
 	if(OUTCONC) {
 		fprintf(fp,"        <DataArray type='Float32' Name='Ci' format='ascii'>\n");
		for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
			fprintf(fp,"%f ",(float)Ci[i]);
		}
 		fprintf(fp,"\n        </DataArray>\n");
 		fprintf(fp,"        <DataArray type='Float32' Name='GradCi' NumberOfComponents='3' format='ascii'>\n");
		for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
			fprintf(fp,"%f %f %f ",(float)gradCi[i*3],(float)gradCi[i*3+1],(float)gradCi[i*3+2]);
		}
 		fprintf(fp,"\n        </DataArray>\n");
 	}
 	
//	fprintf(fp,"        <DataArray type='Float32' Name='F1' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<nP; i++) {
//		fprintf(fp,"%f %f %f ",(float)F1[i*3],(float)F1[i*3+1],(float)F1[i*3+2]);
//	}
 //	fprintf(fp,"\n        </DataArray>\n");

// 	fprintf(fp,"        <DataArray type='Float32' Name='F2' NumberOfComponents='3' format='ascii'>\n");
//	for(int i=0; i<nP; i++) {
//		fprintf(fp,"%f %f %f ",(float)F2[i*3],(float)F2[i*3+1],(float)F2[i*3+2]);
//	}
 //	fprintf(fp,"\n        </DataArray>\n");
 	
	fprintf(fp,"        <DataArray type='Float32' Name='RHO' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
	if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		fprintf(fp,"%f ",(float)RHO[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='ConcVol' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
	if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		fprintf(fp,"%f ",(float)Cv[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='MEU' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
	if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		fprintf(fp,"%f ",(float)MEU[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='MEUy' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
	if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		fprintf(fp,"%f ",(float)MEU_Y[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='II' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
	if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		fprintf(fp,"%f ",(float)II[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Int32' Name='PTYPE' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
	if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		fprintf(fp,"%d ",PTYPE[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	fprintf(fp,"        <DataArray type='Float32' Name='p_smooth' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
	if(Bc[i] == SRF || Typ[i] == WLL)
#endif
		fprintf(fp,"%f ",(float)p_smooth[i]);
	}
	fprintf(fp,"\n        </DataArray>\n");
	
  	if(OUTAUX) {

// 		fprintf(fp,"        <DataArray type='Int32' Name='ParticleType' format='ascii'>\n");
//		for(int i=0; i<nP; i++) {
//#ifdef ONLY_FREE_WALL
//		if(Bc[i] == SRF || Typ[i] == WLL)
//#endif
//			fprintf(fp,"%d ",Typ[i]);
//		}
 //		fprintf(fp,"\n        </DataArray>\n");

// 		fprintf(fp,"        <DataArray type='Float32' Name='mirrorPos' NumberOfComponents='3' format='ascii'>\n");
//		for(int i=0; i<nP; i++) {
//#ifdef ONLY_FREE_WALL
//		if(Bc[i] == SRF || Typ[i] == WLL)
//#endif
//			fprintf(fp,"%f %f %f ",(float)mirrorPos[i*3],(float)mirrorPos[i*3+1],(float)mirrorPos[i*3+2]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");

// 		fprintf(fp,"        <DataArray type='Int32' Name='NearWall' format='ascii'>\n");
//		for(int i=0; i<nP; i++) {
//#ifdef ONLY_FREE_WALL
//		if(Bc[i] == SRF || Typ[i] == WLL)
//#endif
//			fprintf(fp,"%d ",Nw[i]);
//		}
// 		fprintf(fp,"\n        </DataArray>\n");

 		//fprintf(fp,"        <DataArray type='Float32' Name='Fwall' NumberOfComponents='3' format='ascii'>\n");
		//for(int i=0; i<nP; i++) {
		//	fprintf(fp,"%f %f %f ",(float)Fwall[i*3],(float)Fwall[i*3+1],(float)Fwall[i*3+2]);
		//}
 		//fprintf(fp,"\n        </DataArray>\n");
// 		fprintf(fp,"        <DataArray type='Float32' Name='Fwall' NumberOfComponents='3' format='ascii'>\n");
//		for(int i=0; i<nP; i++) {
//			fprintf(fp,"%f %f %f ",(float)forceWall[i*3],(float)forceWall[i*3+1],(float)forceWall[i*3+2]);
//		}
//		fprintf(fp,"\n        </DataArray>\n");
// 		fprintf(fp,"        <DataArray type='Float32' Name='DIV' format='ascii'>\n");
//		for(int i=0; i<nP; i++) {	fprintf(fp,"%f ",(float)DIVi[i]);}
// 		fprintf(fp,"\n        </DataArray>\n");
//		fprintf(fp,"        <DataArray type='Float32' Name='Di' format='ascii'>\n");
//		for(int i=0; i<nP; i++) {	fprintf(fp,"%f ",(float)DIi[i]);}
// 		fprintf(fp,"\n        </DataArray>\n");

 		fprintf(fp,"        <DataArray type='Int32' Name='meshID' format='ascii'>\n");
		for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#endif
			fprintf(fp,"%d ",meshID[i]);
		}
 		fprintf(fp,"\n        </DataArray>\n");

 		//fprintf(fp,"        <DataArray type='Float32' Name='numNeighFree' format='ascii'>\n");
		//for(int i=0; i<nP; i++) {fprintf(fp,"%f ",(float)numNeighFree[i]);}
 		//fprintf(fp,"\n        </DataArray>\n");
 	
 	}

  	fprintf(fp,"      </PointData>\n");

  	// Cells
 	fprintf(fp,"      <Cells>\n");
 	fprintf(fp,"        <DataArray type='Int32' Name='connectivity' format='ascii'>\n");
	for(int i=0, ii = 0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#else
		{
			fprintf(fp,"%d ",ii);
			ii++;
		}
#endif
	}
 	fprintf(fp,"\n        </DataArray>\n");
 	fprintf(fp,"        <DataArray type='Int32' Name='offsets' format='ascii'>\n");
	for(int i=0, ii=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#else
		{
			fprintf(fp,"%d ",ii+1);
			ii++;
		}
#endif
	}
 	fprintf(fp,"\n        </DataArray>\n");
 	fprintf(fp,"        <DataArray type='UInt8' Name='types' format='ascii'>\n");
	for(int i=0; i<nP; i++) {
#ifdef ONLY_FREE_WALL
		if(Bc[i] == SRF || Typ[i] == WLL)
#else
		{
			fprintf(fp,"1 ");
		}
#endif
	}
 	fprintf(fp,"\n        </DataArray>\n");
 	fprintf(fp,"      </Cells>\n");
 	fprintf(fp,"    </Piece>\n");
 	fprintf(fp,"  </UnstructuredGrid>\n");
 	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
	
	iF++;
}

// Write header for vtu files
void WrtPvd(void) {
	char pvd_filename[256];
	sprintf(pvd_filename, OUT_FOLDER".pvd");
	fp = fopen(pvd_filename, "w");
	int nIter = ceil(FIN_TIM/DT);

	//fprintf(fp,"<VTKFile type=""Collection"" version=""0.1"" byte_order=""LittleEndian"">\n");
	fprintf(fp,"<VTKFile type='Collection' version='0.1' byte_order='LittleEndian'>\n");
	fprintf(fp,"  <Collection>\n");
	int j = 0;
	for(int i=0;i<nIter;i++) {
		if(i % OPT_FQC == 0) {
			double timeStep = DT*i;
 			fprintf(fp,"    <DataSet timestep='%.6f' group='A' part='0' file='%s/output%05d.vtu'/>\n",timeStep,OUT_FOLDER,j);
 			j++;
		}
	}
	fprintf(fp,"  </Collection>\n");
	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
}

// Allocation of buckets
void AlcBkt(void) {
	if(DIM == 2) {
		reS = reS2D;
		reL = reL2D;
	}
	else {
		reS = reS3D;
		reL = reL3D;
	}
	reS = PCL_DST*reS;							// Influence radius small
	reL = PCL_DST*reL;							// Influence radius large
	reS2 = reS*reS;
	reL2 = reL*reL;
	eps_reS = EPS_RE*reS2/4;
	re_rep = PCL_DST*REP_RE;					// Influence radius for repulsive force
	DB = reL*(1.0+CRT_NUM);						// Length of one bucket side
	DB2 = DB*DB;	DBinv = 1.0/DB;
	nBx = (int)((MAX_X - MIN_X)*DBinv) + 3;		// Number of buckets in the x direction in the analysis domain
	nBy = (int)((MAX_Y - MIN_Y)*DBinv) + 3;		// Number of buckets in the y direction in the analysis domain
	nBz = (int)((MAX_Z - MIN_Z)*DBinv) + 3;		// Number of buckets in the z direction in the analysis domain
	nBxy = nBx*nBy;
	nBxyz = nBx*nBy*nBz;						// Number of buckets in analysis area
	printf("nBx:%d  nBy:%d  nBz:%d  nBxy:%d  nBxyz:%d\n",nBx,nBy,nBz,nBxy,nBxyz);
	bfst = (int*)malloc(sizeof(int) * nBxyz);	// First particle number stored in the bucket
	blst = (int*)malloc(sizeof(int) * nBxyz);	// Last particle number stored in the bucket
	nxt  = (int*)malloc(sizeof(int) * nP);		// Next particle number in the same bucket
}

// Set parameters
void SetPara(void) {
	n0S = n0L = n0Grad = lmd = nNeigh0 = 0.0;
	int lmin = ceil(reL/PCL_DST) + 1;
	int lmax = ceil(reL/PCL_DST) + 2;
	int flag2D = 0;
	int flag3D = 1;
	if(DIM == 2) {
		flag2D = 1;
		flag3D = 0;
	}
	for(int ix= -lmin;ix<lmax;ix++) {
	for(int iy= -lmin;iy<lmax;iy++) {
	for(int iz= -lmin*flag3D;iz<lmax*flag3D+flag2D;iz++) {
		double x = PCL_DST* (double)ix;
		double y = PCL_DST* (double)iy;
		double z = PCL_DST* (double)iz;
		double dst2 = x*x+y*y+z*z;
		if(dst2 <= reL2) {
			if(dst2 == 0.0)continue;
			double dst = sqrt(dst2);
			n0L += WEI(dst, reL);				// Initial particle number density (large)
			lmd += dst2 * WEI(dst, reL);
			nNeigh0 += 1;						// Initial number of neighboors
			if(dst2 <= reS2) {
				n0S += WEI(dst, reS);			// Initial particle number density (small)
				n0Grad += WEI_GRAD(dst, reS);	// Initial particle number density (gradient operator)
			}
		}
	}}}
	lmd = lmd/n0L;					// Coefficient λ of Laplacian model
	A1 = 2.0*KNM_VS1*DIM/(n0L*lmd);	// Coefficient used to calculate viscosity term
	A1_M = 2.0*DIM/(n0L*lmd);		// Coefficient used to calculate viscosity term Multiphase
	A2 = SND*SND/n0S;				// Coefficient used to calculate pressure E-MPS
	A3 = -DIM/n0Grad;				// Coefficient used to calculate pressure gradient term
	A4 = SND*SND;					// Coefficient used to calculate pressure WC-MPS
	A5 = DRI*PCL_DST/n0S;			// Coefficient used to adjust velocity
	A6 = COEF_A*PCL_DST*PCL_DST*CRT_NUM*MACH;	// Coefficient used to adjust velocity
	Dns[FLD]=DNS_FLD;			Dns[WLL]=DNS_WLL;
	invDns[FLD]=1.0/DNS_FLD;	invDns[WLL]=1.0/DNS_WLL;
	rlim = PCL_DST * DST_LMT_RAT;	// A distance that does not allow further access between particles
	rlim2 = rlim*rlim;
	COL = 1.0 + COL_RAT;
	iLP=iF=0;						// Number of iterations // File number
	TIM=0.0;						// Simulation time

	std::cout << "lo: " << PCL_DST << " m, dt: " << DT << " s, PND0Small: " << n0S << " PND0Large: " << n0L << " PND0Grad: " << n0Grad << " lambda: " << lmd << std::endl;
}

// Initial PND and number of neighbors
void InitPNDnNeigh_omp(void) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double w_sum = 0.0;
		numNeigh[i] = 0;
		dev[i*3] = dev[i*3+1] = dev[i*3+2] = 0;
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {

				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					numNeigh[i] += 1;
					if(dstij2 < reS2) {
						double dst = sqrt(dstij2);
						double wS = WEI(dst, reS);
						pndi[i] += wS;

						dev[i*3  ] += v0ij*wS/PCL_DST;
						dev[i*3+1] += v1ij*wS/PCL_DST;
						dev[i*3+2] += v2ij*wS/PCL_DST;
						w_sum += wS;
				}}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		// Add PND due Wall polygon
		pndi[i] += niw[i];
		if(Typ[i] == WLL)
			pndi[i] = n0S;
		pndS[i] = pndi[i];
		// Add Number of neighboors due Wall polygon
		numNeigh[i] += numNeighw[i];

		if(w_sum > 0.0001) {
			//dev[i*3  ] /= pndS[i];
			//dev[i*3+1] /= pndS[i];
			//dev[i*3+2] /= pndS[i];
			dev[i*3  ] /= w_sum;
			dev[i*3+1] /= w_sum;
			dev[i*3+2] /= w_sum;
		}

		devMod2[i] = dev[i*3]*dev[i*3]+dev[i*3+1]*dev[i*3+1]+dev[i*3+2]*dev[i*3+2];

		//devXnormal[i] = dev[i*3]*NormalWall[i*3]+dev[i*3+1]*NormalWall[i*3+1]+dev[i*3+2]*NormalWall[i*3+2];
		if(dev[i*3]*NormalWall[i*3]+dev[i*3+1]*NormalWall[i*3+1]+dev[i*3+2]*NormalWall[i*3+2] < 0.0)
			devXnormal[i] = 1;
		else
			devXnormal[i] = -1;
	}}
}

// Make buckets
void MkBkt(void) {
	for(int i=0;i< nBxyz ;i++) {	bfst[i] = -1;	}
	for(int i=0;i< nBxyz ;i++) {	blst[i] = -1;	}
	for(int i=0;i< nP ;i++) {	nxt[i] = -1;	}
	for(int i=0; i<nP; i++) {
		if(Typ[i] == GST)continue;
		int ix = (int)((Pos[i*3  ] - MIN_X)*DBinv) +1;
		int iy = (int)((Pos[i*3+1] - MIN_Y)*DBinv) +1;
		int iz = (int)((Pos[i*3+2] - MIN_Z)*DBinv) +1;
		int ib = iz*nBxy + iy*nBx + ix;
		int j = blst[ib];
		blst[ib] = i;
		if(j == -1) {	bfst[ib] = i;	}
		else {			nxt[j] = i;}
	}
}

// Laplacian of viscosity
void VscTrm_omp() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
//	if(Typ[i] == FLD) {
	if(Typ[i] != GST) {
		double meu_i = MEU[i];
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wL = WEI(dst, reL);
					NEU = 2 * meu_i * MEU[j] / (meu_i + MEU[j]);
					//NEU = KNM_VS2 * DNS_FL2;
					if(PTYPE[i] == 1) NEU = NEU/DNS_FL1;
					else NEU = NEU/DNS_FL2;
//					NEU = NEU/RHO[i];
					//if((NEUt[i] + NEUt[j]) > 0) NEU = NEU + (2 * NEUt[i] * RHO[j] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j])) / RHO[i];
					// Original
//					Acc_x +=(Vel[j*3  ]-vec_ix)*w;
//					Acc_y +=(Vel[j*3+1]-vec_iy)*w;
//					Acc_z +=(Vel[j*3+2]-vec_iz)*w;
					// Modified
					Acc_x +=(Vel[j*3  ]-vec_ix)*wL*NEU;
					Acc_y +=(Vel[j*3+1]-vec_iy)*wL*NEU;
					Acc_z +=(Vel[j*3+2]-vec_iz)*wL*NEU;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		// Original
//		Acc[i*3  ]=Acc_x*A1 + G_X;
//		Acc[i*3+1]=Acc_y*A1 + G_Y;
//		Acc[i*3+2]=Acc_z*A1 + G_Z;
		// Modified
		//if(TIM > 0.3) {
		// A1_M = 2.0*DIM/(n0L*lmd);
		Acc[i*3  ]=A1_M*Acc_x + G_X;
		Acc[i*3+1]=A1_M*Acc_y + G_Y;
		Acc[i*3+2]=A1_M*Acc_z + G_Z;
		//}		
		AccStar[i*3  ]=Acc[i*3  ];
		AccStar[i*3+1]=Acc[i*3+1];
		AccStar[i*3+2]=Acc[i*3+2];
	}}
}

// Prediction of pressure gradient
void PrdPrsGrdTrm_omp() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
//	if(Typ[i] == FLD) {
	if(Typ[i] != GST) {
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double Pi = Prs[i];			double ni = pndi[i];		double pre_min = Pi;
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		if(GRD_TYP == 0 || GRD_TYP == 2) {
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					if(pre_min > Prs[j]) pre_min = Prs[j];
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}}
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wS = WEI_GRAD(dst, reS);
					if(GRD_TYP == 0)
						wS *= (Prs[j] - pre_min)/dstij2;
					else if(GRD_TYP == 1)
						wS *= (Prs[j] + Pi)/dstij2;
					else if(GRD_TYP == 2)
						wS *= (Prs[j] + Pi - 2*pre_min)/dstij2;
					else if(GRD_TYP == 3) {
						double nj = pndi[j];
						if(ni > 0.0001 && nj > 0.0001)
							wS *= (ni*Prs[j]/nj + nj*Pi/ni)/dstij2;
					}
					Acc_x += v0ij*wS;	Acc_y += v1ij*wS;	Acc_z += v2ij*wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		// A3 is a negative cte (-DIM/noGrad)
		// Original
//		Acc[i*3  ]+=(1-RLX_PRS)*Acc_x*invDns[FLD]*A3;
//		Acc[i*3+1]+=(1-RLX_PRS)*Acc_y*invDns[FLD]*A3;
//		Acc[i*3+2]+=(1-RLX_PRS)*Acc_z*invDns[FLD]*A3;
		// Modified
		Acc[i*3  ]+=(1-RLX_PRS)*Acc_x*A3/RHO[i];
		Acc[i*3+1]+=(1-RLX_PRS)*Acc_y*A3/RHO[i];
		Acc[i*3+2]+=(1-RLX_PRS)*Acc_z*A3/RHO[i];
	}}
}

// Prediction of pressure gradient (Polygon wall)
void WallPrdPrsGrdTrm_omp() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {

	//Nw[i]=1; // Only to show particles near polygon

	if(Typ[i] == FLD && Nw[i] == 1) {

		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double Pi = Prs[i];			double ni = pndi[i];		double pre_min = Pi;
		
		// Wall gradient Mitsume`s model
	    double Rref_i[9], normaliw[3], normaliwSqrt;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);
	    if(normaliwSqrt > 0.00000001) {
	    	normaliw[0] = normaliw[0]/normaliwSqrt;
	    	normaliw[1] = normaliw[1]/normaliwSqrt;
	    	normaliw[2] = normaliw[2]/normaliwSqrt;
	    }
	    else {
	    	normaliw[0] = 0.0;
	    	normaliw[1] = 0.0;
	    	normaliw[2] = 0.0;
	    }
	    //  Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
	    Rref_i[0] = 1.0 - 2*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2*normaliw[0]*normaliw[2];
		Rref_i[3] = 0.0 - 2*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2*normaliw[1]*normaliw[2];
		Rref_i[6] = 0.0 - 2*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2*normaliw[2]*normaliw[2];
		// Taylor pressure Pj
		double Rai[3];
		Rai[0] = Rref_i[0]*AccStar[i*3] + Rref_i[1]*AccStar[i*3+1] + Rref_i[2]*AccStar[i*3+2];
		Rai[1] = Rref_i[3]*AccStar[i*3] + Rref_i[4]*AccStar[i*3+1] + Rref_i[5]*AccStar[i*3+2];
		Rai[2] = Rref_i[6]*AccStar[i*3] + Rref_i[7]*AccStar[i*3+1] + Rref_i[8]*AccStar[i*3+2];

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		if(GRD_TYP == 0 || GRD_TYP == 2) {
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					if(pre_min > Prs[j]) pre_min = Prs[j];
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}}
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {

					double dst = sqrt(dstimj2);
					double wS = WEI_GRAD(dst, reS);
					
					// Taylor pressure Pj
					double Pj;
//					Pj = Pi + RHO[i]*(Rai[0]*v0 + Rai[1]*v1 + Rai[2]*v2);
					Pj = Prs[j];

					if(GRD_TYP == 0)
						wS *= (Pj - pre_min)/dstimj2;//(Prs[j] - pre_min)/dstimj2;
					else if(GRD_TYP == 1)
						wS *= (Pj + Pi)/dstimj2;//(Prs[j] + Pi)/dstimj2;
					else if(GRD_TYP == 2)
						wS *= (Pj + Pi - 2*pre_min)/dstimj2;//(Prs[j] + Pi - 2*pre_min)/dstimj2;
					else if(GRD_TYP == 3) {
						double nj = pndi[j];
						if(ni > 0.0001 && nj > 0.0001)
							wS *= (ni*Pj/nj + nj*Pi/ni)/dstimj2;
					}
					Acc_x += v0imj*wS;	Acc_y += v1imj*wS;	Acc_z += v2imj*wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		// Add "i" contribution ("i" is a neighboor of "mirror i")
	  	double v0imi = pos_ix - pos_mix;
		double v1imi = pos_iy - pos_miy;
		double v2imi = pos_iz - pos_miz;
		double dstimi2 = v0imi*v0imi+v1imi*v1imi+v2imi*v2imi;
		if(dstimi2 < reS2) {

			double dst = sqrt(dstimi2);
			double wS = WEI_GRAD(dst, reS);

			// Taylor pressure Pj
			double Pj;
//			Pj = Pi + RHO[i]*(Rai[0]*v0 + Rai[1]*v1 + Rai[2]*v2);
			Pj = Pi;

			if(GRD_TYP == 0)
				wS *= (Pj - pre_min)/dstimi2;//(Pi - pre_min)/dstimi2
			else if(GRD_TYP == 1)
				wS *= (Pj + Pi)/dstimi2;//(Pi + Pi)/dstimi2
			else if(GRD_TYP == 2)
				wS *= (Pj + Pi - 2*pre_min)/dstimi2;//(Pi + Pi - 2*pre_min)/dstimi2
			else if(GRD_TYP == 3) {
				double nj = pndi[i];
				if(ni > 0.0001 && nj > 0.0001)
					wS *= (ni*Pj/nj + nj*Pi/ni)/dstimi2;
			}
			Acc_x += v0imi*wS;	Acc_y += v1imi*wS;	Acc_z += v2imi*wS;
	  	}

		// Repulsive force
		double rpsForce[3];
		rpsForce[0]=rpsForce[1]=rpsForce[2] = 0.0;

		if(REP_FOR == 0) {
			// Parallel analysis system for free-surface flow using MPS method with explicitly represented polygon wall boundary model
			// https://doi.org/10.1007/s40571-019-00269-6
			if(normaliwSqrt < re_rep && normaliwSqrt != 0) {
				double wijRep = RHO[i]/(DT*DT)*(re_rep-normaliwSqrt);
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
				//if(i == 6) {
					//printf("x:%lf y:%lf z:%lf x:%lf y:%lf z:%lf\n", Pos[i*3],Pos[i*3+1],Pos[i*3+2],mirrorPos[i*3],mirrorPos[i*3+1],mirrorPos[i*3+2]);
					//printf("nx:%lf ny: %lf nz: %lf nN:%lf\n", normaliw[0],normaliw[1],normaliw[2],normaliwSqrt);
					//printf("Fx:%lf Fy:%lf Fz:%lf\n", rpsForce[0],rpsForce[1],rpsForce[2]);
				//}
			}
		}
		else if(REP_FOR == 1) {
			// Explicitly represented polygon wall boundary model for the explicit MPS method
			// https://doi.org/10.1007/s40571-015-0037-8
			if(normaliwSqrt < re_rep && normaliwSqrt != 0) {
				double wijRep = MIT_REP*WEI_GRAD(normaliwSqrt, re_rep);
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
				//if(i == 6) {
					//printf("x:%lf y:%lf z:%lf x:%lf y:%lf z:%lf\n", Pos[i*3],Pos[i*3+1],Pos[i*3+2],mirrorPos[i*3],mirrorPos[i*3+1],mirrorPos[i*3+2]);
					//printf("nx:%lf ny: %lf nz: %lf nN:%lf\n", normaliw[0],normaliw[1],normaliw[2],normaliwSqrt);
					//printf("Fx:%lf Fy:%lf Fz:%lf\n", rpsForce[0],rpsForce[1],rpsForce[2]);
				//}
			}
		}
		else if(REP_FOR == 2) {
			// Simulating Free Surface Flows with SPH
			// https://doi.org/10.1006/jcph.1994.1034
			if(normaliwSqrt < re_rep && normaliwSqrt != 0) {
				double maxVel2 = vMax*vMax;
				double R1 = (re_rep/normaliwSqrt)*(re_rep/normaliwSqrt);
				double R2 = R1*R1;
				double wijRep = (LNJ_REP*maxVel2/normaliwSqrt)*(R2-R1)*RHO[i];
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
			}
		}
	  	else {
			// SPH particle boundary forces for arbitrary boundaries 
			// https://doi.org/10.1016/j.cpc.2009.05.008
			if(normaliwSqrt < re_rep && normaliwSqrt != 0) {
				double maxVel2 = vMax*vMax;
				double W1 = (1+3/2*normaliwSqrt/(re_rep));
				double W2 = (1-normaliwSqrt/(re_rep))*(1-normaliwSqrt/(re_rep))*(1-normaliwSqrt/(re_rep));
				double wijRep = (KJT_REP*maxVel2/(normaliwSqrt - 0.0*PCL_DST))*(1/8)*(W1)*(W2)*RHO[i];
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
			}
	  	}
		// A3 is a negative cte (-DIM/noGrad)
		// Original
//		Acc[i*3  ] += ((1-RLX_PRS)*(Rref_i[0]*Acc_x + Rref_i[1]*Acc_y + Rref_i[2]*Acc_z)*A3 - rpsForce[0])*invDns[FLD];
//		Acc[i*3+1] += ((1-RLX_PRS)*(Rref_i[3]*Acc_x + Rref_i[4]*Acc_y + Rref_i[5]*Acc_z)*A3 - rpsForce[1])*invDns[FLD];
//		Acc[i*3+2] += ((1-RLX_PRS)*(Rref_i[6]*Acc_x + Rref_i[7]*Acc_y + Rref_i[8]*Acc_z)*A3 - rpsForce[2])*invDns[FLD];
		// Modified
		Acc[i*3  ] += ((1-RLX_PRS)*(Rref_i[0]*Acc_x + Rref_i[1]*Acc_y + Rref_i[2]*Acc_z)*A3 - rpsForce[0])/RHO[i];
		Acc[i*3+1] += ((1-RLX_PRS)*(Rref_i[3]*Acc_x + Rref_i[4]*Acc_y + Rref_i[5]*Acc_z)*A3 - rpsForce[1])/RHO[i];
		Acc[i*3+2] += ((1-RLX_PRS)*(Rref_i[6]*Acc_x + Rref_i[7]*Acc_y + Rref_i[8]*Acc_z)*A3 - rpsForce[2])/RHO[i];

		//Fwall[i*3  ] =  (Rref_i[0]*Acc_x + Rref_i[1]*Acc_y + Rref_i[2]*Acc_z)*invDns[FLD]*A3 - rpsForce[0]*invDns[FLD];
		//Fwall[i*3+1] =  (Rref_i[3]*Acc_x + Rref_i[4]*Acc_y + Rref_i[5]*Acc_z)*invDns[FLD]*A3 - rpsForce[1]*invDns[FLD];
		//Fwall[i*3+2] =  (Rref_i[6]*Acc_x + Rref_i[7]*Acc_y + Rref_i[8]*Acc_z)*invDns[FLD]*A3 - rpsForce[2]*invDns[FLD];
	}}
}

// Update velocity and position
void UpPcl1_omp() {
#pragma omp parallel for
	for(int i=0; i<nP; i++) {
		if(Typ[i] != GST) {
			Vel[i*3  ] +=Acc[i*3  ]*DT;	Vel[i*3+1] +=Acc[i*3+1]*DT;	Vel[i*3+2] +=Acc[i*3+2]*DT;
			if(Typ[i] == FLD) {
				Pos[i*3  ] +=Vel[i*3  ]*DT;	Pos[i*3+1] +=Vel[i*3+1]*DT;	Pos[i*3+2] +=Vel[i*3+2]*DT;
				ChkPcl(i);
			}
		}
		Acc[i*3]=Acc[i*3+1]=Acc[i*3+2]=0.0;
		F1[ i*3]=F1[ i*3+1]=F1[ i*3+2]=0.0;
		F2[ i*3]=F2[ i*3+1]=F2[ i*3+2]=0.0;
		dev[i*3]=dev[i*3+1]=dev[i*3+2]=0.0;
		numNeighw[i]=Nw[i]=0;
		// Set squared distance of particle to triangle mesh to ~infinite
		riw2[i] = 10e8*PCL_DST;
		numNeighFree[i]=0.0;

#ifndef POLYGON_ON
		// Set mirrored particle to ~infinite if wall particles are used
		mirrorPos[i*3  ] = 10e8*PCL_DST; mirrorPos[i*3+1] = 10e8*PCL_DST; mirrorPos[i*3+2] = 10e8*PCL_DST;
#endif

	}
}

// Check collisions between particles
void ChkCol_omp() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] == FLD) {
//		double mi = Dns[Typ[i]];
		double mi;
		if(PTYPE[i] == 1) mi = DNS_FL1;
		else mi = DNS_FL2;

		double pos_ix =  Pos[i*3  ];double pos_iy =  Pos[i*3+1];double pos_iz =  Pos[i*3+2];
		double vec_ix =  Vel[i*3  ];double vec_iy =  Vel[i*3+1];double vec_iz =  Vel[i*3+2];
		double vec_ix2 = Vel[i*3  ];double vec_iy2 = Vel[i*3+1];double vec_iz2 = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < rlim2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double fDT = (vec_ix-Vel[j*3  ])*v0ij+(vec_iy-Vel[j*3+1])*v1ij+(vec_iz-Vel[j*3+2])*v2ij;
					if(fDT > 0.0) {
//						double mj = Dns[Typ[j]];
						double mj;
						if(PTYPE[j] == 1) mj = DNS_FL1;
						else mj = DNS_FL2;

						fDT *= COL*mj/(mi+mj)/dstij2;
						vec_ix2 -= v0ij*fDT;		vec_iy2 -= v1ij*fDT;		vec_iz2 -= v2ij*fDT;
					}
					/*
					double fDT = (Vel[j*3  ]-vec_ix)*v0+(Vel[j*3+1]-vec_iy)*v1+(Vel[j*3+2]-vec_iz)*v2;
					double mj = Dns[Typ[j]];
					fDT *= COL*mj/(mi+mj)/dst2;
					vec_ix2 += v0*fDT;		vec_iy2 += v1*fDT;		vec_iz2 += v2*fDT;
					*/
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		Acc[i*3  ]=vec_ix2;	Acc[i*3+1]=vec_iy2;	Acc[i*3+2]=vec_iz2;

		AccStar[i*3  ]=vec_ix2;	AccStar[i*3+1]=vec_iy2;	AccStar[i*3+2]=vec_iz2;

	}}
#pragma omp parallel for
	for(int i=0; i<nP; i++) {
		// CHANGED !!!
		//Pos[i*3  ]+=(Acc[i*3  ]-Vel[i*3  ])*DT; Pos[i*3+1]+=(Acc[i*3+1]-Vel[i*3+1])*DT; Pos[i*3+2]+=(Acc[i*3+2]-Vel[i*3+2])*DT;
		Vel[i*3  ]=Acc[i*3  ];	Vel[i*3+1]=Acc[i*3+1];	Vel[i*3+2]=Acc[i*3+2];

		//Velk[i*3  ]=Vel[i*3  ];	Velk[i*3+1]=Vel[i*3+1];	Velk[i*3+2]=Vel[i*3+2];
		//Pos[i*3  ]=Posk[i*3  ]+Vel[i*3  ]*DT; Pos[i*3+1]=Posk[i*3+1]+Vel[i*3+1]*DT; Pos[i*3+2]=Posk[i*3+2]+Vel[i*3+2]*DT;
		Acc[i*3  ]=0.0;	Acc[i*3+1]=0.0;	Acc[i*3+2]=0.0;	
	}
}

// Free-surface particles. NPCD (Polygon wall)
void WallNPCD_omp() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
	if(Typ[i] == FLD && Nw[i] == 1) {
		double dev_x = 0.0;	double dev_y = 0.0;	double dev_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		// Wall gradient Mitsume`s model
	    double Rref_i[9], normaliw[3], normaliwSqrt;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);

	    if(normaliwSqrt > 0.00000001) {
	    	normaliw[0] = normaliw[0]/normaliwSqrt;
	    	normaliw[1] = normaliw[1]/normaliwSqrt;
	    	normaliw[2] = normaliw[2]/normaliwSqrt;
	    }
	    else {
	    	normaliw[0] = 0;
	    	normaliw[1] = 0;
	    	normaliw[2] = 0;
	    }

	    //  Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
	    Rref_i[0] = 1.0 - 2*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2*normaliw[0]*normaliw[2];
		Rref_i[3] = 0.0 - 2*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2*normaliw[1]*normaliw[2];
		Rref_i[6] = 0.0 - 2*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2*normaliw[2]*normaliw[2];

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstimj2);
					double wS = WEI_GRAD(dst, reS);
					dev_x += v0imj*wS;
					dev_y += v1imj*wS;
					dev_z += v2imj*wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}

		}}}
		// Add "i" contribution ("i" is a neighboor of "mirror i")
	  	double v0imi = pos_ix - pos_mix;
		double v1imi = pos_iy - pos_miy;
		double v2imi = pos_iz - pos_miz;
		double dstimi2 = v0imi*v0imi+v1imi*v1imi+v2imi*v2imi;
		if(dstimi2 < reS2) {
			double dst = sqrt(dstimi2);
			double wS = WEI_GRAD(dst, reS);
			dev_x += v0imi*wS;
			dev_y += v1imi*wS;
			dev_z += v2imi*wS;
	  	}
	  	dev[i*3  ] += Rref_i[0]*dev_x + Rref_i[1]*dev_y + Rref_i[2]*dev_z;
		dev[i*3+1] += Rref_i[3]*dev_x + Rref_i[4]*dev_y + Rref_i[5]*dev_z;
		dev[i*3+2] += Rref_i[6]*dev_x + Rref_i[7]*dev_y + Rref_i[8]*dev_z;
	}}
}

// Compute PND
void MkPND_omp() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double ni = 0.0; double w_sum = 0.0;
		numNeigh[i] = 0;
		// Add Number of neighboors due Wall polygon
		numNeigh[i] += numNeighw[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					numNeigh[i] += 1;
					//double dst = sqrt(dst2);
					//double wL = WEI(dst, reL/PCL_DST);
					//dev[i*3  ] += v0*wL/PCL_DST;
					//dev[i*3+1] += v1*wL/PCL_DST;
					//dev[i*3+2] += v2*wL/PCL_DST;
					//w_sum += wL;
					if(dstij2 < reS2) {
						double dst = sqrt(dstij2);
						double wS = WEI(dst, reS);
						ni += wS;
						dst = dst/PCL_DST;
						wS = WEI(dst, reS/PCL_DST);
						dev[i*3  ] += v0ij*wS/PCL_DST;
						dev[i*3+1] += v1ij*wS/PCL_DST;
						dev[i*3+2] += v2ij*wS/PCL_DST;
						w_sum += wS;
				}}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

		//		double mi = Dns[Typ[i]];
//		double mi;
//		if(PTYPE[i] == 1) mi = DNS_FL1;
//		else mi = DNS_FL2;

		if(PND_CAL==0 || PND_CAL==1)
//		if(PND_CAL==0)
		{
			// PND due particles and Wall polygon
			pndi[i] = ni + niw[i];
		}

//		if(Typ[i] == WLL) {
			// PND due particles and Wall polygon
//			pndi[i] = ni + niw[i];
//			if(pndi[i] < n0S)
//				pndi[i] = n0S;

//				pndi[i] = n0S*pow((Prs[i]*GAM/(mi*A4)+1),GAM);
//		}

		// Add PND due Wall polygon
		pndS[i] = ni + niw[i];
		// Prevent pndS[i] = 0
//		if(numNeigh[i]>1) {
		if(w_sum>0.0001) {
			//dev[i*3  ] /= pndS[i];
			//dev[i*3+1] /= pndS[i];
			//dev[i*3+2] /= pndS[i];
			dev[i*3  ] /= w_sum;
			dev[i*3+1] /= w_sum;
			dev[i*3+2] /= w_sum;
		}

		devMod2[i] = dev[i*3]*dev[i*3]+dev[i*3+1]*dev[i*3+1]+dev[i*3+2]*dev[i*3+2];

		//devXnormal[i] = dev[i*3]*NormalWall[i*3]+dev[i*3+1]*NormalWall[i*3+1]+dev[i*3+2]*NormalWall[i*3+2];
		if( dev[i*3]*NormalWall[i*3]+dev[i*3+1]*NormalWall[i*3+1]+dev[i*3+2]*NormalWall[i*3+2]< 0.0)
			devXnormal[i] = 1;
		else
			devXnormal[i] = -1;
		
		// First check based on particle number density
//		if(pndS[i] < PND_TRS*n0S)
//			Bc[i] = SRF;
//		else
//			Bc[i] = INR;

		// Boundary particle verification based on relative distance and weight (NPCD)
//		if(Bc[i] == SRF) {
//			double delta2 = DLT_TRS*DLT_TRS*PCL_DST*PCL_DST;
//			if(numNeigh[i] > 4 && devMod2[i] < delta2)
//			{
//				Bc[i] = INR;
				//printf(" INR %d \n", i);
//			}
//		}

//		if(pndS[i] < PND_TRS*n0S && numNeigh[i] < NGH_TRS*nNeigh0)
//			Bc[i] = SRF;
//		else
//			Bc[i] = INR;
	}}
}

// Diffusion term of density/PND (PND_CAL = 2)
// An enhanced weakly-compressible MPS method for free-surface flows
// https://doi.org/10.1016/j.cma.2019.112771
void DifTrm_omp() {
	// A1_M = 2.0*DIM/(n0L*lmd);
	double C1 = DFS_MPS*DT*SND*SND*A1_M/(n0L);
	double C2 = DFS_MPS*PCL_DST*SND*A1_M/(n0L);
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
//	if(Typ[i] == FLD) {
		double Di = 0.0; double DivV = 0.0; double flagDi = 1.0; double pndAux = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double M1[3][3];
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
				M1[i][j] = 0.0;
		double ni = pndi[i];
		if(ni < 0.0001) continue;
		double Pi = Prs[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
//				if(j != i && Typ[j] == FLD) {
					double dst = sqrt(dstij2);
					double wL = WEI(dst, reL);
					double nj = pndi[j];
					if(Typ[i] == FLD && Typ[j] == FLD) {
//					if(Typ[i] == FLD && Typ[j] == FLD && Bc[i] == INR) {
						// A1_M = 2.0*DIM/(n0L*lmd);
//						double pgh = RHO[i]*(G_X*v0+G_Y*v1+G_Z*v2);
//						double CB = RHO[i]*SND*SND/GAM;
//						double nijH = n0S*(pow((pgh+1)/CB,1/GAM)-1);
//						double nijH = n0S*(pow((pgh)/CB+1,1/GAM)-1);
//
					//	pow( ( pgh + 1 ) / CB - 1 , 1  )

						//double CB = SND*SND*RHO[i];
						//double nijH = n0S*((PijH+1)/CB-1);
//						if(isnan(nijH) == 0)
//							Di += C1*(nj - ni - nijH)*wL;
//						else
							////////Di += C1*(nj - ni)*wL;
						//if(isnan(nijH) == 1)
						//if(i == 200)
						//	printf(" pgh %e CB %e ni %e nj %e nijH %e res %e \n", pgh, CB, ni, nj, nijH, n0S*(pow((pgh+1)/CB,1/GAM)-1));
						// PND
//						Di += C1*(nj-ni)*wL;
						//Di += C2*(nj-ni)*wL;
						// Pressure
						// Delta Voronoi smoothed particle hydrodynamics, δ-VSPH
						// https://doi.org/10.1016/j.jcp.2019.109000
						double pgh = -RHO[i]*(G_X*v0ij+G_Y*v1ij+G_Z*v2ij);
						double Pj = Prs[j];
						Di += DT/RHO[i]*A1_M*(Pj-Pi-pgh)*wL;
						//Di += (PCL_DST/SND)/RHO[i]*A1_M*(Pj-Pi+pgh)*wL;
					}
					//else
					//	flagDi = 0.0;
					if(dstij2 < reS2) {
						double vijx = Vel[j*3  ]-vec_ix;
						double vijy = Vel[j*3+1]-vec_iy;
						double vijz = Vel[j*3+2]-vec_iz;
						double wS = WEI(dst, reS);
						if(ni > 0.0001)
							DivV += (DIM/n0S)*(nj/ni)*(vijx*v0ij+vijy*v1ij+vijz*v2ij)*wS/dstij2;

//						M1[0][0] += (DIM/n0S)*(nj/ni)*(v0*vijx)*wS/dst2; M1[0][1] += (DIM/n0S)*(nj/ni)*(v0*vijy)*wS/dst2; M1[0][2] += (DIM/n0S)*(nj/ni)*(v0*vijz)*wS/dst2;
//						M1[1][0] += (DIM/n0S)*(nj/ni)*(v1*vijx)*wS/dst2; M1[1][1] += (DIM/n0S)*(nj/ni)*(v1*vijy)*wS/dst2; M1[1][2] += (DIM/n0S)*(nj/ni)*(v1*vijz)*wS/dst2;
//						M1[2][0] += (DIM/n0S)*(nj/ni)*(v2*vijx)*wS/dst2; M1[2][1] += (DIM/n0S)*(nj/ni)*(v2*vijy)*wS/dst2; M1[2][2] += (DIM/n0S)*(nj/ni)*(v2*vijz)*wS/dst2;

//						if(Typ[i] == WLL)
//						if(Typ[i] == WLL || Bc[i] == SRF)
//							pndAux += wS;
				}}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}


//		DivV = 	CMr1[i*3  ]*M1[0][0] + CMr1[i*3+1]*M1[1][0] + CMr1[i*3+2]*M1[2][0] +
//				CMr2[i*3  ]*M1[0][1] + CMr2[i*3+1]*M1[1][1] + CMr2[i*3+2]*M1[2][1] +
//				CMr3[i*3  ]*M1[0][2] + CMr3[i*3+1]*M1[1][2] + CMr3[i*3+2]*M1[2][2];

		Acc[i*3] = pndi[i]*(1+DT*(-DivV+Di*flagDi));


//		if(isnan(DivV) || isnan(Di))
//			printf(" i %d \n", i);
		DIVi[i] = DivV;
		DIi[i] = Di;

		//Acc[i*3] = pndi[i]*(1+DT*(-(1-DFS_MPS)*DivV+Di*flagDi));
		//Acc[i*3] = pndS[i]*(1+DT*(-DivV+Di*flagDi));
		//Acc[i*3] =     n0S*(1+DT*(-DivV+Di*flagDi)); // Ruim
//		if(Typ[i] == WLL)
//		{
//			if(pndAux < n0S)
//				pndAux = n0S;
//			Acc[i*3] = pndAux;
//		}
//		if(Bc[i] == SRF)
//		{
//			Acc[i*3] = pndAux;
//		}
	}}
//#pragma omp parallel for
//	for(int i=0; i<nP; i++) {
//	if(Typ[i] != GST) {
/////	if(Typ[i] == FLD) {
//		pndi[i] = Acc[i*3];
//		Acc[i*3]=0.0;
//	}}
}

// Diffusion term of density/PND (Polygon wall) - Free-slip (PND_CAL = 2)
void WallSlipDifTrm_omp() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
	double ni = pndi[i];
	//if(Typ[i] == FLD && ni > 0.0001) {
	if(Typ[i] == FLD && Nw[i] == 1 && ni > 0.0001) {
//	if(Typ[i] == FLD) {
		double DivV = 0.0;
		//double Pi = Prs[i];
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
//		double vec_ix = Velk[i*3  ];	double vec_iy = Velk[i*3+1];	double vec_iz = Velk[i*3+2

		// Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
		double Rref_i[9], normaliw[3], normalMod2;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];
	    if(normalMod2 > 0.00000001) {
	    	double normalMod = sqrt(normalMod2);
	    	normaliw[0] = normaliw[0]/normalMod;
	    	normaliw[1] = normaliw[1]/normalMod;
	    	normaliw[2] = normaliw[2]/normalMod;
	    }
	    else {
	    	normaliw[0] = 0;
	    	normaliw[1] = 0;
	    	normaliw[2] = 0;
	    }

	    //  Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
	    Rref_i[0] = 1.0 - 2*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2*normaliw[0]*normaliw[2];
		Rref_i[3] = 0.0 - 2*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2*normaliw[1]*normaliw[2];
		Rref_i[6] = 0.0 - 2*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2*normaliw[2]*normaliw[2];

		// Mirror particle velocity vi' = Rref_i * vi
      	double vec_mix = (Rref_i[0]*vec_ix + Rref_i[1]*vec_iy + Rref_i[2]*vec_iz);
		double vec_miy = (Rref_i[3]*vec_ix + Rref_i[4]*vec_iy + Rref_i[5]*vec_iz);
		double vec_miz = (Rref_i[6]*vec_ix + Rref_i[7]*vec_iy + Rref_i[8]*vec_iz);

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
					if(j != i && Typ[j] != GST) {
						double dst = sqrt(dstimj2);
						double wS = WEI(dst, reS);
						double nj = pndi[j];
						double vijx = Vel[j*3  ]-vec_mix;
						double vijy = Vel[j*3+1]-vec_miy;
						double vijz = Vel[j*3+2]-vec_miz;
						DivV += (DIM/n0S)*(nj/ni)*(vijx*v0imj+vijy*v1imj+vijz*v2imj)*wS/dstimj2;
					}}
				j = nxt[j];
				if(j == -1) break;
				}
		}}}

		// Add "i" contribution ("i" is a neighboor of "mirror i")
		double v0imi = pos_ix - pos_mix;
		double v1imi = pos_iy - pos_miy;
		double v2imi = pos_iz - pos_miz;
		double dstimi2 = v0imi*v0imi+v1imi*v1imi+v2imi*v2imi;
		if(dstimi2 < reS2) {
			double dst = sqrt(dstimi2);
			double wS = WEI(dst, reS);
			double nj = pndi[i];
			double vijx = vec_ix-vec_mix;
			double vijy = vec_iy-vec_miy;
			double vijz = vec_iz-vec_miz;
			DivV += (DIM/n0S)*(nj/ni)*(vijx*v0imi+vijy*v1imi+vijz*v2imi)*wS/dstimi2;
	  	}

		Acc[i*3] += -pndi[i]*DT*DivV;
		//Acc[i*3] = n0S*(1+DT*(-DivV+Di*flagDi));
	}}
//#pragma omp parallel for
//	for(int i=0; i<nP; i++) {
/////	if(Typ[i] != GST) {
//	if(Typ[i] == FLD) {
//		pndi[i] += Acc[i*3];
//		Acc[i*3]=0.0;
//	}}
}

// Diffusion term of density/PND (Polygon wall) - No-slip (PND_CAL = 2)
void WallNoSlipDifTrm_omp() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
	double ni = pndi[i];
	//if(Typ[i] == FLD && ni > 0.0001) {
	if(Typ[i] == FLD && Nw[i] == 1 && ni > 0.0001) {
//	if(Typ[i] == FLD) {
		double DivV = 0.0;
		//double Pi = Prs[i];
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
//		double vec_ix = Velk[i*3  ];	double vec_iy = Velk[i*3+1];	double vec_iz = Velk[i*3+2

		// Inverse matrix Rinv_i = - I
		double Rinv_i[9], Rref_i[9], normaliw[3], normalMod2;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];
	    if(normalMod2 > 0.00000001) {
	    	double normalMod = sqrt(normalMod2);
	    	normaliw[0] = normaliw[0]/normalMod;
	    	normaliw[1] = normaliw[1]/normalMod;
	    	normaliw[2] = normaliw[2]/normalMod;
	    }
	    else {
	    	normaliw[0] = 0;
	    	normaliw[1] = 0;
	    	normaliw[2] = 0;
	    }

	    //  Inverse transformation matrix Rinv_i = - I
	    Rinv_i[0] = -1.0; Rinv_i[1] =  0.0; Rinv_i[2] =  0.0;
		Rinv_i[3] =  0.0; Rinv_i[4] = -1.0; Rinv_i[5] =  0.0;
		Rinv_i[6] =  0.0; Rinv_i[7] =  0.0; Rinv_i[8] = -1.0;

	    //  Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
	    Rref_i[0] = 1.0 - 2*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2*normaliw[0]*normaliw[2];
		Rref_i[3] = 0.0 - 2*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2*normaliw[1]*normaliw[2];
		Rref_i[6] = 0.0 - 2*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2*normaliw[2]*normaliw[2];

		double viwall[3], vtil[3];
		// Wall velocity (0 if fixed)
		viwall[0]=viwall[1]=viwall[2]=0.0;

#ifdef FORCED_ON
		if(meshID[i] == 2) {
			viwall[0] = velVWall[0];
			viwall[1] = velVWall[1];
			viwall[2] = velVWall[2];
		}
#endif

		// normal_iwall*v_iwall
		double dotnv = normaliw[0]*viwall[0] + normaliw[1]*viwall[1] + normaliw[2]*viwall[2];
		// vtil = vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}
		vtil[0] = vec_ix - 2*(viwall[0] - dotnv*normaliw[0]);
		vtil[1] = vec_iy - 2*(viwall[1] - dotnv*normaliw[1]);
		vtil[2] = vec_iz - 2*(viwall[2] - dotnv*normaliw[2]);
		// Mirror particle velocity vi' = Rinv_i * [vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}] 
      	double vec_mix = (Rinv_i[0]*vtil[0] + Rinv_i[1]*vtil[1] + Rinv_i[2]*vtil[2]);
		double vec_miy = (Rinv_i[3]*vtil[0] + Rinv_i[4]*vtil[1] + Rinv_i[5]*vtil[2]);
		double vec_miz = (Rinv_i[6]*vtil[0] + Rinv_i[7]*vtil[1] + Rinv_i[8]*vtil[2]);

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
					if(j != i && Typ[j] != GST) {
						double dst = sqrt(dstimj2);
						double wS = WEI(dst, reS);
						double nj = pndi[j];
						double vijx = -(Vel[j*3  ]-vec_mix);
						double vijy = -(Vel[j*3+1]-vec_miy);
						double vijz = -(Vel[j*3+2]-vec_miz);
						// Refelected rij' = Rref_i * ri'j
      					double v0m = (Rref_i[0]*v0imj + Rref_i[1]*v1imj + Rref_i[2]*v2imj);
						double v1m = (Rref_i[3]*v0imj + Rref_i[4]*v1imj + Rref_i[5]*v2imj);
						double v2m = (Rref_i[6]*v0imj + Rref_i[7]*v1imj + Rref_i[8]*v2imj);
						DivV += (DIM/n0S)*(nj/ni)*(vijx*v0m+vijy*v1m+vijz*v2m)*wS/dstimj2;
					}}
				j = nxt[j];
				if(j == -1) break;
				}
		}}}

		// Add "i" contribution ("i" is a neighboor of "mirror i")
		double v0imi = pos_ix - pos_mix;
		double v1imi = pos_iy - pos_miy;
		double v2imi = pos_iz - pos_miz;
		double dstimi2 = v0imi*v0imi+v1imi*v1imi+v2imi*v2imi;
		if(dstimi2 < reS2) {
			double dst = sqrt(dstimi2);
			double wS = WEI(dst, reS);
			double nj = pndi[i];
			double vijx = -(vec_ix-vec_mix);
			double vijy = -(vec_iy-vec_miy);
			double vijz = -(vec_iz-vec_miz);
			// Refelected rij' = Rref_i * ri'j
			double v0m = (Rref_i[0]*v0imi + Rref_i[1]*v1imi + Rref_i[2]*v2imi);
			double v1m = (Rref_i[3]*v0imi + Rref_i[4]*v1imi + Rref_i[5]*v2imi);
			double v2m = (Rref_i[6]*v0imi + Rref_i[7]*v1imi + Rref_i[8]*v2imi);
			DivV += (DIM/n0S)*(nj/ni)*(vijx*v0m+vijy*v1m+vijz*v2m)*wS/dstimi2;
	  	}

		Acc[i*3] += -pndi[i]*DT*DivV;
		//Acc[i*3] = n0S*(1+DT*(-DivV+Di*flagDi));
	}}
//#pragma omp parallel for
//	for(int i=0; i<nP; i++) {
/////	if(Typ[i] != GST) {
//	if(Typ[i] == FLD) {
//		pndi[i] += Acc[i*3];
//		Acc[i*3]=0.0;
//	}}
}

// Update PND (PND_CAL = 2)
void UpPND_omp(void) {
#pragma omp parallel for
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
//		if(Typ[i] == FLD)
			pndi[i] = Acc[i*3];
//		else
//		{
//			pndi[i] = Acc[i*3];
			/*
			double mi;
			if(PTYPE[i] == 1) 
				mi = DNS_FL1;
			else 
				mi = DNS_FL2;
			if(MPS_TYP == 0)
				pndi[i] = n0S*(Prs[i]/(mi*A4)+1);
			else if(MPS_TYP == 1)
				pndi[i] = n0S*pow(Prs[i]*GAM/(mi*A4)+1,GAM);
				*/
//		}
		Acc[i*3]=0.0;
	}}
}

// Mean PND at wall and dummy particles (PND_CAL = 2)
void MeanPNDWallDummy_omp(void) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
//	if(Typ[i] != GST) {
//	if(Typ[i] == WLL) {
	if(Typ[i] == WLL || Bc[i] == SRF) {
//	if(Bc[i] == SRF) {
		double PNDup = 0.0;
		double PNDdo = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wS = WEI(dst, reS);
					PNDup += pndi[j]*wS;
					PNDdo += wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		Acc[i*3  ] = PNDup;
		Acc[i*3+1] = PNDdo;
		//Acc[i*3] = n0S*(1+DT*(-DivV+Di*flagDi));
//	}}}
	}}
#pragma omp parallel for
	for(int i=0; i<nP; i++) {
//	if(Typ[i] != GST) {
//	if(Typ[i] == WLL) {
	if(Typ[i] == WLL || Bc[i] == SRF) {
//	if(Bc[i] == SRF) {
		// Prevent PNDdo = 0
//		if(numNeigh[i] < 1)
		if(Acc[i*3+1] < 0.0001)
			pndi[i] = Acc[i*3];
		else
			pndi[i] = Acc[i*3]/(Acc[i*3+1]);
//			pndi[i] = Acc[i*3]/(Acc[i*3+1] + 0.01*reS2/4);
	}
//	}
		Acc[i*3]=0.0;Acc[i*3+1]=0.0;
	}
}

// Mean PND (Polygon wall) (PND_CAL = 1)
void WallMeanPND_omp() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
//	if(Typ[i] == FLD) {
	//if(Typ[i] != GST) {
	if(Typ[i] != GST && Nw[i] == 1) {
		double PNDup = 0.0;
		double PNDdo = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstimj2);
					double wS = WEI(dst, reS);
					PNDup += pndi[j]*wS;
					PNDdo += wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}

		}}}
		// Add "i" contribution ("i" is a neighboor of "mirror i")
	  	double v0imi = pos_ix - pos_mix;
		double v1imi = pos_iy - pos_miy;
		double v2imi = pos_iz - pos_miz;
		double dstimi2 = v0imi*v0imi+v1imi*v1imi+v2imi*v2imi;
		if(dstimi2 < reS2) {
			double dst = sqrt(dstimi2);
			double wS = WEI(dst, reS);
			PNDup += pndi[i]*wS;
			PNDdo += wS;
	  	}
	  	Acc[i*3  ] = PNDup;
	  	Acc[i*3+1] = PNDdo;
	}}
}

// Mean PND (PND_CAL = 1)
void MeanPND_omp(void) {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
		double PNDup = 0.0;
		double PNDdo = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wS = WEI(dst, reS);
					PNDup += pndi[j]*wS;
					PNDdo += wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		Acc[i*3  ] += PNDup;
		Acc[i*3+1] += PNDdo;
		//Acc[i*3] = n0S*(1+DT*(-DivV+Di*flagDi));
	}}
#pragma omp parallel for
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
//	if(Typ[i] == FLD) {
		// Prevent PNDdo = 0
		if(numNeigh[i] < 1)
			pndi[i] = Acc[i*3];
		else
			pndi[i] = Acc[i*3]/Acc[i*3+1];

		Acc[i*3]=0.0;Acc[i*3+1]=0.0;
	}}
}

// Compute pressure EMPS (MPS_TYP = 0) and type of particle
void MkPrs_omp() {
// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
#pragma omp parallel for
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
/*
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double ni = 0.0;
		numNeigh[i] = 0;
		// Add Number of neighboors due Wall polygon
		numNeigh[i] += numNeighw[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					numNeigh[i] += 1;
					if(dstij2 < reS2) {
						double dst = sqrt(dstij2);
						double wS = WEI(dst, reS);
						ni += wS;
						dev[i*3  ] += v0ij*wS;
						dev[i*3+1] += v1ij*wS;
						dev[i*3+2] += v2ij*wS;
				}}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

		//		double mi = Dns[Typ[i]];
		double mi;
		if(PTYPE[i] == 1) mi = DNS_FL1;
		else mi = DNS_FL2;

		if(PND_CAL == 0 || PND_CAL == 1)
//		if(PND_CAL == 0)
		{
			// PND due particles and Wall polygon
			pndi[i] = ni + niw[i];
		}
//		if(Typ[i] == WLL) {
			// PND due particles and Wall polygon
//			pndi[i] = ni + niw[i];
//			if(pndi[i] < n0S)
//				pndi[i] = n0S;

//				pndi[i] = n0S*(Prs[i]/(mi*A4)+1);
//		}
		// Add PND due Wall polygon
		pndS[i] = ni + niw[i];
		// Prevent pndS[i] = 0
		if(numNeigh[i] > 1) {
			dev[i*3  ] /= pndS[i];
			dev[i*3+1] /= pndS[i];
			dev[i*3+2] /= pndS[i];
		}
		devMod2[i] = dev[i*3]*dev[i*3]+dev[i*3+1]*dev[i*3+1]+dev[i*3+2]*dev[i*3+2];

		
		if(dev[i*3]*NormalWall[i*3]+dev[i*3+1]*NormalWall[i*3+1]+dev[i*3+2]*NormalWall[i*3+2] < 0.0)
			devXnormal[i] = 1;
		//A2 = SND*SND/n0S
//		double pressure = 0.0;
*/
		// First check based on particle number density
		if(pndS[i] < PND_TRS*n0S)
			Bc[i] = SRF;
		else
			Bc[i] = INR;

		// Boundary particle verification based on relative distance and weight (NPCD)
		if(Bc[i] == SRF) {
			double delta2 = DLT_TRS*DLT_TRS*PCL_DST*PCL_DST;
			if(numNeigh[i] > 4 && devMod2[i] < delta2)
				Bc[i] = INR;
		}

		//		double mi = Dns[Typ[i]];
		double mi;
		if(PTYPE[i] == 1) mi = DNS_FL1;
		else mi = DNS_FL2;

		double pressure = 0.0;
		if(Bc[i] == INR)
			pressure = (pndi[i] - n0S) * A2 * mi;

//		if(pndS[i] < PND_TRS*n0S && numNeigh[i] < NGH_TRS*nNeigh0)
//			Bc[i] = SRF;
//		else
//		{
//			Bc[i] = INR;
//			pressure = (pndi[i] - n0S) * A2 * mi;
//		}
//#ifdef POLYGON_ON
//		if(pndi[i] > PND_TRS*n0S)
//#else
//		if(pndi[i] > PND_TRS*n0S || numNeigh[i] > NGH_TRS*nNeigh0)
//#endif

//		pressure = -mi*G_Z*(0.3-pos_iz);

		if(pressure < 0.0)
			pressure = 0.0;
		Prs[i] = pressure;
	}}
}

// Compute pressure WCMPS (MPS_TYP = 1) and type of particle
void MkPrsWc_omp() {
// Use #pragma omp parallel for schedule(dynamic,64) if there are "for" inside the main "for"
//#pragma omp parallel for schedule(dynamic,64)
#pragma omp parallel for
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
/*
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double ni = 0.0; double w_sum = 0.0;
		numNeigh[i] = 0;
		// Add Number of neighboors due Wall polygon
		numNeigh[i] += numNeighw[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					numNeigh[i] += 1;
					//double dst = sqrt(dstij2);
					//double wL = WEI(dst, reL/PCL_DST);
					//dev[i*3  ] += v0ij*wL/PCL_DST;
					//dev[i*3+1] += v1ij*wL/PCL_DST;
					//dev[i*3+2] += v2ij*wL/PCL_DST;
					//w_sum += wL;
					if(dstij2 < reS2) {
						double dst = sqrt(dstij2);
						double wS = WEI(dst, reS);
						ni += wS;
						dst = dst/PCL_DST;
						wS = WEI(dst, reS/PCL_DST);
						dev[i*3  ] += v0ij*wS/PCL_DST;
						dev[i*3+1] += v1ij*wS/PCL_DST;
						dev[i*3+2] += v2ij*wS/PCL_DST;
						w_sum += wS;
				}}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

		//		double mi = Dns[Typ[i]];
		double mi;
		if(PTYPE[i] == 1) mi = DNS_FL1;
		else mi = DNS_FL2;

		if(PND_CAL == 0 || PND_CAL == 1)
//		if(PND_CAL == 0)
		{
			// PND due particles and Wall polygon
			pndi[i] = ni + niw[i];
		}

//		if(Typ[i] == WLL) {
			// PND due particles and Wall polygon
//			pndi[i] = ni + niw[i];
//			if(pndi[i] < n0S)
//				pndi[i] = n0S;

//				pndi[i] = n0S*pow((Prs[i]*GAM/(mi*A4)+1),GAM);
//		}

		// Add PND due Wall polygon
		pndS[i] = ni + niw[i];
		// Prevent pndS[i] = 0
//		if(numNeigh[i] > 1) {
		if(w_sum > 0.0001) {
			//dev[i*3  ] /= pndS[i];
			//dev[i*3+1] /= pndS[i];
			//dev[i*3+2] /= pndS[i];
			dev[i*3  ] /= w_sum;
			dev[i*3+1] /= w_sum;
			dev[i*3+2] /= w_sum;
		}

		devMod2[i] = dev[i*3]*dev[i*3]+dev[i*3+1]*dev[i*3+1]+dev[i*3+2]*dev[i*3+2];

		//devXnormal[i] = dev[i*3]*NormalWall[i*3]+dev[i*3+1]*NormalWall[i*3+1]+dev[i*3+2]*NormalWall[i*3+2];
		if(dev[i*3]*NormalWall[i*3]+dev[i*3+1]*NormalWall[i*3+1]+dev[i*3+2]*NormalWall[i*3+2] < 0.0)
			devXnormal[i] = 1;
		else
			devXnormal[i] = -1;

		// A4 = SND*SND
		double pressure = 0.0;
*/		
		// First check based on particle number density
		if(pndS[i] < PND_TRS*n0S)
			Bc[i] = SRF;
		else
			Bc[i] = INR;

		// Boundary particle verification based on relative distance and weight (NPCD)
		if(Bc[i] == SRF) {
			double delta2 = DLT_TRS*DLT_TRS*PCL_DST*PCL_DST;
			if(numNeigh[i] > 4 && devMod2[i] < delta2)
			{
				Bc[i] = INR;
				//printf(" INR %d \n", i);
			}
		}


		if(pndS[i] < PND_TRS*n0S && numNeigh[i] < NGH_TRS*nNeigh0)
//		if(pndi[i] < PND_TRS*n0S && numNeigh[i] < NGH_TRS*nNeigh0)
			Bc[i] = SRF;
		else
			Bc[i] = INR;

		//		double mi = Dns[Typ[i]];
		double mi;
		if(PTYPE[i] == 1) mi = DNS_FL1;
		else mi = DNS_FL2;

		double pressure = 0.0;
		if(Bc[i] == INR)
//		if(Bc[i] == INR && Typ[i] == FLD)
			pressure = (mi*A4/GAM)*(pow(pndi[i]/n0S,GAM)-1);
		
		/*
		if(pndS[i] < PND_TRS*n0S && numNeigh[i] < NGH_TRS*nNeigh0)
			Bc[i] = SRF;
		else
		{
			Bc[i] = INR;
			pressure = (mi*A4/GAM)*(pow(pndi[i]/n0S,GAM)-1);
			//pressure = (ni - n0) * A2 * mi;
//			pressure = (pndi[i] - n0S) * A2 * mi;
		}
		*/
//#ifdef POLYGON_ON
//		if(pndi[i] > PND_TRS*n0S)
//#else
//		if(pndi[i] > PND_TRS*n0S || numNeigh[i] > NGH_TRS*nNeigh0)
//#endif

//		pressure = -mi*G_Z*(0.20 - 0.5*PCL_DST -pos_iz);
//		pressure = -mi*G_Z*(0.18 - 0.5*PCL_DST -pos_iz); // lat
		
		if(pressure < 0.0)
			pressure = 0.0;
		Prs[i] = pressure;
	}}
}

// Compute pressure. Wall and dummy particles
void MkPrsWallDummy_omp() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] == WLL) {
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double ni = 0.0;
		double pressure = 0.0;
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2 < reS2) {
//				if(j != i && Typ[j] != GST) {
				if(j != i && Typ[j] == FLD) {
					double dst = sqrt(dst2);
					double wS = WEI(dst, reS);
					ni += wS;
					pressure += (Prs[j] - RHO[j]*(G_X*v0+G_Y*v1+G_Z*v2))*wS;
					//pressure += (Prs[j] + RHO[j]*(G_X*v0+G_Y*v1+G_Z*v2))*wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		if(pressure < 0.0)
			pressure = 0.0;
		if(ni > 0)
			Prs[i] = pressure/ni;
		else
			Prs[i] = pressure;
	}}
}

// Compute pressure. Inner particles near polygon walls
void MkPrsNearWall_omp() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
	//if(Typ[i] == FLD && Prs[i] == 0 /*&& Nw[i] == 1*/) {
	if(Typ[i] == FLD && Prs[i] == 0 && Nw[i] == 1) {
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double pressure = 0.0;
		double sumWij = 0.0;
		int nTotal = 0;
		int nFree = 0;
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
//				if(dstij2 < 1.2*PCL_DST) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wS = WEI(dst, reS);
//					double wS = WEI_WEND(dst, 1.2*PCL_DST);
					sumWij += wS;
//					pressure += Prs[j]*wS;
					pressure += (Prs[j] - RHO[j]*(G_X*v0ij+G_Y*v1ij+G_Z*v2ij))*wS;
					nTotal += 1;
					if(Bc[j] == SRF)
						nFree += 1;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		
		if(nTotal > 0)
			numNeighFree[i] = double(nFree)/nTotal;
		else
			numNeighFree[i] = 1.0;

//		if(numNeighFree[i]<=0.5){
		if(pressure > 0){
	  		if(sumWij > 0.0001)
				Prs[i] = pressure/sumWij;
			else
				Prs[i] = pressure;
		}
//		}
	}}
}

// Determinant of matrix
double detMatrix (double M11, double M12, double M13, double M21, double M22, double M23, double M31, double M32, double M33)
{
    return (M11*M22*M33 + M12*M23*M31 + M13*M21*M32)
          - (M13*M22*M31 + M12*M21*M33 + M11*M23*M32);
}

// Inverse of matrix
int inverseMatrix (int dim, double &M11, double &M12, double &M13, double &M21, double &M22, double &M23, double &M31, double &M32, double &M33) {
	double M[3][3], Maux[3][3];

   	Maux[0][0] = M11;	Maux[0][1] = M12;	Maux[0][2] = M13;
   	Maux[1][0] = M21;	Maux[1][1] = M22;	Maux[1][2] = M23;
   	Maux[2][0] = M31;	Maux[2][1] = M32;	Maux[2][2] = M33;

   	if(dim == 2)
   		Maux[2][2] = 1.0;

   	// Convert matrix to identity
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++) {
            if(i == j) M[i][j] = 1.0;
            else M[i][j] = 0.0;
        }

    for(int k = 0; k < dim; k++) {
        
        if(fabs(Maux[k][k]) <= 1e-6) {
            M11 = 1.0;	M12 = 0.0;	M13 = 0.0;
   			M21 = 0.0;	M22 = 1.0;	M23 = 0.0;
   			M31 = 0.0;	M32 = 0.0;	M33 = 1.0;
            return 0;
        }
        
        for(int i = 0; i < dim; i++) {
            if(i == k)
                continue;
            
            double m = Maux[i][k]/Maux[k][k];
            
            for(int j = 0; j < dim; j++) {
                Maux[i][j] -= m*Maux[k][j];
                M[i][j] -= m*M[k][j];
            }
        }
    }
        
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            M[i][j] /= Maux[i][i];
    
    M11 = M[0][0];	M12 = M[0][1];	M13 = M[0][2];
   	M21 = M[1][0];	M22 = M[1][1];	M23 = M[1][2];
   	M31 = M[2][0];	M32 = M[2][1];	M33 = M[2][2];

    return 1;
}

// Correction matrix
void CorMtxTrm_omp() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
//	if(Typ[i] == FLD) {
	if(Typ[i] != GST) {
		CMr1[i*3  ] = 0.0; CMr1[i*3+1] = 0.0; CMr1[i*3+2] = 0.0;
		CMr2[i*3  ] = 0.0; CMr2[i*3+1] = 0.0; CMr2[i*3+2] = 0.0;
		CMr3[i*3  ] = 0.0; CMr3[i*3+1] = 0.0; CMr3[i*3+2] = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wS = WEI_GRAD(dst, reS);
					CMr1[i*3  ] += wS*v0ij*v0ij/dstij2;	CMr1[i*3+1] += wS*v0ij*v1ij/dstij2;	CMr1[i*3+2] += wS*v0ij*v2ij/dstij2;
					CMr2[i*3  ] += wS*v1ij*v0ij/dstij2;	CMr2[i*3+1] += wS*v1ij*v1ij/dstij2;	CMr2[i*3+2] += wS*v1ij*v2ij/dstij2;
					CMr3[i*3  ] += wS*v2ij*v0ij/dstij2;	CMr3[i*3+1] += wS*v2ij*v1ij/dstij2;	CMr3[i*3+2] += wS*v2ij*v2ij/dstij2;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

		// A3 is a negative cte (-DIM/noGrad)
		CMr1[i*3  ] *= -A3;		CMr1[i*3+1] *= -A3;		CMr1[i*3+2] *= -A3;
		CMr2[i*3  ] *= -A3;		CMr2[i*3+1] *= -A3;		CMr2[i*3+2] *= -A3;
		CMr3[i*3  ] *= -A3;		CMr3[i*3+1] *= -A3;		CMr3[i*3+2] *= -A3;

		// Inverse of the matrix
		int rcv = inverseMatrix(DIM,CMr1[i*3],CMr1[i*3+1],CMr1[i*3+2],CMr2[i*3],CMr2[i*3+1],CMr2[i*3+2],CMr3[i*3],CMr3[i*3+1],CMr3[i*3+2]);

//		if(i == 200) {
//			printf("\n X %e %e %e ", CMr1[i*3  ], CMr1[i*3+1], CMr1[i*3+2]);
//			printf("\n Y %e %e %e ", CMr2[i*3  ], CMr2[i*3+1], CMr2[i*3+2]);
//			printf("\n Z %e %e %e \n", CMr3[i*3  ], CMr3[i*3+1], CMr3[i*3+2]);
//		}
	}}
}

// Pressure gradient
void PrsGrdTrm_omp() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
//	if(Typ[i] == FLD) {
	if(Typ[i] != GST) {
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double pre_min = Prs[i];
		double Pi = Prs[i];
		double ni = pndi[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		if(GRD_TYP == 0 || GRD_TYP == 2) {
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					if(pre_min > Prs[j]) pre_min = Prs[j];
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}}
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wS = WEI_GRAD(dst, reS);
					if(GRD_TYP == 0)
						wS *= (Prs[j] - pre_min)/dstij2;
					else if(GRD_TYP == 1)
						wS *= (Prs[j] + Pi)/dstij2;
					else if(GRD_TYP == 2)
						wS *= (Prs[j] + Pi - 2*pre_min)/dstij2;
					else if(GRD_TYP == 3) {
						double nj = pndi[j];
						if(ni > 0.0001 && nj > 0.0001)
							wS *= (ni*Prs[j]/nj + nj*Pi/ni)/dstij2;
					}
					Acc_x += v0ij*wS;	Acc_y += v1ij*wS;	Acc_z += v2ij*wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		// A3 is a negative cte (-DIM/noGrad)
		// Original
//		Acc[i*3  ]=RLX_PRS*Acc_x*invDns[FLD]*A3;
//		Acc[i*3+1]=RLX_PRS*Acc_y*invDns[FLD]*A3;
//		Acc[i*3+2]=RLX_PRS*Acc_z*invDns[FLD]*A3;
		// Modified
		if(GRD_COR == 0) {
			Acc[i*3  ]=RLX_PRS*Acc_x*A3/RHO[i];
			Acc[i*3+1]=RLX_PRS*Acc_y*A3/RHO[i];
			Acc[i*3+2]=RLX_PRS*Acc_z*A3/RHO[i];
		}
		else {
		//	if(CMr1[1*3] > 1.0) {
		//		printf("\n X %e %e %e ", CMr1[i*3  ], CMr1[i*3+1], CMr1[i*3+2]);
		//		printf("\n Y %e %e %e ", CMr2[i*3  ], CMr2[i*3+1], CMr2[i*3+2]);
		//		printf("\n Z %e %e %e \n", CMr3[i*3  ], CMr3[i*3+1], CMr3[i*3+2]);
			//}
			Acc[i*3  ]=(RLX_PRS*A3/RHO[i])*(Acc_x*CMr1[i*3] + Acc_y*CMr1[i*3+1] + Acc_z*CMr1[i*3+2]);
			Acc[i*3+1]=(RLX_PRS*A3/RHO[i])*(Acc_x*CMr2[i*3] + Acc_y*CMr2[i*3+1] + Acc_z*CMr2[i*3+2]);
			Acc[i*3+2]=(RLX_PRS*A3/RHO[i])*(Acc_x*CMr3[i*3] + Acc_y*CMr3[i*3+1] + Acc_z*CMr3[i*3+2]);
		}
	}}
}

// Pressure gradient (Polygon wall)
void WallPrsGrdTrm_omp() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {

	//Nw[i]=1; // Only to show particles near polygon
	
	//if(Typ[i] == FLD) {
	if(Typ[i] == FLD && Nw[i] == 1) {
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double pre_min = Prs[i];
		double Pi = Prs[i];
		double ni = pndi[i];
		// Wall gradient Mitsume`s model
	    double Rref_i[9], normaliw[3], normaliwSqrt;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);

	    if(normaliwSqrt > 0.00000001) {
	    	normaliw[0] = normaliw[0]/normaliwSqrt;
	    	normaliw[1] = normaliw[1]/normaliwSqrt;
	    	normaliw[2] = normaliw[2]/normaliwSqrt;
	    }
	    else {
	    	normaliw[0] = 0;
	    	normaliw[1] = 0;
	    	normaliw[2] = 0;
	    }

	    //  Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
	    Rref_i[0] = 1.0 - 2*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2*normaliw[0]*normaliw[2];
		Rref_i[3] = 0.0 - 2*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2*normaliw[1]*normaliw[2];
		Rref_i[6] = 0.0 - 2*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2*normaliw[2]*normaliw[2];

		// Taylor pressure Pj
		double Rai[3];
		Rai[0] = Rref_i[0]*AccStar[i*3] + Rref_i[1]*AccStar[i*3+1] + Rref_i[2]*AccStar[i*3+2];
		Rai[1] = Rref_i[3]*AccStar[i*3] + Rref_i[4]*AccStar[i*3+1] + Rref_i[5]*AccStar[i*3+2];
		Rai[2] = Rref_i[6]*AccStar[i*3] + Rref_i[7]*AccStar[i*3+1] + Rref_i[8]*AccStar[i*3+2];

		// if(i == 16107)
		// {
		// 	printf("\ni:%5d TIM: %lf / Rref_i: ", i, TIM);
		// 	for(int rr=0; rr<9; rr++)
		// 		printf("%lf ", Rref_i[rr]);
		// 	printf("\ni:%5d TIM: %lf / Rai: ", i, TIM);
		// 	for(int rr=0; rr<3; rr++)
		// 		printf("%lf ", Rai[rr]);
		// 	printf("\ni:%5d TIM: %lf / ai: %lf, %lf, %lf", i, TIM, AccStar[i*3], AccStar[i*3+1], AccStar[i*3+2]);
			
		// }
			
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		if(GRD_TYP == 0 || GRD_TYP == 2) {
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;

				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					if(pre_min > Prs[j])pre_min = Prs[j];
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}}
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {

					double dst = sqrt(dstimj2);
					double wS = WEI_GRAD(dst, reS);

					// Taylor pressure Pj
					double Pj = Pi + RHO[i]*(Rai[0]*v0imj + Rai[1]*v1imj + Rai[2]*v2imj);
					Pj = Prs[j];

					// if(i == 16107)
					// 	printf("\ni:%5d j:%5d TIM: %lf / Pj: %lf / Pi: %lf / Zj: %lf / Zi: %lf", i, j, TIM, Pj, Pi, Pos[j*3+2], pos_miz);

					if(GRD_TYP == 0)
						wS *= (Pj - pre_min)/dstimj2;//(Prs[j] - pre_min)/dstimj2
					else if(GRD_TYP == 1)
						wS *= (Pj + Pi)/dstimj2;//(Prs[j] + Prs[i])/dstimj2
					else if(GRD_TYP == 2)
						wS *= (Pj + Pi - 2*pre_min)/dstimj2;//(Prs[j] + Prs[i] - 2*pre_min)/dstimj2
					else if(GRD_TYP == 3) {
						double nj = pndi[j];
						if(ni > 0.0001 && nj > 0.0001)
							wS *= (ni*Pj/nj + nj*Pi/ni)/dstimj2;
					}
					Acc_x += v0imj*wS;	Acc_y += v1imj*wS;	Acc_z += v2imj*wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		// Add "i" contribution ("i" is a neighboor of "mirror i")
	  	double v0imi = pos_ix - pos_mix;
		double v1imi = pos_iy - pos_miy;
		double v2imi = pos_iz - pos_miz;
		double dstimi2 = v0imi*v0imi+v1imi*v1imi+v2imi*v2imi;
		if(dstimi2 < reS2) {

			double dst = sqrt(dstimi2);
			double wS = WEI_GRAD(dst, reS);

			// Taylor pressure Pj
			double Pj = Pi + RHO[i]*(Rai[0]*v0imi + Rai[1]*v1imi + Rai[2]*v2imi);
			Pj = Prs[i];

			//if(i == 16107)
			//	printf("\ni:%5d TIM: %lf / Pj: %lf / Pi: %lf / Zj: %lf / Zi: %lf", i, TIM, Pj, Pi, pos_iz, pos_miz);

			if(GRD_TYP == 0)
				wS *= (Pj - pre_min)/dstimi2;//(Prs[i] - pre_min)/dstimi2
			else if(GRD_TYP == 1)
				wS *= (Pj + Pi)/dstimi2;//(Prs[i] + Prs[i])/dstimi2;
			else if(GRD_TYP == 2)
				wS *= (Pj + Pi - 2*pre_min)/dstimi2;//(Prs[i] + Prs[i] - 2*pre_min)/dstimi2;
			else if(GRD_TYP == 3) {
				double nj = pndi[i];
				if(ni > 0.0001 && nj > 0.0001)
					wS *= (ni*Pj/nj + nj*Pi/ni)/dstimi2;
			}
			Acc_x += v0imi*wS;	Acc_y += v1imi*wS;	Acc_z += v2imi*wS;
	  	}

		// Repulsive force
	  	double rpsForce[3];
	  	rpsForce[0]=rpsForce[1]=rpsForce[2] = 0.0;

	  	if(REP_FOR == 0) {
			// Parallel analysis system for free-surface flow using MPS method with explicitly represented polygon wall boundary model
			// https://doi.org/10.1007/s40571-019-00269-6
			if(normaliwSqrt < re_rep && normaliwSqrt != 0) {
				double wijRep = RHO[i]/(DT*DT)*(re_rep-normaliwSqrt);
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
			}
		}
		else if(REP_FOR == 1) {
			// Explicitly represented polygon wall boundary model for the explicit MPS method
			// https://doi.org/10.1007/s40571-015-0037-8
			if(normaliwSqrt < re_rep && normaliwSqrt != 0) {
				double wijRep = MIT_REP*WEI_GRAD(normaliwSqrt, re_rep);
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
			}
		}
		else if(REP_FOR == 2) {
			// Simulating Free Surface Flows with SPH
			// https://doi.org/10.1006/jcph.1994.1034
			if(normaliwSqrt < re_rep && normaliwSqrt != 0) {
				double maxVel2 = vMax*vMax;
				double R1 = (re_rep/normaliwSqrt)*(re_rep/normaliwSqrt);
				double R2 = R1*R1;
				double wijRep = (LNJ_REP*maxVel2/normaliwSqrt)*(R2-R1)*RHO[i];
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
			}
		}
		else {
			// SPH particle boundary forces for arbitrary boundaries 
			// https://doi.org/10.1016/j.cpc.2009.05.008
			if(normaliwSqrt < re_rep && normaliwSqrt != 0) {
				double maxVel2 = vMax*vMax;
				double W1 = (1+3/2*normaliwSqrt/(re_rep));
				double W2 = (1-normaliwSqrt/(re_rep))*(1-normaliwSqrt/(re_rep))*(1-normaliwSqrt/(re_rep));
				double wijRep = (KJT_REP*maxVel2/(normaliwSqrt - 0.0*PCL_DST))*(1/8)*(W1)*(W2)*RHO[i];
				rpsForce[0] = - wijRep*normaliw[0];
				rpsForce[1] = - wijRep*normaliw[1];
				rpsForce[2] = - wijRep*normaliw[2];
			}
	  	}
		// A3 is a negative cte (-DIM/noGrad)
		// Original
//		Acc[i*3  ] += (RLX_PRS*(Rref_i[0]*Acc_x + Rref_i[1]*Acc_y + Rref_i[2]*Acc_z)*A3 - rpsForce[0])*invDns[FLD];
//		Acc[i*3+1] += (RLX_PRS*(Rref_i[3]*Acc_x + Rref_i[4]*Acc_y + Rref_i[5]*Acc_z)*A3 - rpsForce[1])*invDns[FLD];
//		Acc[i*3+2] += (RLX_PRS*(Rref_i[6]*Acc_x + Rref_i[7]*Acc_y + Rref_i[8]*Acc_z)*A3 - rpsForce[2])*invDns[FLD];
		// Modified
		Acc[i*3  ] += (RLX_PRS*(Rref_i[0]*Acc_x + Rref_i[1]*Acc_y + Rref_i[2]*Acc_z)*A3 - rpsForce[0])/RHO[i];
		Acc[i*3+1] += (RLX_PRS*(Rref_i[3]*Acc_x + Rref_i[4]*Acc_y + Rref_i[5]*Acc_z)*A3 - rpsForce[1])/RHO[i];
		Acc[i*3+2] += (RLX_PRS*(Rref_i[6]*Acc_x + Rref_i[7]*Acc_y + Rref_i[8]*Acc_z)*A3 - rpsForce[2])/RHO[i];

	}}
}



// Calculation of the volume of fraction if phase II in the mixture
void VolFract_omp()
{
	if(Fraction_method == 1) {   //Linear distribution
#pragma omp parallel for schedule(dynamic,64)
		for(int i=0; i<nP; i++) {
		if(Typ[i] == FLD) {
			double sum1 = 0, sum2 = 0;
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
			int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
			int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
			for(int jz=iz-1;jz<=iz+1;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*nBxy + jy*nBx + jx;
				int j = bfst[jb];
				if(j == -1) continue;
				for(;;) {
					double v0 = Pos[j*3  ] - pos_ix;
					double v1 = Pos[j*3+1] - pos_iy;
					double v2 = Pos[j*3+2] - pos_iz;
					double dst2 = v0*v0+v1*v1+v2*v2;
					if(dst2 < reS2) {
					if(j != i && Typ[j] == FLD) {
						sum1 = sum1 + 1;
						if(PTYPE[j] >= 2) sum2 = sum2 + 1;
						}}
					j = nxt[j];
					if(j == -1) break;
				}
			}}}
			if(sum1 == 0)
				Cv[i] = 0.0;
			else 
				Cv[i] = sum2 / sum1;
		}}
	}
	else if(Fraction_method == 2) {   //Non linear :  Smoothed using the weight funtion
#pragma omp parallel for schedule(dynamic,64)
		for(int i=0; i<nP; i++) {
		if(Typ[i] == FLD) {
			double sum1 = 0, sum2 = 0;
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
			int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
			int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
			for(int jz=iz-1;jz<=iz+1;jz++) {
			for(int jy=iy-1;jy<=iy+1;jy++) {
			for(int jx=ix-1;jx<=ix+1;jx++) {
				int jb = jz*nBxy + jy*nBx + jx;
				int j = bfst[jb];
				if(j == -1) continue;
				for(;;) {
					double v0 = Pos[j*3  ] - pos_ix;
					double v1 = Pos[j*3+1] - pos_iy;
					double v2 = Pos[j*3+2] - pos_iz;
					double dst2 = v0*v0+v1*v1+v2*v2;
					if(dst2 < reS2) {
					if(j != i && Typ[j] == FLD) {
						double dst = sqrt(dst2);
						double wS = WEI(dst, reS);
						sum1 = sum1 + wS;
						if(PTYPE[j] >= 2) sum2 = sum2 + wS;
						}}
					j = nxt[j];
					if(j == -1) break;
				}
			}}}
			if(sum1 == 0)
				Cv[i] = 0.0;
			else 
				Cv[i] = sum2 / sum1;
		}}
	}
}

//void NonNwtVscTrm_omp(double *x_vel, double *y_vel, double *z_vel) {
// Viscosity interaction values for "real" fluid particles
void VscIntVal_omp() {

	//double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;//,grain_VF, *p_smooth;
	double d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
//	Xmin = 0.0 - PCL_DST*3; Xmax = 1.65 + PCL_DST*3;
//	Ymin = 0.0 - PCL_DST*3; Ymax = 0.15 + PCL_DST*3;
	Zmin = 0.0 - PCL_DST*3; //Zmax = 0.7 + PCL_DST*30;
	// damErosion3D
	Xmin = 0.0 - PCL_DST*3; Xmax = 2.00 + PCL_DST*3;
	Ymin = 0.0 - PCL_DST*3; Ymax = 0.10 + PCL_DST*3;
	// Changed !!!

	// Search free-surface particles for each interval of aa = 2 particles in wall
	int aa = 2, kx, ky;
	int kx_max = int((Xmax - Xmin) / aa / PCL_DST) + 1;
	int ky_max = int((Ymax - Ymin) / aa / PCL_DST) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
/*
	S11 = new double[nP + 1];
	S22 = new double[nP + 1];
	S33 = new double[nP + 1];
	S12 = new double[nP + 1];
	S13 = new double[nP + 1];
	S23 = new double[nP + 1];
*/
	//p_smooth = new double[nP + 1];
	BL = new double*[kx_max + 1];  // bed level
	WL = new double*[kx_max + 1];  // water level
	PS = new double*[kx_max + 1];  // pressure sediment

#pragma omp parallel for
	for(int m = 1; m <= kx_max; m++)
	{
		BL[m] = new double[ky_max + 1];
		WL[m] = new double[ky_max + 1];
		PS[m] = new double[ky_max + 1];
	}

	// Determining the bed level
#pragma omp parallel for schedule(dynamic,64)
	for(kx = 1; kx <= kx_max; kx++)
	{
		for(ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0; i<nP; i++) {
		if(Typ[i] == FLD) {
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		
			kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
			ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;


			//if(pos_iz > BL[kx][ky] && Cv[i] > 0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = pnew[i]; }
			if(pos_iz > BL[kx][ky] && Cv[i] > 0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = Prs[i]; }
			if(pos_iz > WL[kx][ky] && PTYPE[i] == 1) { WL[kx][ky] = pos_iz; }
		}
	}

	// Strain rate calculation
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] == FLD) {
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;

//		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wS = WEI(dst, reS);
					double vec_ijx = Vel[j*3  ]-vec_ix;	
					double vec_ijy = Vel[j*3+1]-vec_iy;	
					double vec_ijz = Vel[j*3+2]-vec_iz;

					sum1 += vec_ijx*v0ij*wS/dstij2;
					sum2 += vec_ijx*v1ij*wS/dstij2;
					sum3 += vec_ijx*v2ij*wS/dstij2;
					
					sum4 += vec_ijy*v0ij*wS/dstij2;
					sum5 += vec_ijy*v1ij*wS/dstij2;
					sum6 += vec_ijy*v2ij*wS/dstij2;
					
					sum7 += vec_ijz*v0ij*wS/dstij2;
					sum8 += vec_ijz*v1ij*wS/dstij2;
					sum9 += vec_ijz*v2ij*wS/dstij2;
					
					sum10 += Prs[j]*wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

		// A3 is a negative cte (-DIM/n0Grad)
		Uxx = -A3*sum1; Uxy = -A3*sum2; Uxz = -A3*sum3;
		Uyx = -A3*sum4; Uyy = -A3*sum5; Uyz = -A3*sum6;
		Uzx = -A3*sum7; Uzy = -A3*sum8; Uzz = -A3*sum9;

		p_smooth[i] = sum10 / n0Grad;
		if(p_smooth[i] < 0) p_smooth[i] = 0;

		S11[i] = 0.5*(Uxx + Uxx);
		S12[i] = 0.5*(Uxy + Uyx);
		S13[i] = 0.5*(Uxz + Uzx);
		S22[i] = 0.5*(Uyy + Uyy);
		S23[i] = 0.5*(Uyz + Uzy);
		S33[i] = 0.5*(Uzz + Uzz);

		//II[i] = 0.5*Uxx*Uxx + 0.5*Uyy*Uyy + 0.25*(Uxy + Uyx)*(Uxy + Uyx);
		//II[i] = 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		II[i] = 0.5*(S11[i] * S11[i] + 2 * S12[i] * S12[i] + 2 * S13[i] * S13[i] + S22[i] * S22[i] + 2 * S23[i] * S23[i] + S33[i] * S33[i]);
		//II[i]= S11[i]*S22[i] +S22[i]*S33[i]+ S11[i]*S33[i] - S12[i]*S12[i] -S13[i]*S13[i]- S23[i]*S23[i];
		if(II[i] < 0 || II[i] * 0 != 0) II[i] = 0;
		//II=fabs(S11[i]*S22[i]-S12[i]*S12[i]);
		
		//std::cout << " II: " << II[i] << std::endl;
	}}

	// Newtonian viscosity
	if(Fluid2_type == 0)
	{
#pragma omp parallel for
		for(int i=0; i<nP; i++) {
		if(Typ[i] == FLD) {
			if(PTYPE[i] <= 1)MEU[i] = KNM_VS1 * DNS_FL1;
			if(PTYPE[i] != 1)MEU[i] = KNM_VS2 * DNS_FL2;
		}}

//		if(TURB > 0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2*sqrt(II[i]);

//			if(NEUt[i] * 0 != 0)  NEUt[i] = 0;
//			if(NEUt[i] > 1)     NEUt[i] = 1;
//		}
	}

	// Granular Fluid
	if(Fluid2_type == 1)
	{
	// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
//#pragma omp parallel for
		for(int i=0; i<nP; i++) {
		if(Typ[i] == FLD) {

//			double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];

			if(PTYPE[i] == 1) {
				MEU[i] = KNM_VS1 * DNS_FL1;
			}
			else if(PTYPE[i] == 2) {
				// phi: internal friction angle
				// phi2: maximum friction angle
				phi = (Cv[i] - 0.25)*PHI / (1 - 0.25);
				phi2 = (Cv[i] - 0.25)*PHI_2 / (1 - 0.25);
				if(Cv[i] <= 0.25) { phi = 0.00001; phi2 = 0.00001; } // phi close to zero
				if(PTYPE[i] <= 0) phi = PHI_BED; // ghost

				// Normal stress calculation (mechanical pressure)
				p_rheo_new[i] = p_smooth[i];

				kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
				ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress=(BL[k]-pos_iy+DL/2)*(DNS_FL2)*9.81;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - pos_iz + PCL_DST / 2)*(DNS_FL2 - DNS_FL1)*9.81 - (vec_ix*vec_ix + vec_iy*vec_iy + vec_iz*vec_iz)*(DNS_FL2 - DNS_FL1) / 2.0;	// normal_stress= Gama.H

				if(p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81<0) p_smooth[i] = (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
				if(TIM <= 1) normal_stress = 1.0*(1 - TIM)*(p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81) + 1.0*(TIM)*normal_stress;

//				normal_stress = p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
//				normal_stress = p_smooth[i];
				//normal_stress=normal_stress*0.61*1500/DNS_FL2;
				if(normal_stress < 1 || Cv[i] < 0.5) normal_stress = 1;

				p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/DNS_SDT);		// Free-fall (dry granular material)
				Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/(DNS_FL1*Cd));	// Grain inertia (submerged)
				//Inertia[i] = sqrt(II[i])*(KNM_VS1*DNS_FL1)/normal_stress ;	// Viscous regime

//				Inertia[i] = 1.0;
				// VF_max VF_min
				VF[i] = VF_max - (VF_max - VF_min)*Inertia[i];
				if(VF[i] < VF_min) VF[i] = VF_min;
				RHO[i] = DNS_SDT * VF[i] + (1-VF[i])*DNS_FL1;
				phi = phi * VF[i] / VF_max;

				// Mohr-Coulomb
				double yield_stress = cohes * cos(phi) + normal_stress * sin(phi);

				if(yield_stress < 0) yield_stress = 0;

				double visc_max = (yield_stress*mm*0.5 + MEU0);

				if(II[i] > 0)
					MEU_Y[i] = yield_stress * (1 - exp(-mm * sqrt(II[i]))) / 2.0 / sqrt(II[i]);
				else
					MEU_Y[i] = visc_max;

				// H-B rheology

				//meu_0 = MEU0;

				// Non-linear Meu(I) rheology
				//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/DNS_FL2)+sqrt(II[i])*DG);			//free fall
				meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/(DNS_FL1*Cd))+sqrt(II[i])*DG);		//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(KNM_VS1*DNS_FL1)/(I0*normal_stress+sqrt(II[i])*(KNM_VS1*DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL2)/I0;		//free fall
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL1*Cd)/I0;	//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*(KNM_VS1*DNS_FL1)/I0;					//viscous

				if(II[i] <= 0 || (meu_0 * 0) != 0) meu_0 = MEU0;

				visc_max = (yield_stress*mm*0.5 + meu_0);

				//if(isnan(II[i]) || isinf(II[i])) {
					//std::cout << " viscmax: " << II[i] << std::endl;
				//	assert(visc_max >= 0 || visc_max <= 0);
				//}
				
				// Herschel bulkley papanastasiou
				MEU[i] = MEU_Y[i] + MEU0 * pow(4 * II[i], (N - 1) / 2);

				// MEU_Y rheological model
				//MEU[i] = MEU_Y[i] + meu_0;
				
				if(II[i] == 0 || MEU[i]>visc_max) {
					//std::cout << " MEU>viscmax: " << yield_stress*mm*0.5 << " meu0: " << meu_0 << " II: " << II[i] << std::endl;
					MEU[i] = visc_max;
				}
				if(PTYPE[i] <= 0) MEU[i] = MEU[i] * Cv[i] + DNS_FL1*KNM_VS1*(1 - Cv[i]);

				//if(MEU[i]/RHO[i] > maxVIS) maxVIS = MEU[i]/RHO[i];
			}
			
			if(PTYPE[i] >= 2) {
				if(Cv[i] > 0.5) RHO[i] = DNS_FL2;
				else RHO[i] = Cv[i] * DNS_FL2 + (1 - Cv[i]) * DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if(stress_cal_method == 2)
//		{
//			for(i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
//				for(l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if(i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2 * MEU[i] * MEU[j] / (MEU[i] + MEU[j]);
//						if((NEUt[i] + NEUt[j])>0) meuij = meuij + 2 * NEUt[i] * RHO[i] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j]);

//						sum1 = sum1 + meuij * (x_vel[j] - x_vel[i])*DX(i, j)*w / d / d;
//						sum2 = sum2 + meuij * (x_vel[j] - x_vel[i])*DY(i, j)*w / d / d;
//						sum3 = sum3 + meuij * (x_vel[j] - x_vel[i])*DZ(i, j)*w / d / d;

//						sum4 = sum4 + meuij * (y_vel[j] - y_vel[i])*DX(i, j)*w / d / d;
//						sum5 = sum5 + meuij * (y_vel[j] - y_vel[i])*DY(i, j)*w / d / d;
//						sum6 = sum6 + meuij * (y_vel[j] - y_vel[i])*DZ(i, j)*w / d / d;

//						sum7 = sum7 + meuij * (z_vel[j] - z_vel[i])*DX(i, j)*w / d / d;
//						sum8 = sum8 + meuij * (z_vel[j] - z_vel[i])*DY(i, j)*w / d / d;
//						sum9 = sum9 + meuij * (z_vel[j] - z_vel[i])*DZ(i, j)*w / d / d;
//					}
//				}

//				Tau_xx[i] = (DIM / n0) * 2 * sum1;
//				Tau_yy[i] = (DIM / n0) * 2 * sum5;
//				Tau_zz[i] = (DIM / n0) * 2 * sum9;

//				Tau_xy[i] = (DIM / n0)*(sum2 + sum4);
//				Tau_xz[i] = (DIM / n0)*(sum3 + sum7);
//				Tau_yz[i] = (DIM / n0)*(sum6 + sum8);
//			}
//		}

	} // if(Fluid2_type == 1)

	//---------------------------------------------------------------

//	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS; //delete[]p_smooth;
//	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL; //p_smooth = NULL;

	delete[]BL; delete[]WL; delete[]PS;
	BL = NULL; WL = NULL; PS = NULL;
}

// Free-slip condition. Viscosity interaction values
void WallSlipVscIntVal_omp() {

	//double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;//, grain_VF, *p_smooth;
	double d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
//	Xmin = 0.0 - PCL_DST*3; Xmax = 1.65 + PCL_DST*3;
//	Ymin = 0.0 - PCL_DST*3; Ymax = 0.15 + PCL_DST*3;
	Zmin = 0.0 - PCL_DST*3; //Zmax = 0.7 + PCL_DST*30;
	// damErosion3D
	Xmin = 0.0 - PCL_DST*3; Xmax = 2.00 + PCL_DST*3;
	Ymin = 0.0 - PCL_DST*3; Ymax = 0.10 + PCL_DST*3;
	// Changed !!!

	// Search free-surface particles for each interval of aa = 2 particles in wall
	int aa = 2, kx, ky;
	int kx_max = int((Xmax - Xmin) / aa / PCL_DST) + 1;
	int ky_max = int((Ymax - Ymin) / aa / PCL_DST) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
	double aUxx, aUxy, aUxz, aUyx, aUyy, aUyz, aUzx, aUzy, aUzz;
/*
	S11 = new double[nP + 1];
	S22 = new double[nP + 1];
	S33 = new double[nP + 1];
	S12 = new double[nP + 1];
	S13 = new double[nP + 1];
	S23 = new double[nP + 1];
*/
	//p_smooth = new double[nP + 1];
	BL = new double*[kx_max + 1];  // bed level
	WL = new double*[kx_max + 1];  // water level
	PS = new double*[kx_max + 1];  // pressure sediment

#pragma omp parallel for
	for(int m = 1; m <= kx_max; m++)
	{
		BL[m] = new double[ky_max + 1];
		WL[m] = new double[ky_max + 1];
		PS[m] = new double[ky_max + 1];
	}

	// Determining the bed level
#pragma omp parallel for schedule(dynamic,64)
	for(kx = 1; kx <= kx_max; kx++)
	{
		for(ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0; i<nP; i++) {
		if(Typ[i] == FLD) {
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			
			kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
			ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;

			//if(pos_iz>BL[kx][ky] && Cv[i]>0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = pnew[i]; }
			if(pos_iz>BL[kx][ky] && Cv[i]>0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = Prs[i]; }
			if(pos_iz>WL[kx][ky] && PTYPE[i] == 1) { WL[kx][ky] = pos_iz; }
		}
	}

	// Strain rate calculation
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
	//if(Typ[i] == FLD) {
	if(Typ[i] == FLD && Nw[i] == 1) {
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;

//		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
//		double vec_ix = Velk[i*3  ];	double vec_iy = Velk[i*3+1];	double vec_iz = Velk[i*3+2

		// Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
		double Rref_i[9], normaliw[3], normalMod2;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

	    if(normalMod2 > 0.00000001) {
	    	double normalMod = sqrt(normalMod2);
	    	normaliw[0] = normaliw[0]/normalMod;
	    	normaliw[1] = normaliw[1]/normalMod;
	    	normaliw[2] = normaliw[2]/normalMod;
	    }
	    else {
	    	normaliw[0] = 0;
	    	normaliw[1] = 0;
	    	normaliw[2] = 0;
	    }

	    //  Transformation matrix R_i = I - 2*normal_iwall*normal_iwall
	    Rref_i[0] = 1.0 - 2*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2*normaliw[0]*normaliw[2];
		Rref_i[3] = 0.0 - 2*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2*normaliw[1]*normaliw[2];
		Rref_i[6] = 0.0 - 2*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2*normaliw[2]*normaliw[2];

		// Mirror particle velocity vi' = Ri * vi
      	double vec_mix = (Rref_i[0]*vec_ix + Rref_i[1]*vec_iy + Rref_i[2]*vec_iz);
		double vec_miy = (Rref_i[3]*vec_ix + Rref_i[4]*vec_iy + Rref_i[5]*vec_iz);
		double vec_miz = (Rref_i[6]*vec_ix + Rref_i[7]*vec_iy + Rref_i[8]*vec_iz);

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstimj2);
					double wS = WEI(dst, reS);

					double vec_mijx = Vel[j*3  ]-vec_mix;	
					double vec_mijy = Vel[j*3+1]-vec_miy;	
					double vec_mijz = Vel[j*3+2]-vec_miz;

					sum1 += vec_mijx*v0imj*wS/dstimj2;
					sum2 += vec_mijx*v1imj*wS/dstimj2;
					sum3 += vec_mijx*v2imj*wS/dstimj2;
					
					sum4 += vec_mijy*v0imj*wS/dstimj2;
					sum5 += vec_mijy*v1imj*wS/dstimj2;
					sum6 += vec_mijy*v2imj*wS/dstimj2;
					
					sum7 += vec_mijz*v0imj*wS/dstimj2;
					sum8 += vec_mijz*v1imj*wS/dstimj2;
					sum9 += vec_mijz*v2imj*wS/dstimj2;
					
					sum10 += Prs[j]*wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

		// Rref_i * gradU
		// A3 is a negative cte (-DIM/n0Grad)
		aUxx = -A3*(Rref_i[0]*sum1 + Rref_i[1]*sum4 + Rref_i[2]*sum7);
		aUxy = -A3*(Rref_i[0]*sum2 + Rref_i[1]*sum5 + Rref_i[2]*sum8);
		aUxz = -A3*(Rref_i[0]*sum3 + Rref_i[1]*sum6 + Rref_i[2]*sum9);

		aUyx = -A3*(Rref_i[3]*sum1 + Rref_i[4]*sum4 + Rref_i[5]*sum7);
		aUyy = -A3*(Rref_i[3]*sum2 + Rref_i[4]*sum5 + Rref_i[5]*sum8);
		aUyz = -A3*(Rref_i[3]*sum3 + Rref_i[4]*sum6 + Rref_i[5]*sum9);

		aUzx = -A3*(Rref_i[6]*sum1 + Rref_i[7]*sum4 + Rref_i[8]*sum7);
		aUzy = -A3*(Rref_i[6]*sum2 + Rref_i[7]*sum5 + Rref_i[8]*sum8);
		aUzz = -A3*(Rref_i[6]*sum3 + Rref_i[7]*sum6 + Rref_i[8]*sum9);

		// Rref_i * gradU * Rref_i
		Uxx = aUxx*Rref_i[0] + aUxy*Rref_i[3] + aUxz*Rref_i[6];
		Uxy = aUxx*Rref_i[1] + aUxy*Rref_i[4] + aUxz*Rref_i[7];
		Uxz = aUxx*Rref_i[2] + aUxy*Rref_i[5] + aUxz*Rref_i[8];

		Uyx = aUyx*Rref_i[0] + aUyy*Rref_i[3] + aUyz*Rref_i[6];
		Uyy = aUyx*Rref_i[1] + aUyy*Rref_i[4] + aUyz*Rref_i[7];
		Uyz = aUyx*Rref_i[2] + aUyy*Rref_i[5] + aUyz*Rref_i[8];

		Uzx = aUzx*Rref_i[0] + aUzy*Rref_i[3] + aUzz*Rref_i[6];
		Uzy = aUzx*Rref_i[1] + aUzy*Rref_i[4] + aUzz*Rref_i[7];
		Uzz = aUzx*Rref_i[2] + aUzy*Rref_i[5] + aUzz*Rref_i[8];

		// Addition of smoothed pressure for particles near mesh
		p_smooth[i] += sum10 / n0Grad;
		if(p_smooth[i]<0) p_smooth[i] = 0;

		// 0.5*(gradU + gradUt)
		S11[i] = 0.5*(Uxx + Uxx);
		S12[i] = 0.5*(Uxy + Uyx);
		S13[i] = 0.5*(Uxz + Uzx);
		S22[i] = 0.5*(Uyy + Uyy);
		S23[i] = 0.5*(Uyz + Uzy);
		S33[i] = 0.5*(Uzz + Uzz);

		//II[i] = 0.5*Uxx*Uxx + 0.5*Uyy*Uyy + 0.25*(Uxy + Uyx)*(Uxy + Uyx);

		// Addition of II for particles near mesh
		//II[i] += 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		II[i] += 0.5*(S11[i] * S11[i] + 2 * S12[i] * S12[i] + 2 * S13[i] * S13[i] + S22[i] * S22[i] + 2 * S23[i] * S23[i] + S33[i] * S33[i]);
//		II[i] = 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		//II[i]= S11[i]*S22[i] +S22[i]*S33[i]+ S11[i]*S33[i] - S12[i]*S12[i] -S13[i]*S13[i]- S23[i]*S23[i] ;
		if(II[i]<0 || II[i] * 0 != 0) II[i] = 0;
		//II=fabs(S11[i]*S22[i]-S12[i]*S12[i]);
	}}
	
	// Newtonian viscosity
	if(Fluid2_type == 0)
	{
	// Loop only for particles near mesh
#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<nP; i++) {
		//if(Typ[i] == FLD) {
		if(Typ[i] == FLD && Nw[i] == 1) {
			if(PTYPE[i] <= 1)MEU[i] = KNM_VS1 * DNS_FL1;
			if(PTYPE[i] != 1)MEU[i] = KNM_VS2 * DNS_FL2;
		}}

//		if(TURB>0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2*sqrt(II[i]);

//			if(NEUt[i] * 0 != 0)  NEUt[i] = 0;
//			if(NEUt[i]>1)     NEUt[i] = 1;
//		}
	}

	// Granular Fluid
	if(Fluid2_type == 1)
	{
		// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
		// Loop only for particles near mesh
//#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<nP; i++) {
		//if(Typ[i] == FLD) {
		if(Typ[i] == FLD && Nw[i] == 1) {
//			double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];

			if(PTYPE[i] == 1) {
				MEU[i] = KNM_VS1 * DNS_FL1;
			}
			else if(PTYPE[i] == 2) {
				phi = (Cv[i] - 0.25)*PHI / (1 - 0.25);
				phi2 = (Cv[i] - 0.25)*PHI_2 / (1 - 0.25);
				if(Cv[i] <= 0.25) { phi = 0.00001; phi2 = 0.00001; } // phi close to zero
				if(PTYPE[i] <= 0) phi = PHI_BED;

				// Normal stress calculation (mehcanical pressure)
				p_rheo_new[i] = p_smooth[i];

				kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
				ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress=(BL[k]-pos_iy+DL/2)*(DNS_FL2)*9.81;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - pos_iz + PCL_DST / 2)*(DNS_FL2 - DNS_FL1)*9.81 - (vec_ix*vec_ix + vec_iy*vec_iy + vec_iz*vec_iz)*(DNS_FL2 - DNS_FL1) / 2.0;	// normal_stress= Gama.H

				if(p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81<0) p_smooth[i] = (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
				if(TIM <= 1) normal_stress = 1.0*(1 - TIM)*(p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81) + 1.0*(TIM)*normal_stress;

//				normal_stress = p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
//				normal_stress = p_smooth[i];
				//normal_stress=normal_stress*0.61*1500/DNS_FL2;
				if(normal_stress < 1 || Cv[i] < 0.5) normal_stress = 1;

				p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/DNS_SDT);		// Free-fall (dry granular material)
				Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/(DNS_FL1*Cd));	// Grain inertia (submerged)
				//Inertia[i] = sqrt(II[i])*(KNM_VS1*DNS_FL1)/normal_stress ;	// Viscous regime

				// VF_max VF_min
				VF[i] = VF_max - (VF_max - VF_min)*Inertia[i];
				if(VF[i] < VF_min) VF[i] = VF_min;
				RHO[i] = DNS_SDT * VF[i] + (1-VF[i])*DNS_FL1;
				phi = phi * VF[i] / VF_max;

				double yield_stress = cohes * cos(phi) + normal_stress * sin(phi);

				if(yield_stress < 0) yield_stress = 0;

				double visc_max = (yield_stress*mm*0.5 + MEU0);

				if(II[i]>0)
					MEU_Y[i] = yield_stress * (1 - exp(-mm * sqrt(II[i]))) / 2.0 / sqrt(II[i]);
				else
					MEU_Y[i] = visc_max;

				// H-B rheology

				//meu_0 = MEU0;

				// Non-linear Meu(I) rheology
				//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/DNS_FL2)+sqrt(II[i])*DG);			//free fall
				meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/(DNS_FL1*Cd))+sqrt(II[i])*DG);		//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(KNM_VS1*DNS_FL1)/(I0*normal_stress+sqrt(II[i])*(KNM_VS1*DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL2)/I0;		//free fall
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL1*Cd)/I0;	//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*(KNM_VS1*DNS_FL1)/I0;					//viscous

				if(II[i] <= 0 || (meu_0 * 0) != 0) meu_0 = MEU0;

				visc_max = (yield_stress*mm*0.5 + meu_0);

				// Herschel bulkley papanastasiou
				MEU[i] = MEU_Y[i] + MEU0 * pow(4 * II[i], (N - 1) / 2);

				// MEU_Y rheological model
				//MEU[i] = MEU_Y[i] + meu_0;
				
				if(II[i] == 0 || MEU[i]>visc_max) MEU[i] = visc_max;
				if(PTYPE[i] <= 0) MEU[i] = MEU[i] * Cv[i] + DNS_FL1*KNM_VS1*(1 - Cv[i]);
			}
			
			if(PTYPE[i] >= 2) {
				if(Cv[i] > 0.5) RHO[i] = DNS_FL2;
				else RHO[i] = Cv[i] * DNS_FL2 + (1 - Cv[i]) * DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if(stress_cal_method == 2)
//		{
//			for(i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
//				for(l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if(i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2 * MEU[i] * MEU[j] / (MEU[i] + MEU[j]);
//						if((NEUt[i] + NEUt[j])>0) meuij = meuij + 2 * NEUt[i] * RHO[i] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j]);

//						sum1 = sum1 + meuij * (x_vel[j] - x_vel[i])*DX(i, j)*w / d / d;
//						sum2 = sum2 + meuij * (x_vel[j] - x_vel[i])*DY(i, j)*w / d / d;
//						sum3 = sum3 + meuij * (x_vel[j] - x_vel[i])*DZ(i, j)*w / d / d;

//						sum4 = sum4 + meuij * (y_vel[j] - y_vel[i])*DX(i, j)*w / d / d;
//						sum5 = sum5 + meuij * (y_vel[j] - y_vel[i])*DY(i, j)*w / d / d;
//						sum6 = sum6 + meuij * (y_vel[j] - y_vel[i])*DZ(i, j)*w / d / d;

//						sum7 = sum7 + meuij * (z_vel[j] - z_vel[i])*DX(i, j)*w / d / d;
//						sum8 = sum8 + meuij * (z_vel[j] - z_vel[i])*DY(i, j)*w / d / d;
//						sum9 = sum9 + meuij * (z_vel[j] - z_vel[i])*DZ(i, j)*w / d / d;
//					}
//				}

//				Tau_xx[i] = (DIM / n0) * 2 * sum1;
//				Tau_yy[i] = (DIM / n0) * 2 * sum5;
//				Tau_zz[i] = (DIM / n0) * 2 * sum9;

//				Tau_xy[i] = (DIM / n0)*(sum2 + sum4);
//				Tau_xz[i] = (DIM / n0)*(sum3 + sum7);
//				Tau_yz[i] = (DIM / n0)*(sum6 + sum8);
//			}
//		}

	} // if(Fluid2_type == 1)

	//---------------------------------------------------------------

//	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS;// delete[]p_smooth;
//	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL;// p_smooth = NULL;

	delete[]BL; delete[]WL; delete[]PS;
	BL = NULL; WL = NULL; PS = NULL;
}

// No-Slip condition. Viscosity interaction values
void WallNoSlipVscIntVal_omp() {

	//double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;//, grain_VF, *p_smooth;
	double d, phi = 0.0, phi2 = 0.0, meu_0, normal_stress;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
//	Xmin = 0.0 - PCL_DST*3; Xmax = 1.65 + PCL_DST*3;
//	Ymin = 0.0 - PCL_DST*3; Ymax = 0.15 + PCL_DST*3;
	Zmin = 0.0 - PCL_DST*3; //Zmax = 0.7 + PCL_DST*30;
	// damErosion3D
	Xmin = 0.0 - PCL_DST*3; Xmax = 2.00 + PCL_DST*3;
	Ymin = 0.0 - PCL_DST*3; Ymax = 0.10 + PCL_DST*3;
	// Changed !!!

	// Search free-surface particles for each interval of aa = 2 particles in wall
	int aa = 2, kx, ky;
	int kx_max = int((Xmax - Xmin) / aa / PCL_DST) + 1;
	int ky_max = int((Ymax - Ymin) / aa / PCL_DST) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
/*
	S11 = new double[nP + 1];
	S22 = new double[nP + 1];
	S33 = new double[nP + 1];
	S12 = new double[nP + 1];
	S13 = new double[nP + 1];
	S23 = new double[nP + 1];
*/
	//p_smooth = new double[nP + 1];
	BL = new double*[kx_max + 1];  // bed level
	WL = new double*[kx_max + 1];  // water level
	PS = new double*[kx_max + 1];  // pressure sediment

#pragma omp parallel for
	for(int m = 1; m <= kx_max; m++)
	{
		BL[m] = new double[ky_max + 1];
		WL[m] = new double[ky_max + 1];
		PS[m] = new double[ky_max + 1];
	}

	// Determining the bed level
#pragma omp parallel for schedule(dynamic,64)
	for(kx = 1; kx <= kx_max; kx++)
	{
		for(ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0; i<nP; i++) {
		if(Typ[i] == FLD) {
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		
			kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
			ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;

			//if(pos_iz>BL[kx][ky] && Cv[i]>0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = pnew[i]; }
			if(pos_iz>BL[kx][ky] && Cv[i]>0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = Prs[i]; }
			if(pos_iz>WL[kx][ky] && PTYPE[i] == 1) { WL[kx][ky] = pos_iz; }
		}
	}

	// Strain rate calculation
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
	//if(Typ[i] == FLD) {
	if(Typ[i] == FLD && Nw[i] == 1) {
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;

//		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
//		double vec_ix = Velk[i*3  ];	double vec_iy = Velk[i*3+1];	double vec_iz = Velk[i*3+2];

		// Transformation matrix R_i = I
		double Rref_i[9], Rinv_i[9], normaliw[3], normalMod2;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

	    if(normalMod2 > 0.00000001) {
	    	double normalMod = sqrt(normalMod2);
	    	normaliw[0] = normaliw[0]/normalMod;
	    	normaliw[1] = normaliw[1]/normalMod;
	    	normaliw[2] = normaliw[2]/normalMod;
	    }
	    else {
	    	normaliw[0] = 0;
	    	normaliw[1] = 0;
	    	normaliw[2] = 0;
	    }

	    //  Inverse transformation matrix Rinv_i = - I
	    Rinv_i[0] = -1.0; Rinv_i[1] =  0.0; Rinv_i[2] =  0.0;
		Rinv_i[3] =  0.0; Rinv_i[4] = -1.0; Rinv_i[5] =  0.0;
		Rinv_i[6] =  0.0; Rinv_i[7] =  0.0; Rinv_i[8] = -1.0;

		//  Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
	    Rref_i[0] = 1.0 - 2*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2*normaliw[0]*normaliw[2];
		Rref_i[3] = 0.0 - 2*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2*normaliw[1]*normaliw[2];
		Rref_i[6] = 0.0 - 2*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2*normaliw[2]*normaliw[2];

		double viwall[3], vtil[3];
		// Wall velocity (0 if fixed)
		viwall[0]=viwall[1]=viwall[2]=0.0;

#ifdef FORCED_ON
		if(meshID[i] == 2) {
			viwall[0] = velVWall[0];
			viwall[1] = velVWall[1];
			viwall[2] = velVWall[2];
		}
#endif

		// normal_iwall*v_iwall
		double dotnv = normaliw[0]*viwall[0] + normaliw[1]*viwall[1] + normaliw[2]*viwall[2];
		// vtil = vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}
		vtil[0] = vec_ix - 2*(viwall[0] - dotnv*normaliw[0]);
		vtil[1] = vec_iy - 2*(viwall[1] - dotnv*normaliw[1]);
		vtil[2] = vec_iz - 2*(viwall[2] - dotnv*normaliw[2]);
		// Mirror particle velocity vi' = Ri_inv * [vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}] 
      	double vec_mix = (Rinv_i[0]*vtil[0] + Rinv_i[1]*vtil[1] + Rinv_i[2]*vtil[2]);
		double vec_miy = (Rinv_i[3]*vtil[0] + Rinv_i[4]*vtil[1] + Rinv_i[5]*vtil[2]);
		double vec_miz = (Rinv_i[6]*vtil[0] + Rinv_i[7]*vtil[1] + Rinv_i[8]*vtil[2]);

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstimj2);
					double wS = WEI(dst, reS);

					double vec_mijx = Vel[j*3  ]-vec_mix;	
					double vec_mijy = Vel[j*3+1]-vec_miy;	
					double vec_mijz = Vel[j*3+2]-vec_miz;

					sum1 += vec_mijx*v0imj*wS/dstimj2;
					sum2 += vec_mijx*v1imj*wS/dstimj2;
					sum3 += vec_mijx*v2imj*wS/dstimj2;
					
					sum4 += vec_mijy*v0imj*wS/dstimj2;
					sum5 += vec_mijy*v1imj*wS/dstimj2;
					sum6 += vec_mijy*v2imj*wS/dstimj2;
					
					sum7 += vec_mijz*v0imj*wS/dstimj2;
					sum8 += vec_mijz*v1imj*wS/dstimj2;
					sum9 += vec_mijz*v2imj*wS/dstimj2;
					
					sum10 += Prs[j]*wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

		// Rinv_i * gradU * Rref_i = - gradU * Rref_i
		// A3 is a negative cte (-DIM/n0Grad)
		Uxx = A3*(sum1*Rref_i[0] + sum2*Rref_i[3] + sum3*Rref_i[6]);
		Uxy = A3*(sum1*Rref_i[1] + sum2*Rref_i[4] + sum3*Rref_i[7]);
		Uxz = A3*(sum1*Rref_i[2] + sum2*Rref_i[5] + sum3*Rref_i[8]);

		Uyx = A3*(sum4*Rref_i[0] + sum5*Rref_i[3] + sum6*Rref_i[6]);
		Uyy = A3*(sum4*Rref_i[1] + sum5*Rref_i[4] + sum6*Rref_i[7]);
		Uyz = A3*(sum4*Rref_i[2] + sum5*Rref_i[5] + sum6*Rref_i[8]);

		Uzx = A3*(sum7*Rref_i[0] + sum8*Rref_i[3] + sum9*Rref_i[6]);
		Uzy = A3*(sum7*Rref_i[1] + sum8*Rref_i[4] + sum9*Rref_i[7]);
		Uzz = A3*(sum7*Rref_i[2] + sum8*Rref_i[5] + sum9*Rref_i[8]);

		// Addition of smoothed pressure for particles near mesh
		p_smooth[i] += sum10 / n0Grad;
		if(p_smooth[i]<0) p_smooth[i] = 0;

		// - (Rref_i * gradU) - (Rref_i * gradU)t
		S11[i] = 0.5*(Uxx + Uxx);
		S12[i] = 0.5*(Uxy + Uyx);
		S13[i] = 0.5*(Uxz + Uzx);
		S22[i] = 0.5*(Uyy + Uyy);
		S23[i] = 0.5*(Uyz + Uzy);
		S33[i] = 0.5*(Uzz + Uzz);

		//II[i] = 0.5*Uxx*Uxx + 0.5*Uyy*Uyy + 0.25*(Uxy + Uyx)*(Uxy + Uyx);
		
		// Addition of II for particles near mesh
		//II[i] += 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		II[i] += 0.5*(S11[i] * S11[i] + 2 * S12[i] * S12[i] + 2 * S13[i] * S13[i] + S22[i] * S22[i] + 2 * S23[i] * S23[i] + S33[i] * S33[i]);
//		II[i] = 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		//II[i]= S11[i]*S22[i] +S22[i]*S33[i]+ S11[i]*S33[i] - S12[i]*S12[i] -S13[i]*S13[i]- S23[i]*S23[i] ;
		if(II[i]<0 || II[i] * 0 != 0) II[i] = 0;
		//II=fabs(S11[i]*S22[i]-S12[i]*S12[i]);
	}}
	
	// Newtonian viscosity
	if(Fluid2_type == 0)
	{
	// Loop only for particles near mesh
#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<nP; i++) {
		//if(Typ[i] == FLD) {
		if(Typ[i] == FLD && Nw[i] == 1) {
			if(PTYPE[i] <= 1)MEU[i] = KNM_VS1 * DNS_FL1;
			if(PTYPE[i] != 1)MEU[i] = KNM_VS2 * DNS_FL2;
		}}

//		if(TURB>0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2*sqrt(II[i]);

//			if(NEUt[i] * 0 != 0)  NEUt[i] = 0;
//			if(NEUt[i]>1)     NEUt[i] = 1;
//		}
	}

	// Granular Fluid
	if(Fluid2_type == 1)
	{
		// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
		// Loop only for particles near mesh
//#pragma omp parallel for
		//for(int im=0;im<nPartNearMesh;im++) {
		//int i = partNearMesh[im];
		for(int i=0; i<nP; i++) {
		//if(Typ[i] == FLD) {
		if(Typ[i] == FLD && Nw[i] == 1) {
//			double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];

			if(PTYPE[i] == 1) 
				MEU[i] = KNM_VS1 * DNS_FL1;
			else if(PTYPE[i] == 2) {
				phi = (Cv[i] - 0.25)*PHI / (1 - 0.25);
				phi2 = (Cv[i] - 0.25)*PHI_2 / (1 - 0.25);
				if(Cv[i] <= 0.25) { phi = 0.00001; phi2 = 0.00001; } // phi close to zero
				if(PTYPE[i] <= 0) phi = PHI_BED;

				// Normal stress calculation (mechanical pressure)
				p_rheo_new[i] = p_smooth[i];

				kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
				ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress=(BL[k]-pos_iy+DL/2)*(DNS_FL2)*9.81;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - pos_iz + PCL_DST / 2)*(DNS_FL2 - DNS_FL1)*9.81 - (vec_ix*vec_ix + vec_iy*vec_iy + vec_iz*vec_iz)*(DNS_FL2 - DNS_FL1) / 2.0;	// normal_stress= Gama.H

				if(p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81<0) p_smooth[i] = (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
				if(TIM <= 1) normal_stress = 1.0*(1 - TIM)*(p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81) + 1.0*(TIM)*normal_stress;

//				normal_stress = p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
//				normal_stress = p_smooth[i];
				//normal_stress=normal_stress*0.61*1500/DNS_FL2;
				if(normal_stress < 1 || Cv[i] < 0.5) normal_stress = 1;

				p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/DNS_SDT);		// Free-fall (dry granular material)
				Inertia[i] = sqrt(II[i])*DG/sqrt(normal_stress/(DNS_FL1*Cd));	// Grain inertia (submerged)
				//Inertia[i] = sqrt(II[i])*(KNM_VS1*DNS_FL1)/normal_stress ;	// Viscous regime

				// VF_max VF_min
				VF[i] = VF_max - (VF_max - VF_min)*Inertia[i];
				if(VF[i] < VF_min) VF[i] = VF_min;
				RHO[i] = DNS_SDT * VF[i] + (1-VF[i])*DNS_FL1;
				phi = phi * VF[i] / VF_max;

				double yield_stress = cohes * cos(phi) + normal_stress * sin(phi);

				if(yield_stress < 0) yield_stress = 0;

				double visc_max = (yield_stress*mm*0.5 + MEU0);

				if(II[i]>0)
					MEU_Y[i] = yield_stress * (1 - exp(-mm * sqrt(II[i]))) / 2.0 / sqrt(II[i]);
				else
					MEU_Y[i] = visc_max;

				// H-B rheology

				//meu_0 = MEU0;

				// Non-linear Meu(I) rheology
				//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/DNS_FL2)+sqrt(II[i])*DG);			//free fall
				meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*DG/(I0*sqrt(normal_stress/(DNS_FL1*Cd))+sqrt(II[i])*DG);		//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(KNM_VS1*DNS_FL1)/(I0*normal_stress+sqrt(II[i])*(KNM_VS1*DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL2)/I0;		//free fall
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*DG*sqrt(normal_stress*DNS_FL1*Cd)/I0;	//grain inertia
			   	//meu_0 = 0.5*(tan(phi2) - tan(phi))*(KNM_VS1*DNS_FL1)/I0;					//viscous

				if(II[i] <= 0 || (meu_0 * 0) != 0) meu_0 = MEU0;

				visc_max = (yield_stress*mm*0.5 + meu_0);

				// Herschel bulkley papanastasiou
				MEU[i] = MEU_Y[i] + MEU0 * pow(4 * II[i], (N - 1) / 2);

				// MEU_Y rheological model
				//MEU[i] = MEU_Y[i] + meu_0;
				
				if(II[i] == 0 || MEU[i]>visc_max) MEU[i] = visc_max;
				if(PTYPE[i] <= 0) MEU[i] = MEU[i] * Cv[i] + DNS_FL1*KNM_VS1*(1 - Cv[i]);
			}
			
			if(PTYPE[i] >= 2) {
				if(Cv[i] > 0.5) RHO[i] = DNS_FL2;
				else RHO[i] = Cv[i] * DNS_FL2 + (1 - Cv[i]) * DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if(stress_cal_method == 2)
//		{
//			for(i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
//				for(l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if(i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2 * MEU[i] * MEU[j] / (MEU[i] + MEU[j]);
//						if((NEUt[i] + NEUt[j])>0) meuij = meuij + 2 * NEUt[i] * RHO[i] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j]);

//						sum1 = sum1 + meuij * (x_vel[j] - x_vel[i])*DX(i, j)*w / d / d;
//						sum2 = sum2 + meuij * (x_vel[j] - x_vel[i])*DY(i, j)*w / d / d;
//						sum3 = sum3 + meuij * (x_vel[j] - x_vel[i])*DZ(i, j)*w / d / d;

//						sum4 = sum4 + meuij * (y_vel[j] - y_vel[i])*DX(i, j)*w / d / d;
//						sum5 = sum5 + meuij * (y_vel[j] - y_vel[i])*DY(i, j)*w / d / d;
//						sum6 = sum6 + meuij * (y_vel[j] - y_vel[i])*DZ(i, j)*w / d / d;

//						sum7 = sum7 + meuij * (z_vel[j] - z_vel[i])*DX(i, j)*w / d / d;
//						sum8 = sum8 + meuij * (z_vel[j] - z_vel[i])*DY(i, j)*w / d / d;
//						sum9 = sum9 + meuij * (z_vel[j] - z_vel[i])*DZ(i, j)*w / d / d;
//					}
//				}

//				Tau_xx[i] = (DIM / n0) * 2 * sum1;
//				Tau_yy[i] = (DIM / n0) * 2 * sum5;
//				Tau_zz[i] = (DIM / n0) * 2 * sum9;

//				Tau_xy[i] = (DIM / n0)*(sum2 + sum4);
//				Tau_xz[i] = (DIM / n0)*(sum3 + sum7);
//				Tau_yz[i] = (DIM / n0)*(sum6 + sum8);
//			}
//		}

	} // if(Fluid2_type == 1)

	//---------------------------------------------------------------

//	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS; //delete[]p_smooth;
//	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL;// p_smooth = NULL;
	
	delete[]BL; delete[]WL; delete[]PS;
	BL = NULL; WL = NULL; PS = NULL;
}

// Free-Slip condition. Viscosity interaction values (Polygon wall)
void WallSlipVscTrm_omp() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
	//if(Typ[i] == FLD) {
	if(Typ[i] == FLD && Nw[i] == 1) {
		
		double meu_i = MEU[i];
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
//		double vec_ix = Velk[i*3  ];	double vec_iy = Velk[i*3+1];	double vec_iz = Velk[i*3+2

		// Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
		double Rref_i[9], normaliw[3], normalMod2;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

	    if(normalMod2 > 0.00000001) {
	    	double normalMod = sqrt(normalMod2);
	    	normaliw[0] = normaliw[0]/normalMod;
	    	normaliw[1] = normaliw[1]/normalMod;
	    	normaliw[2] = normaliw[2]/normalMod;
	    }
	    else {
	    	normaliw[0] = 0;
	    	normaliw[1] = 0;
	    	normaliw[2] = 0;
	    }

	    //  Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
	    Rref_i[0] = 1.0 - 2*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2*normaliw[0]*normaliw[2];
		Rref_i[3] = 0.0 - 2*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2*normaliw[1]*normaliw[2];
		Rref_i[6] = 0.0 - 2*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2*normaliw[2]*normaliw[2];

		// Mirror particle velocity vi' = Rref_i * vi
      	double vec_mix = (Rref_i[0]*vec_ix + Rref_i[1]*vec_iy + Rref_i[2]*vec_iz);
		double vec_miy = (Rref_i[3]*vec_ix + Rref_i[4]*vec_iy + Rref_i[5]*vec_iz);
		double vec_miz = (Rref_i[6]*vec_ix + Rref_i[7]*vec_iy + Rref_i[8]*vec_iz);

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && dstimj2 < reL2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstimj2);
					double wL = WEI(dst, reL);

					NEU = 2 * meu_i * MEU[j] / (meu_i + MEU[j]);

//NEU = KNM_VS2 * DNS_FL2;

					if(PTYPE[i] == 1) NEU = NEU/DNS_FL1;
					else NEU = NEU/DNS_FL2;

					//if((NEUt[i] + NEUt[j])>0) NEU = NEU + (2 * NEUt[i] * RHO[j] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j])) / RHO[i];

					// Original
//					Acc_x +=(Vel[j*3  ]-vec_mix)*w;
//					Acc_y +=(Vel[j*3+1]-vec_miy)*w;
//					Acc_z +=(Vel[j*3+2]-vec_miz)*w;
					// Modified
					Acc_x +=(Vel[j*3  ]-vec_mix)*wL*NEU;
					Acc_y +=(Vel[j*3+1]-vec_miy)*wL*NEU;
					Acc_z +=(Vel[j*3+2]-vec_miz)*wL*NEU;

					//Acc_x +=(Velk[j*3  ]-vec_mix)*w;
					//Acc_y +=(Velk[j*3+1]-vec_miy)*w;
					//Acc_z +=(Velk[j*3+2]-vec_miz)*w;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

		// Add "i" contribution ("i" is a neighboor of "mirror i")
		double v0imi = pos_ix - pos_mix;
		double v1imi = pos_iy - pos_miy;
		double v2imi = pos_iz - pos_miz;
		double dstimi2 = v0imi*v0imi+v1imi*v1imi+v2imi*v2imi;
		if(dstimi2 < reL2) {
			double dst = sqrt(dstimi2);
			double wL = WEI(dst, reL);

			NEU = 2 * meu_i * meu_i / (meu_i + meu_i);

//NEU = KNM_VS2 * DNS_FL2;

			if(PTYPE[i] == 1) NEU = NEU/DNS_FL1;
			else NEU = NEU/DNS_FL2;

			// Original
//			Acc_x +=(vec_ix-vec_mix)*w;
//			Acc_y +=(vec_iy-vec_miy)*w;
//			Acc_z +=(vec_iz-vec_miz)*w;

			// Modified
			Acc_x +=(vec_ix-vec_mix)*wL*NEU;
			Acc_y +=(vec_iy-vec_miy)*wL*NEU;
			Acc_z +=(vec_iz-vec_miz)*wL*NEU;
	  	}

	  	// Wall laplacian Mitsume`s model
      	// Correction of velocity
      	// Original
//      Acc[i*3  ] += (Rref_i[0]*Acc_x + Rref_i[1]*Acc_y + Rref_i[2]*Acc_z)*A1;
//		Acc[i*3+1] += (Rref_i[3]*Acc_x + Rref_i[4]*Acc_y + Rref_i[5]*Acc_z)*A1;
//		Acc[i*3+2] += (Rref_i[6]*Acc_x + Rref_i[7]*Acc_y + Rref_i[8]*Acc_z)*A1;

	  	// A1_M = 2.0*DIM/(n0L*lmd);
		// Modified
		Acc[i*3  ] += (Rref_i[0]*Acc_x + Rref_i[1]*Acc_y + Rref_i[2]*Acc_z)*A1_M;
		Acc[i*3+1] += (Rref_i[3]*Acc_x + Rref_i[4]*Acc_y + Rref_i[5]*Acc_z)*A1_M;
		Acc[i*3+2] += (Rref_i[6]*Acc_x + Rref_i[7]*Acc_y + Rref_i[8]*Acc_z)*A1_M;
	}}
}

// No-Slip condition. Viscosity interaction values (Polygon wall)
void WallNoSlipVscTrm_omp() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
	//if(Typ[i] == FLD) {
	if(Typ[i] == FLD && Nw[i] == 1) {

		double meu_i = MEU[i];
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
//		double vec_ix = Velk[i*3  ];	double vec_iy = Velk[i*3+1];	double vec_iz = Velk[i*3+2];

		// Inverse matrix Rinv_i = - I
		double Rinv_i[9], normaliw[3], normalMod2;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

	    if(normalMod2 > 0.00000001) {
	    	double normalMod = sqrt(normalMod2);
	    	normaliw[0] = normaliw[0]/normalMod;
	    	normaliw[1] = normaliw[1]/normalMod;
	    	normaliw[2] = normaliw[2]/normalMod;
	    }
	    else {
	    	normaliw[0] = 0;
	    	normaliw[1] = 0;
	    	normaliw[2] = 0;
	    }

	    //  Inverse transformation matrix Rinv_i = - I
	    Rinv_i[0] = -1.0; Rinv_i[1] =  0.0; Rinv_i[2] =  0.0;
		Rinv_i[3] =  0.0; Rinv_i[4] = -1.0; Rinv_i[5] =  0.0;
		Rinv_i[6] =  0.0; Rinv_i[7] =  0.0; Rinv_i[8] = -1.0;

		double viwall[3], vtil[3];
		// Wall velocity (0 if fixed)
		viwall[0]=viwall[1]=viwall[2]=0.0;

#ifdef FORCED_ON
		if(meshID[i] == 2) {
			viwall[0] = velVWall[0];
			viwall[1] = velVWall[1];
			viwall[2] = velVWall[2];
		}
#endif

		// normal_iwall*v_iwall
		double dotnv = normaliw[0]*viwall[0] + normaliw[1]*viwall[1] + normaliw[2]*viwall[2];
		// vtil = vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}
		vtil[0] = vec_ix - 2*(viwall[0] - dotnv*normaliw[0]);
		vtil[1] = vec_iy - 2*(viwall[1] - dotnv*normaliw[1]);
		vtil[2] = vec_iz - 2*(viwall[2] - dotnv*normaliw[2]);
		// Mirror particle velocity vi' = Rinv_i * [vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}] 
      	double vec_mix = (Rinv_i[0]*vtil[0] + Rinv_i[1]*vtil[1] + Rinv_i[2]*vtil[2]);
		double vec_miy = (Rinv_i[3]*vtil[0] + Rinv_i[4]*vtil[1] + Rinv_i[5]*vtil[2]);
		double vec_miz = (Rinv_i[6]*vtil[0] + Rinv_i[7]*vtil[1] + Rinv_i[8]*vtil[2]);

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reL2 && dstimj2 < reL2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstimj2);
					double wL = WEI(dst, reL);

					NEU = 2 * meu_i * MEU[j] / (meu_i + MEU[j]);

//NEU = KNM_VS2 * DNS_FL2;


					if(PTYPE[i] == 1) NEU = NEU/DNS_FL1;
					else NEU = NEU/DNS_FL2;
					
					//if((NEUt[i] + NEUt[j])>0) NEU = NEU + (2 * NEUt[i] * RHO[j] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j])) / RHO[i];

					// Original
//					Acc_x +=(Vel[j*3  ]-vec_mix)*w;
//					Acc_y +=(Vel[j*3+1]-vec_miy)*w;
//					Acc_z +=(Vel[j*3+2]-vec_miz)*w;
					// Modified
					Acc_x +=(Vel[j*3  ]-vec_mix)*wL*NEU;
					Acc_y +=(Vel[j*3+1]-vec_miy)*wL*NEU;
					Acc_z +=(Vel[j*3+2]-vec_miz)*wL*NEU;
					
					//Acc_x +=(Velk[j*3  ]-vec_mix)*w;
					//Acc_y +=(Velk[j*3+1]-vec_miy)*w;
					//Acc_z +=(Velk[j*3+2]-vec_miz)*w;

					//if(i==2817) {
            		//	Fwall[i*3  ] += 1;//AA[0];
					//	Fwall[i*3+1] += 1;//AA[1];
					//	Fwall[i*3+2] += 1;//AA[2];
					//	std::cout << j << " " << Velk[j*3] << " " << Velk[j*3+1] << " " << Velk[j*3+2] << " Pj " << Prs[j] << std::endl;
          			//}
					//Acc_x += (Rinv_i[0]*(Velk[j*3  ]-vec_mix)+ Rinv_i[1]*(Velk[j*3+1]-vec_miy) + Rinv_i[2]*(Velk[j*3+2]-vec_miz))*w;
					//Acc_y += (Rinv_i[3]*(Velk[j*3  ]-vec_mix)+ Rinv_i[4]*(Velk[j*3+1]-vec_miy) + Rinv_i[5]*(Velk[j*3+2]-vec_miz))*w;
					//Acc_z += (Rinv_i[6]*(Velk[j*3  ]-vec_mix)+ Rinv_i[7]*(Velk[j*3+1]-vec_miy) + Rinv_i[8]*(Velk[j*3+2]-vec_miz))*w;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

		// Add "i" contribution ("i" is a neighboor of "mirror i")
	  	double v0imi = pos_ix - pos_mix;
		double v1imi = pos_iy - pos_miy;
		double v2imi = pos_iz - pos_miz;
		double dstimi2 = v0imi*v0imi+v1imi*v1imi+v2imi*v2imi;
		if(dstimi2 < reL2) {
			double dst = sqrt(dstimi2);
			double wL = WEI(dst, reL);

			NEU = 2 * meu_i * meu_i / (meu_i + meu_i);

//NEU = KNM_VS2 * DNS_FL2;

			if(PTYPE[i] == 1) NEU = NEU/DNS_FL1;
			else NEU = NEU/DNS_FL2;

			// Original
//			Acc_x +=(vec_ix-vec_mix)*w;
//			Acc_y +=(vec_iy-vec_miy)*w;
//			Acc_z +=(vec_iz-vec_miz)*w;

			// Modified
			Acc_x +=(vec_ix-vec_mix)*wL*NEU;
			Acc_y +=(vec_iy-vec_miy)*wL*NEU;
			Acc_z +=(vec_iz-vec_miz)*wL*NEU;

			//Acc_x += (Rinv_i[0]*(vec_ix-vec_mix)+ Rinv_i[1]*(vec_iy-vec_miy) + Rinv_i[2]*(vec_iz-vec_miz))*w;
			//Acc_y += (Rinv_i[3]*(vec_ix-vec_mix)+ Rinv_i[4]*(vec_iy-vec_miy) + Rinv_i[5]*(vec_iz-vec_miz))*w;
			//Acc_z += (Rinv_i[6]*(vec_ix-vec_mix)+ Rinv_i[7]*(vec_iy-vec_miy) + Rinv_i[8]*(vec_iz-vec_miz))*w;
	  	}

	  	//if(i==6) {
	  	//	printf("Vx:%lf Vy:%lf Vz:%lf Vx:%lf Vy:%lf Vz:%lf\n",vec_ix,vec_mix,vec_iy,vec_miy,vec_iz,vec_miz);
		//	printf("Accx:%lf Bccy:%lf Bccz:%lf \n", Acc[i*3],Acc[i*3+1],Acc[i*3+2]);
			//printf("Bccx:%lf Bccy:%lf Bccz:%lf \n", Acc[i*3],Acc[i*3+1],Acc[i*3+2]);
		//}
		// Wall laplacian Mitsume`s model
      	// Correction of velocity
      	// A1_M = 2.0*DIM/(n0L*lmd);
      	// Original
//     	Acc[i*3  ] += (Rinv_i[0]*Acc_x + Rinv_i[1]*Acc_y + Rinv_i[2]*Acc_z)*A1;
//		Acc[i*3+1] += (Rinv_i[3]*Acc_x + Rinv_i[4]*Acc_y + Rinv_i[5]*Acc_z)*A1;
//		Acc[i*3+2] += (Rinv_i[6]*Acc_x + Rinv_i[7]*Acc_y + Rinv_i[8]*Acc_z)*A1;
		// Modified
		Acc[i*3  ] += (Rinv_i[0]*Acc_x + Rinv_i[1]*Acc_y + Rinv_i[2]*Acc_z)*A1_M;
		Acc[i*3+1] += (Rinv_i[3]*Acc_x + Rinv_i[4]*Acc_y + Rinv_i[5]*Acc_z)*A1_M;
		Acc[i*3+2] += (Rinv_i[6]*Acc_x + Rinv_i[7]*Acc_y + Rinv_i[8]*Acc_z)*A1_M;
		//Acv[i*3  ] = (Rinv_i[0]*Acc_x + Rinv_i[1]*Acc_y + Rinv_i[2]*Acc_z)*A1;
		//Acv[i*3+1] = (Rinv_i[3]*Acc_x + Rinv_i[4]*Acc_y + Rinv_i[5]*Acc_z)*A1;
		//Acv[i*3+2] = (Rinv_i[6]*Acc_x + Rinv_i[7]*Acc_y + Rinv_i[8]*Acc_z)*A1;

		//double AA[3];
		
		//AA[0] = (Rinv_i[0]*Acc_x + Rinv_i[1]*Acc_y + Rinv_i[2]*Acc_z)*A1;
		//AA[1] = (Rinv_i[3]*Acc_x + Rinv_i[4]*Acc_y + Rinv_i[5]*Acc_z)*A1;
		//AA[2] = (Rinv_i[6]*Acc_x + Rinv_i[7]*Acc_y + Rinv_i[8]*Acc_z)*A1;

		//F1[i*3  ] = AA[0];
		//F1[i*3+1] = AA[1];
		//F1[i*3+2] = AA[2];

		//if(i==2817) {
			//Fwall[i*3  ] = AA[0];
			//Fwall[i*3+1] = AA[1];
			//Fwall[i*3+2] = AA[2];
			//std::cout << "t: " << TIM << std::endl;
			//std::cout << "Fwall " << AA[0] << " " << AA[1] << " " << AA[2] << std::endl;
			//std::cout << "Veli " << vec_ix << " " << vec_iy << " " << vec_iz << std::endl;
			//std::cout << "Posmi " << pos_mix << " " << pos_miy << " " << pos_miz << std::endl;
			//std::cout << "pndi " << pndi[i] << " Pi " << Prs[i] << std::endl;
			//printf("Accx:%lf Accy:%lf Accz:%lf \n", AA[0],AA[1],AA[2]);
		//}
		//if(i==6) {
			//printf("Accx:%lf Accy:%lf Accz:%lf \n", Acc[i*3],Acc[i*3+1],Acc[i*3+2]);
		//	printf("Time:%e\n", TIM);
		//	printf("Xi:%e %e %e Xm:%e %e %e\n", pos_ix,pos_iy,pos_iz,pos_mix,pos_miy,pos_miz);
		//	printf("Vi:%e %e %e Vm:%e %e %e\n", vec_ix,vec_iy,vec_iz,vec_mix,vec_miy,vec_miz);
		//	printf("Acc:%e %e %e\n", AA[0],AA[1],AA[2]);
		//}
	}}
}

// Update velocity and positions
void UpPcl2_omp(void) {
	vMax = 0.0;						// Maximum velocity
	double auxiliar[5] = {1.2, -3.3, 4.3, -0.3, 5.6};

	// https://stackoverflow.com/questions/39989473/use-openmp-in-c11-to-find-the-maximum-of-the-calculated-values
#pragma omp parallel
{
	double local_vMax = 0.0;
#pragma omp for
	for(int i=0; i<nP; i++) {
//		if(Typ[i] == FLD) {
		if(Typ[i] != GST) {
			Vel[i*3  ]+=Acc[i*3  ]*DT;	Vel[i*3+1]+=Acc[i*3+1]*DT;	Vel[i*3+2]+=Acc[i*3+2]*DT;
			if(Typ[i] == FLD) {
				Pos[i*3  ]+=Acc[i*3  ]*DT*DT;	Pos[i*3+1]+=Acc[i*3+1]*DT*DT;	Pos[i*3+2]+=Acc[i*3+2]*DT*DT;
				ChkPcl(i);
			}
			Acc[i*3]=Acc[i*3+1]=Acc[i*3+2]=0.0;

			//Pos[i*3  ]=Posk[i*3 ]+Vel[i*3  ]*DT;	Pos[i*3+1]=Posk[i*3+1]+Vel[i*3+1]*DT;	Pos[i*3+2]=Posk[i*3+2]+Vel[i*3+2]*DT;
			//Posk[i*3  ]=Pos[i*3  ];	Posk[i*3+1]=Pos[i*3+1];	Posk[i*3+2]=Pos[i*3+2];
			//Velk[i*3  ]=Vel[i*3  ];	Velk[i*3+1]=Vel[i*3+1];	Velk[i*3+2]=Vel[i*3+2];

			//F1[i*3  ]=Acc[i*3  ];	F1[i*3+1]=Acc[i*3+1];	F1[i*3+2]=Acc[i*3+2];
			//F2[i*3  ]=Acv[i*3  ];	F2[i*3+1]=Acv[i*3+1];	F2[i*3+2]=Acv[i*3+2];
			//Acc[i*3]=Acc[i*3+1]=Acc[i*3+2]=0.0;
			//Acv[i*3]=Acv[i*3+1]=Acv[i*3+2]=0.0;

			double vMod2 = Vel[i*3  ]*Vel[i*3  ] + Vel[i*3+1]*Vel[i*3+1] + Vel[i*3+2]*Vel[i*3+2];
			if(vMod2 > local_vMax*local_vMax)
				local_vMax = sqrt(vMod2);
		}
	}

#pragma omp critical
	{
		if (local_vMax > vMax)
			vMax = local_vMax;
	}
}
	Crt = DT*vMax/PCL_DST;
}

// Shifting technique
// Improvements for accuracy and stability in a weakly-compressible particle method
// https://www.sciencedirect.com/science/article/pii/S0045793016302250
void AdjVel1_omp() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double du_ix = 0.0;	double du_iy = 0.0;	double du_iz = 0.0;
		//double ni = 0.0;
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double dw = DEL_WEI(dst, reS);
					du_ix += dw*(Vel[j*3  ]-vec_ix);
					du_iy += dw*(Vel[j*3+1]-vec_iy);
					du_iz += dw*(Vel[j*3+2]-vec_iz);
					//double w = WEI(dst, r);
					//ni += w;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

//#ifdef POLYGON_ON
//		if(pndi[i] > PND_TRS*n0S)
//#else
//		if(pndi[i] > PND_TRS*n0S || numNeigh[i] > NGH_TRS*nNeigh0)
//#endif
//		if(pndS[i] > PND_TRS*n0S || numNeigh[i] > NGH_TRS*nNeigh0)
		if(Bc[i] == INR) {
			Vel[i*3  ] -= A5*du_ix;
			Vel[i*3+1] -= A5*du_iy;
			Vel[i*3+2] -= A5*du_iz;
		}
		// else {
		// 	double Inn[9], duAux[3];
		// 	// I - nxn
		// 	Inn[0] = 1.0 - Normal[i*3  ]*Normal[i*3  ]; Inn[1] = 0.0 - Normal[i*3  ]*Normal[i*3+1]; Inn[2] = 0.0 - Normal[i*3  ]*Normal[i*3+2];
		// 	Inn[3] = 0.0 - Normal[i*3+1]*Normal[i*3  ]; Inn[4] = 1.0 - Normal[i*3+1]*Normal[i*3+1]; Inn[5] = 0.0 - Normal[i*3+1]*Normal[i*3+2];
		// 	Inn[6] = 0.0 - Normal[i*3+2]*Normal[i*3  ]; Inn[7] = 0.0 - Normal[i*3+2]*Normal[i*3+1]; Inn[8] = 1.0 - Normal[i*3+2]*Normal[i*3+2];
		// 	// (I - nxn)dr
		// 	duAux[0] = Inn[0]*du_ix + Inn[1]*du_iy + Inn[2]*du_iz;
		// 	duAux[1] = Inn[3]*du_ix + Inn[4]*du_iy + Inn[5]*du_iz;
		// 	duAux[2] = Inn[6]*du_ix + Inn[7]*du_iy + Inn[8]*du_iz;
		// 	Vel[i*3  ] -= A5*duAux[0];
		// 	Vel[i*3+1] -= A5*duAux[1];
		// 	Vel[i*3+2] -= A5*duAux[2];
		// }
	}}
}

// Normal vector on the fluid
// An accurate and stable multiphase moving particle semi-implicit method based on a corrective matrix for all particle interaction models
// https://onlinelibrary.wiley.com/doi/full/10.1002/nme.5844
void Normal_omp() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double dr_ix = 0.0;	double dr_iy = 0.0;	double dr_iz = 0.0;
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wS = WEI(dst, reS);
					dr_ix += v0ij*wS/dst;
					dr_iy += v1ij*wS/dst;
					dr_iz += v2ij*wS/dst;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		double norm2 = dr_ix*dr_ix + dr_iy*dr_iy + dr_iz*dr_iz;
		if(norm2 > 0.0) {
			double norm = sqrt(norm2);
			Normal[i*3  ] = dr_ix/norm;
			Normal[i*3+1] = dr_iy/norm;
			Normal[i*3+2] = dr_iz/norm;
		}
		else {
			Normal[i*3  ] = 0.0;
			Normal[i*3+1] = 0.0;
			Normal[i*3+2] = 0.0;
		}
	}}
}

// Shifting technique (Polygon wall)
void WallAdjVel1_omp() {
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
	//if(Typ[i] == FLD) {
	if(Typ[i] == FLD && Nw[i] == 1) {
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double du_ix = 0.0;	double du_iy = 0.0;	double du_iz = 0.0;

		// No-slip

		// Inverse matrix Rinv_i = - I
		double Rinv_i[9], normaliw[3], normalMod2;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normalMod2 = normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2];

	    if(normalMod2 > 0.00000001) {
	    	double normalMod = sqrt(normalMod2);
	    	normaliw[0] = normaliw[0]/normalMod;
	    	normaliw[1] = normaliw[1]/normalMod;
	    	normaliw[2] = normaliw[2]/normalMod;
	    }
	    else {
	    	normaliw[0] = 0;
	    	normaliw[1] = 0;
	    	normaliw[2] = 0;
	    }

	    // Inverse transformation matrix Rinv_i = - I
	    Rinv_i[0] = -1.0; Rinv_i[1] =  0.0; Rinv_i[2] =  0.0;
		Rinv_i[3] =  0.0; Rinv_i[4] = -1.0; Rinv_i[5] =  0.0;
		Rinv_i[6] =  0.0; Rinv_i[7] =  0.0; Rinv_i[8] = -1.0;

		double viwall[3], vtil[3];
		// Wall velocity (0 if fixed)
		viwall[0]=viwall[1]=viwall[2]=0.0;

#ifdef FORCED_ON
		if(meshID[i] == 2) {
			viwall[0] = velVWall[0];
			viwall[1] = velVWall[1];
			viwall[2] = velVWall[2];
		}
#endif

		// normal_iwall*v_iwall
		double dotnv = normaliw[0]*viwall[0] + normaliw[1]*viwall[1] + normaliw[2]*viwall[2];
		// vtil = vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}
		vtil[0] = vec_ix - 2*(viwall[0] - dotnv*normaliw[0]);
		vtil[1] = vec_iy - 2*(viwall[1] - dotnv*normaliw[1]);
		vtil[2] = vec_iz - 2*(viwall[2] - dotnv*normaliw[2]);
		// Mirror particle velocity vi' = Rinv_i * [vi - 2 {v_iwall - (normal_iwall*v_iwall)normal_iwall}] 
      	double vec_mix = (Rinv_i[0]*vtil[0] + Rinv_i[1]*vtil[1] + Rinv_i[2]*vtil[2]);
		double vec_miy = (Rinv_i[3]*vtil[0] + Rinv_i[4]*vtil[1] + Rinv_i[5]*vtil[2]);
		double vec_miz = (Rinv_i[6]*vtil[0] + Rinv_i[7]*vtil[1] + Rinv_i[8]*vtil[2]);

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstimj2);
					double wS = DEL_WEI(dst, reS);
					du_ix += wS*(Vel[j*3  ]-vec_mix);
					du_iy += wS*(Vel[j*3+1]-vec_miy);
					du_iz += wS*(Vel[j*3+2]-vec_miz);
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		// Add "i" contribution ("i" is a neighboor of "mirror i")
	  	double v0imi = pos_ix - pos_mix;
		double v1imi = pos_iy - pos_miy;
		double v2imi = pos_iz - pos_miz;
		double dstimi2 = v0imi*v0imi+v1imi*v1imi+v2imi*v2imi;
		if(dstimi2 < reS2) {
			double dst = sqrt(dstimi2);
			double wS = DEL_WEI(dst, reS);
			du_ix += wS*(vec_ix-vec_mix);
			du_iy += wS*(vec_iy-vec_miy);
			du_iz += wS*(vec_iz-vec_miz);
	  	}

//#ifdef POLYGON_ON
//		if(pndi[i] > PND_TRS*n0S)
//#else
//		if(pndi[i] > PND_TRS*n0S || numNeigh[i] > NGH_TRS*nNeigh0)
//#endif
//		if(pndS[i] > PND_TRS*n0S || numNeigh[i] > NGH_TRS*nNeigh0)
		if(Bc[i] == INR) {
			double dux = Rinv_i[0]*du_ix + Rinv_i[1]*du_iy + Rinv_i[2]*du_iz;
	  		double duy = Rinv_i[3]*du_ix + Rinv_i[4]*du_iy + Rinv_i[5]*du_iz;
	  		double duz = Rinv_i[6]*du_ix + Rinv_i[7]*du_iy + Rinv_i[8]*du_iz;
			Vel[i*3  ] -= A5*dux;
			Vel[i*3+1] -= A5*duy;
			Vel[i*3+2] -= A5*duz;
		}
		// else {
		// 	double Inn[9], duAux[3];
		// 	// I - nxn
		// 	Inn[0] = 1.0 - Normal[i*3  ]*Normal[i*3  ]; Inn[1] = 0.0 - Normal[i*3  ]*Normal[i*3+1]; Inn[2] = 0.0 - Normal[i*3  ]*Normal[i*3+2];
		// 	Inn[3] = 0.0 - Normal[i*3+1]*Normal[i*3  ]; Inn[4] = 1.0 - Normal[i*3+1]*Normal[i*3+1]; Inn[5] = 0.0 - Normal[i*3+1]*Normal[i*3+2];
		// 	Inn[6] = 0.0 - Normal[i*3+2]*Normal[i*3  ]; Inn[7] = 0.0 - Normal[i*3+2]*Normal[i*3+1]; Inn[8] = 1.0 - Normal[i*3+2]*Normal[i*3+2];
		// 	double dux = Rinv_i[0]*du_ix + Rinv_i[1]*du_iy + Rinv_i[2]*du_iz;
	 //  		double duy = Rinv_i[3]*du_ix + Rinv_i[4]*du_iy + Rinv_i[5]*du_iz;
	 //  		double duz = Rinv_i[6]*du_ix + Rinv_i[7]*du_iy + Rinv_i[8]*du_iz;
		// 	// (I - nxn)dr
		// 	duAux[0] = Inn[0]*dux + Inn[1]*duy + Inn[2]*duz;
		// 	duAux[1] = Inn[3]*dux + Inn[4]*duy + Inn[5]*duz;
		// 	duAux[2] = Inn[6]*dux + Inn[7]*duy + Inn[8]*duz;
		// 	Vel[i*3  ] -= A5*duAux[0];
		// 	Vel[i*3+1] -= A5*duAux[1];
		// 	Vel[i*3+2] -= A5*duAux[2];
		// }
	}}
}

// Concentration and Gradient of concentration
void GradConc_omp() {
	// Concentration
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
		Ci[i] = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wS = WEI(dst, reS);
					Ci[i] += wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}

		// Add PND due Wall polygon
		Ci[i] += niw[i];

		Ci[i] /= n0S;
	}}
	// Gradient of concentration
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] == FLD) {
		gradCi[i*3  ] = 0.0;	gradCi[i*3+1] = 0.0;	gradCi[i*3+2] = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double du_ix = 0.0;	double du_iy = 0.0;	double du_iz = 0.0;
		double C_i = Ci[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;
				// If j is inside the neighborhood of i and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstij2);
					double wS = WEI(dst, reS);
					gradCi[i*3  ] += (C_i + Ci[j])*v0ij*wS/dstij2;
					gradCi[i*3+1] += (C_i + Ci[j])*v1ij*wS/dstij2;
					gradCi[i*3+2] += (C_i + Ci[j])*v2ij*wS/dstij2;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		//A3 = -DIM/n0Grad
		gradCi[3*i  ] *= -A3;
		gradCi[3*i+1] *= -A3;
		gradCi[3*i+2] *= -A3;

//#ifdef POLYGON_ON
//		if(pndi[i] > PND_TRS*n0S)
//#else
//		if(pndi[i] > PND_TRS*n0S || numNeigh[i] > NGH_TRS*nNeigh0)
//#endif
//		if(pndS[i] > PND_TRS*n0S || numNeigh[i] > NGH_TRS*nNeigh0)
		if(Bc[i] == INR) {
			// A6 = COEF_A*PCL_DST*PCL_DST*CRT_NUM*MACH;	// Coefficient used to adjust velocity
			Pos[i*3  ] -= A6*gradCi[3*i  ];
			Pos[i*3+1] -= A6*gradCi[3*i+1];
			Pos[i*3+2] -= A6*gradCi[3*i+2];
		}
		/*
		else
		{
			double Inn[9], drAux[3];
		 	// I - nxn
		 	Inn[0] = 1.0 - Normal[i*3  ]*Normal[i*3  ]; Inn[1] = 0.0 - Normal[i*3  ]*Normal[i*3+1]; Inn[2] = 0.0 - Normal[i*3  ]*Normal[i*3+2];
		 	Inn[3] = 0.0 - Normal[i*3+1]*Normal[i*3  ]; Inn[4] = 1.0 - Normal[i*3+1]*Normal[i*3+1]; Inn[5] = 0.0 - Normal[i*3+1]*Normal[i*3+2];
		 	Inn[6] = 0.0 - Normal[i*3+2]*Normal[i*3  ]; Inn[7] = 0.0 - Normal[i*3+2]*Normal[i*3+1]; Inn[8] = 1.0 - Normal[i*3+2]*Normal[i*3+2];
		 	// (I - nxn)dr
		 	drAux[0] = Inn[0]*gradCi[3*i  ] + Inn[1]*gradCi[3*i+1] + Inn[2]*gradCi[3*i+2];
		 	drAux[1] = Inn[3]*gradCi[3*i  ] + Inn[4]*gradCi[3*i+1] + Inn[5]*gradCi[3*i+2];
		 	drAux[2] = Inn[6]*gradCi[3*i  ] + Inn[7]*gradCi[3*i+1] + Inn[8]*gradCi[3*i+2];
		 	// A6 = COEF_A*PCL_DST*PCL_DST*CRT_NUM*MACH;	// Coefficient used to adjust velocity
			Pos[i*3  ] -= A6*drAux[0];
			Pos[i*3+1] -= A6*drAux[1];
			Pos[i*3+2] -= A6*drAux[2];
		}
		*/
		// else {
		// 	double Inn[9], duAux[3];
		// 	// I - nxn
		// 	Inn[0] = 1.0 - Normal[i*3  ]*Normal[i*3  ]; Inn[1] = 0.0 - Normal[i*3  ]*Normal[i*3+1]; Inn[2] = 0.0 - Normal[i*3  ]*Normal[i*3+2];
		// 	Inn[3] = 0.0 - Normal[i*3+1]*Normal[i*3  ]; Inn[4] = 1.0 - Normal[i*3+1]*Normal[i*3+1]; Inn[5] = 0.0 - Normal[i*3+1]*Normal[i*3+2];
		// 	Inn[6] = 0.0 - Normal[i*3+2]*Normal[i*3  ]; Inn[7] = 0.0 - Normal[i*3+2]*Normal[i*3+1]; Inn[8] = 1.0 - Normal[i*3+2]*Normal[i*3+2];
		// 	// (I - nxn)dr
		// 	duAux[0] = Inn[0]*du_ix + Inn[1]*du_iy + Inn[2]*du_iz;
		// 	duAux[1] = Inn[3]*du_ix + Inn[4]*du_iy + Inn[5]*du_iz;
		// 	duAux[2] = Inn[6]*du_ix + Inn[7]*du_iy + Inn[8]*du_iz;
		// 	Vel[i*3  ] -= A5*duAux[0];
		// 	Vel[i*3+1] -= A5*duAux[1];
		// 	Vel[i*3+2] -= A5*duAux[2];
		// }
	}}
}

// Concentration and Gradient of concentration (Polygon wall)
void WallGradConc_omp() {
	// Gradient of concentration due Polygon wall
	//int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	//for(int im=0;im<nPartNearMesh;im++) {
	//int i = partNearMesh[im];
	for(int i=0; i<nP; i++) {
	//if(Typ[i] == FLD) {
	if(Typ[i] == FLD && Nw[i] == 1) {
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double dr_x = 0.0;			double dr_y = 0.0;			double dr_z = 0.0;
		double gradCiWall_x = 0.0;			double gradCiWall_y = 0.0;			double gradCiWall_z = 0.0;
		double ni = pndi[i];
		double C_i = Ci[i];

		// Wall gradient Mitsume`s model
	    double Rref_i[9], normaliw[3], normaliwSqrt;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);

	    if(normaliwSqrt > 0.00000001) {
	    	normaliw[0] = normaliw[0]/normaliwSqrt;
	    	normaliw[1] = normaliw[1]/normaliwSqrt;
	    	normaliw[2] = normaliw[2]/normaliwSqrt;
	    }
	    else {
	    	normaliw[0] = 0;
	    	normaliw[1] = 0;
	    	normaliw[2] = 0;
	    }

	    //  Transformation matrix Rref_i = I - 2*normal_iwall*normal_iwall
	    Rref_i[0] = 1.0 - 2*normaliw[0]*normaliw[0]; Rref_i[1] = 0.0 - 2*normaliw[0]*normaliw[1]; Rref_i[2] = 0.0 - 2*normaliw[0]*normaliw[2];
		Rref_i[3] = 0.0 - 2*normaliw[1]*normaliw[0]; Rref_i[4] = 1.0 - 2*normaliw[1]*normaliw[1]; Rref_i[5] = 0.0 - 2*normaliw[1]*normaliw[2];
		Rref_i[6] = 0.0 - 2*normaliw[2]*normaliw[0]; Rref_i[7] = 0.0 - 2*normaliw[2]*normaliw[1]; Rref_i[8] = 1.0 - 2*normaliw[2]*normaliw[2];

		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0imj = Pos[j*3  ] - pos_mix;
				double v1imj = Pos[j*3+1] - pos_miy;
				double v2imj = Pos[j*3+2] - pos_miz;

				double dstimj2 = v0imj*v0imj+v1imj*v1imj+v2imj*v2imj;

				// If j is inside the neighborhood of i and im (intersection) and 
				// is not at the same side of im (avoid real j in the virtual neihborhood)
				if(dstij2 < reS2 && dstimj2 < reS2 && dstij2 < dstimj2) {
				if(j != i && Typ[j] != GST) {
					double dst = sqrt(dstimj2);
					double wS = WEI_GRAD(dst, reS);

					dr_x += (C_i + Ci[j])*v0imj*wS/dstimj2;
					dr_y += (C_i + Ci[j])*v1imj*wS/dstimj2;
					dr_z += (C_i + Ci[j])*v2imj*wS/dstimj2;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		// Add "i" contribution ("i" is a neighboor of "mirror i")
	  	double v0imi = pos_ix - pos_mix;
		double v1imi = pos_iy - pos_miy;
		double v2imi = pos_iz - pos_miz;
		double dstimi2 = v0imi*v0imi+v1imi*v1imi+v2imi*v2imi;
		if(dstimi2 < reS2) {
			double dst = sqrt(dstimi2);
			double wS = WEI_GRAD(dst, reS);

			dr_x += (C_i + C_i)*v0imi*wS/dstimi2;
			dr_y += (C_i + C_i)*v1imi*wS/dstimi2;
			dr_z += (C_i + C_i)*v2imi*wS/dstimi2;
	  	}

		// A3 is a negative cte (-DIM/noGrad)
		// Original
//		Acc[i*3  ] += (RLX_PRS*(Rref_i[0]*Acc_x + Rref_i[1]*Acc_y + Rref_i[2]*Acc_z)*A3 - rpsForce[0])*invDns[FLD];
//		Acc[i*3+1] += (RLX_PRS*(Rref_i[3]*Acc_x + Rref_i[4]*Acc_y + Rref_i[5]*Acc_z)*A3 - rpsForce[1])*invDns[FLD];
//		Acc[i*3+2] += (RLX_PRS*(Rref_i[6]*Acc_x + Rref_i[7]*Acc_y + Rref_i[8]*Acc_z)*A3 - rpsForce[2])*invDns[FLD];
		// Modified
		gradCiWall_x = -(Rref_i[0]*dr_x + Rref_i[1]*dr_y + Rref_i[2]*dr_z)*A3;
		gradCiWall_y = -(Rref_i[3]*dr_x + Rref_i[4]*dr_y + Rref_i[5]*dr_z)*A3;
		gradCiWall_z = -(Rref_i[6]*dr_x + Rref_i[7]*dr_y + Rref_i[8]*dr_z)*A3;

		gradCi[i*3  ] += gradCiWall_x;
		gradCi[i*3+1] += gradCiWall_y;
		gradCi[i*3+2] += gradCiWall_z;

		if(Bc[i] == INR) {
			// A6 = COEF_A*PCL_DST*PCL_DST*CRT_NUM*MACH;	// Coefficient used to adjust velocity
			Pos[i*3  ] -= A6*gradCiWall_x;
			Pos[i*3+1] -= A6*gradCiWall_y;
			Pos[i*3+2] -= A6*gradCiWall_z;
		}
	}}
}

// Normal vector on the fluid
void NormalConc_omp() {
#pragma omp parallel for
	for(int i=0; i<nP; i++) {
	if(Typ[i] != GST) {
		double norm2GradCi = gradCi[3*i]*gradCi[3*i] + gradCi[3*i+1]*gradCi[3*i+1] + gradCi[3*i+2]*gradCi[3*i+2];
		
		if(norm2GradCi > 0.0) {
			double norm = sqrt(norm2GradCi);
			Normal[i*3  ] = -gradCi[i*3  ]/norm;
			Normal[i*3+1] = -gradCi[i*3+1]/norm;
			Normal[i*3+2] = -gradCi[i*3+2]/norm;
		}
		else {
			Normal[i*3  ] = 0.0;
			Normal[i*3+1] = 0.0;
			Normal[i*3+2] = 0.0;
		}
	}}
}

// Update velocity at wall and dummy particles
void UpVelWalLDummy_omp() {
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0; i<nP; i++) {
	if(Typ[i] == WLL) {
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		double du_ix = 0.0;	double du_iy = 0.0;	double du_iz = 0.0;
		double ni = 0.0;
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++) {
		for(int jy=iy-1;jy<=iy+1;jy++) {
		for(int jx=ix-1;jx<=ix+1;jx++) {
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;) {
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2 < reS2) {
//				if(j != i && Typ[j] != GST) {
				if(j != i && Typ[j] == FLD) {
					double dst = sqrt(dst2);
					double wS = WEI(dst, reS);
					ni += wS;
					du_ix += Vel[j*3  ]*wS;
					du_iy += Vel[j*3+1]*wS;
					du_iz += Vel[j*3+2]*wS;
				}}
				j = nxt[j];
				if(j == -1) break;
			}
		}}}
		if(ni > 0.0) {
			Vel[i*3  ] = 2*Vel[i*3  ] - du_ix/ni;
			Vel[i*3+1] = 2*Vel[i*3+1] - du_iy/ni;
			Vel[i*3+2] = 2*Vel[i*3+2] - du_iz/ni;
		}
		else {
			Vel[i*3  ] = 2*Vel[i*3  ] - du_ix;
			Vel[i*3+1] = 2*Vel[i*3+1] - du_iy;
			Vel[i*3+2] = 2*Vel[i*3+2] - du_iz;
		}
	}}
}

// Main loop
void ClcEMPS(mesh* &mesh) {
	while(1) {
		if(iLP%10 == 0) {
			int p_num = 0;
#pragma omp parallel for reduction(+: p_num)
			for(int i=0; i<nP; i++) {
				if(Typ[i] != GST)
					p_num++;
			}
			timer_end = get_dtime();
			printf("%5d th TIM: %lf / p_num: %d / Vmax: %lf / Courant: %lf", iLP,TIM,p_num, vMax, Crt);

			//printf("/ NearMesh %d ", partNearMesh.size());

			int seconds, hours, minutes;
			seconds = int(timer_end -timer_sta);
			minutes = seconds / 60;
			hours = minutes / 60;
			printf(" / RunTime: %02dh%02dm%02dsec\n", int(hours), int(minutes%60), int(seconds%60));
		}
		if(iLP%OPT_FQC == 0 ) {
			//WrtDat();
#ifdef BINARY
			WrtVtuBinary(); // - NOT WORKING !!!
#else
			WrtVtuAscii();
#endif

#ifdef FEM_ON
			// Write deformable mesh
			mesh[FEM_ID].writeMeshFile(FEM_ID, OUT_FOLDER, iF-1);
#endif
#ifdef FORCED_ON
			// Write forced rigid mesh
			mesh[FRW_ID].writeMeshFile(FRW_ID, OUT_FOLDER, iF-1);
#endif

			if(TIM >= FIN_TIM ) { break;}
		}

		// Small subdomains (buckets)
		MkBkt();
		
//		if(Fluid2_type==1) {
//			VolFract_omp(); // Calculation of the volume of fraction if phase II in the mixture
//			VscIntVal_omp();// Calculation of viscosity
//		}
		// Calculation of acceleration due laplacian of viscosity and gravity
		VscTrm_omp();
		// Add acceleration due pressure gradient (Prediction)
		if(RLX_PRS < 1.0) {
			PrdPrsGrdTrm_omp();
#ifdef POLYGON_ON
			// Pressure gradient on polygon wall
			WallPrdPrsGrdTrm_omp();
#endif
		}
		// Update velocity and positions. Set some variables to zero or inf
		UpPcl1_omp();
		// Check collision
		ChkCol_omp();

#ifdef POLYGON_ON
		// Fluid particles: Calculation of PND due wall and number of neighboors
		// Positions of wall and mirror particles
		for(int me = 0; me < NUM_MESHS; me++) {
			mesh[me].closestPointPNDBoundaryAABB(reS2, reL2, nP, wijType, Typ, FLD, GST, me, STA_ID, FEM_ID, FRW_ID, Pos, wallPos, mirrorPos, riw2, elementID, meshID, NormalWall);
//			mesh[me].closestPointPNDBoundaryAABB(reS2, reL2, nP, wijType, Typ, FLD, GST, Pos, wallPos, mirrorPos, riw2, niw, numNeighw, elementID, partNearMesh);
		}
		partNearMesh.clear();
		//for(int me = 0; me < NUM_MESHS; me++) {
		//	mesh[me].updateparticlesNearMesh(reS2, reL2, nP, wijType, Typ, FLD, GST, riw2, niw, numNeighw, partNearMesh);
		//}
		// Only call once since the distance riw2 have the values from all the meshes
		mesh[0].updateparticlesNearMesh(reS2, reL2, nP, wijType, Typ, FLD, GST, riw2, niw, numNeighw, partNearMesh, Nw);
#endif
#ifdef POLYGON_ON
		// NPCD PND due polygon wall
		WallNPCD_omp();
#endif
		// PND calculation
		MkPND_omp();
		// Diffusion term
		if(PND_CAL == 2) {
			DifTrm_omp();
#ifdef POLYGON_ON
			// Add contribution of PND due polygon wall
			if(SLP)
				WallSlipDifTrm_omp();
			else
				WallNoSlipDifTrm_omp();
#endif
			// Update PND
			UpPND_omp();
			// PND wall and dummy
			MeanPNDWallDummy_omp();
		}
		// Mean PND
		if(PND_CAL == 1) {
#ifdef POLYGON_ON
			// Contribution due polygon wall
			WallMeanPND_omp();
#endif
			MeanPND_omp();
		}
		// Pressure calculation and type of particle
		if(MPS_TYP == 0)
			MkPrs_omp();
		else if(MPS_TYP == 1)
			MkPrsWc_omp();
		// Wall and dummy pressure
		MkPrsWallDummy_omp();
#ifdef POLYGON_ON
		// Pressure at inner particles near corners
//		MkPrsNearWall_omp();
#endif
		// Correction matrix
		if(GRD_COR == 1)
			CorMtxTrm_omp();
		// Calculation of acceleration due pressure gradient
		PrsGrdTrm_omp();
#ifdef POLYGON_ON
		// Add acceleration due pressure gradient on polygon wall
		WallPrsGrdTrm_omp();
#endif
		// Add acceleration due laplacian of viscosity on wall
		if(SLP) {
			if(Fluid2_type == 1) {
				VolFract_omp(); // Calculation of the volume of fraction if phase II in the mixture

				VscIntVal_omp();// Calculation of viscosity due fluid particles
#ifdef POLYGON_ON
				WallSlipVscIntVal_omp();// Calculation of viscosity due polygon wall
#endif
			}
#ifdef POLYGON_ON
			WallSlipVscTrm_omp(); // Slip
#endif
		}
		else {
			if(Fluid2_type == 1) {
				VolFract_omp(); // Calculation of the volume of fraction if phase II in the mixture

				VscIntVal_omp();// Calculation of viscosity due fluid particles
#ifdef POLYGON_ON
				WallNoSlipVscIntVal_omp();// Calculation of viscosity due polygon wall
#endif
			}
#ifdef POLYGON_ON
			WallNoSlipVscTrm_omp(); // No-slip
#endif
		}

#ifdef FORCED_ON
		// Update forced mesh
		mesh[FRW_ID].updateForcedMesh(nodeFRWX, nodeFRWY, nodeFRWZ, velVWall, DT, TIM);
#endif

		// Update velocity and positions
		UpPcl2_omp();
		// Adjust velocity
		if(ADJ_VEL == 1) {
			// MatrixCorr_omp();
			//Normal_omp();
			AdjVel1_omp();
			//WallAdjVel1_omp();
		}
		else if(ADJ_VEL == 2) {
//			NormalConc_omp();
			GradConc_omp();
#ifdef POLYGON_ON
			WallGradConc_omp();
#endif
		}
		// Wall and dummy velocity
		UpVelWalLDummy_omp();
		// Pressure calculation
//		if(MPS_TYP == 0)
//			MkPrs_omp();
//		else if(MPS_TYP == 1)
//			MkPrsWc_omp();

		for(int i=0; i<nP; i++) {pav[i] += Prs[i];}
		iLP++;
		TIM += DT;
	}
}

// Main
int main( int argc, char** argv) {
	printf("start emps_omp.\n");
    // Creating a directory 
    if(mkdir(OUT_FOLDER, 0777) == -1) {
    	printf("Unable to create directory.\n");
    }
    else {
    	printf("Directory created.\n");
    }
	// Read and allocte memory for data
	RdDat();

	// Weight function
	wijType = WGT_TYP;

	// Read mesh
	//RdMes();

	mesh* solidMesh = NULL;
	solidMesh = new mesh[NUM_MESHS];
	

	//mesh mesh[NUM_MESHS];
	//meshs = new mesh[2];

	for(int me = 0; me < NUM_MESHS; me++) {
		solidMesh[me].initMesh(nP);
		if(me == 0)
			solidMesh[0].readMeshFile(IN_MESH_0);
#ifdef FEM_ON
		if(me == 1)
			solidMesh[1].readMeshFile(IN_MESH_1);
#endif
#ifdef FORCED_ON
		if(me == 2)
			solidMesh[2].readMeshFile(IN_MESH_2);
#endif
//		if(me == 3)
//			mesh[3].readMeshFile(IN_MESH_3);
	}

	
	// Allocation of memory
	AlcBkt();
	// Setting parameters
	SetPara();
	// Small subdomains (buckets)
	MkBkt();

	timer_sta = get_dtime();
	
	std::cout << " Number of meshes: " << NUM_MESHS << std::endl;
	for(int me = 0; me < NUM_MESHS; me++) {
  		std::cout << " after removing duplicates, mesh containts " << solidMesh[me].NV.rows() << " vertices and " << solidMesh[me].NF.rows() << " faces" << std::endl;
  	}

  	// Creation of the wall weight (Zij) and number of neighboors functions (numNeighWall)
	std::cout << " reS: " << reS << " reL: " << reL << std::endl;
	for(int me = 0; me < NUM_MESHS; me++) {
		solidMesh[me].initWijnNeigh(DIM, wijType, PCL_DST, reL, reS);
	}

#ifdef FORCED_ON
	
	nodeFRWX = (double*)malloc(sizeof(double)*solidMesh[FRW_ID].NV.rows());	// Node position
	nodeFRWY = (double*)malloc(sizeof(double)*solidMesh[FRW_ID].NV.rows());	// Node position
	nodeFRWZ = (double*)malloc(sizeof(double)*solidMesh[FRW_ID].NV.rows());	// Node position

	for(int nn=0;nn<solidMesh[FRW_ID].NV.rows();nn++) {
		nodeFRWX[nn] = solidMesh[FRW_ID].NV(nn,0);
		nodeFRWY[nn] = solidMesh[FRW_ID].NV(nn,1);
		nodeFRWZ[nn] = solidMesh[FRW_ID].NV(nn,2);
	//	std::cout << "N: " << nn << " X: " << nodeFRWX[nn] << " Y: " << nodeFRWY[nn] << " Z: " << nodeFRWZ[nn] << std::endl;
	}

	velVWall[0]=velVWall[1]=velVWall[2]=0.0;

#endif

#ifdef POLYGON_ON
	// Positions of wall and mirror particles
	for(int me = 0; me < NUM_MESHS; me++) {
		solidMesh[me].closestPointPNDBoundaryAABB(reS2, reL2, nP, wijType, Typ, FLD, GST, me, STA_ID, FEM_ID, FRW_ID, Pos, wallPos, mirrorPos, riw2, elementID, meshID, NormalWall);
//		solidMesh[me].closestPointPNDBoundaryAABB(reS2, reL2, nP, wijType, Typ, FLD, GST, Pos, wallPos, mirrorPos, riw2, niw, numNeighw, elementID, partNearMesh);
	}
	partNearMesh.clear();
	//for(int me = 0; me < NUM_MESHS; me++) {
	//	solidMesh[me].updateparticlesNearMesh(reS2, reL2, nP, wijType, Typ, FLD, GST, riw2, niw, numNeighw, partNearMesh);
	//}
	// Only call once since the distance riw2 have the values from all the meshes
	solidMesh[0].updateparticlesNearMesh(reS2, reL2, nP, wijType, Typ, FLD, GST, riw2, niw, numNeighw, partNearMesh, Nw);
#endif

#ifdef POLYGON_ON
	// NPCD PND due polygon wall
	WallNPCD_omp();
#endif
	// Initial PND
	InitPNDnNeigh_omp();
#ifdef POLYGON_ON
	// Contribution due polygon wall
	WallMeanPND_omp();
#endif
	MeanPND_omp();
	// Pressure calculation and type of particle
	if(MPS_TYP == 0)
		MkPrs_omp();
	else if(MPS_TYP == 1)
		MkPrsWc_omp();

	// PVD file
	WrtPvd();

	// Main loop
	ClcEMPS(solidMesh);

	// Deallocate memory block
	free(Acc); free(Pos); free(Vel); free(Prs);	free(pav); free(Typ); free(Normal); free(dev); free(devMod2);
	free(bfst);	free(blst);	free(nxt); free(Ci); free(gradCi); free(CMr1); free(CMr2); free(CMr3); free(DIVi);
	free(DIi); free(numNeigh); free(pndi); free(pndS); free(Bc); free(AccStar); free(NormalWall); free(devXnormal);
	free(wallPos); free(mirrorPos); free(niw); free(Nw); free(numNeighw); free(F1); free(F2); free(riw2); free(meshID);
	free(Cv); free(II); free(PTYPE); free(MEU); free(MEU_Y); free(Inertia);	free(pnew);	free(p_rheo_new); 
	free(RHO); free(p_smooth); free(VF); free(S12);	free(S13); free(S23); free(S11); free(S22); free(S33);
	free(elementID);
	free(nodeFRWX); free(nodeFRWY); free(nodeFRWZ);
	free(numNeighFree);
	//	free(Posk); free(Velk); free(Acv);

	delete[] solidMesh;

	printf("end emps_omp.\n");
	return 0;
}
