#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "polygon.h"
#include <sys/time.h>
#include <bits/stdc++.h> 
#include <iostream> 
#include <sys/stat.h> 
#include <sys/types.h> 

// Grid file
#define IN_FILE "input/dambreak_fluid.prof"
//#define IN_FILE "input/brumadinho_fluid_lo10p00.prof"
//#define IN_FILE "input/brumadinho_fluid_lo05p00_desl.prof"

// Mesh file
#define IN_MESH "input/dam1610.stl"
//#define IN_MESH "input/BRUMADINHO_space10_model_lucas_top_0p50.stl"

// Output folder
//#define OUT_FOLDER "BRUMADINHO_lo10p00_rho1500_v04p00e-01"
//#define OUT_FOLDER "BRUMADINHO_lo05p00_rho1500_v04p00e-01"
#define OUT_FOLDER "dam1610_NONEW_02"

// Output
#define OUTPND 1
#define OUTAUX 1

// Geometry
#define PCL_DST 0.01				// Average particle distance (m) (10/5) (0.01)
#define MIN_X  (0.0 - PCL_DST*3)	// Minimum value in the x direction of the analysis domain (m) (0)
#define MIN_Y  (0.0 - PCL_DST*3)	// Minimum value in the y direction of the analysis domain (m) (0)
#define MIN_Z  (0.0 - PCL_DST*3)	// Minimum value in the z direction of the analysis domain (m) (700) (0)
#define MAX_X  (1.65 + PCL_DST*3)	// Maximum value in the x direction of the analysis domain (m) (6000) (1.65)
#define MAX_Y  (0.15 + PCL_DST*3)	// Maximum value in the y direction of the analysis domain (m) (5300) (0.15)
#define MAX_Z  (0.70 + PCL_DST*30)	// Maximum value in the z direction of the analysis domain (m) (1300) (0.70)
// Model
#define DIM 3				// Dimension
#define GST -1				// Ghost particle ID
#define FLD 0				// Fluid particle ID
#define WLL 1				// Wal particle ID
#define SRF 1				// Surface particle BC
#define INR 0				// Inner particle BC
#define OTR -1				// Other particle BC
#define NUM_TYP 2			// Number of particle types	
// Physical
#define DNS_FLD 1000		// Fluid particle density (kg/m3)
#define DNS_WLL 1000		// Wall particle density (kg/m3)
#define KNM_VS1 0.000001			// Kinematic viscosity phase 1 (m2/s)
#define KNM_VS2 0.0001		// Kinematic viscosity phase 2 (m2/s)
#define DNS_FL1 1000.0		// Fluid particle density phase 1 (kg/m3)
#define DNS_FL2 2000.0		// bulk particle density phase 2 (kg/m3)
#define DNS_SDT	2000.0		// sediment density (kg/m3)
#define G_X 0.0				// Gravity acceleration x component (m/s2)
#define G_Y 0.0				// Gravity acceleration y component (m/s2)
#define G_Z -9.81			// Gravity acceleration z component (m/s2)
// Rheological parameters
#define Fluid2_type 1		// Newtonian:0  , Non Newtonian:1 
#define N 1.0				// flow behaviour (power law) index
#define MEU0 0.03			// consistency index
#define PHI 0.05			// friction angle (RAD)
#define PHI_WAL 0.005		// friction angle (RAD)
#define PHI_BED 0.005		// friction angle (RAD)
#define PHI_2 0.2			// second friction angle Values are based on  Minatti & Paris (2015)
#define cohes 0.0				// cohesiveness coefficient
#define Fraction_method 2   // Method of calculation of volume of fraction. 1: Linear dist across the interface, 2: smoothed value
//#define visc_max 20			// maximum viscosity uses to avoid singularity
#define dg 0.01			// grain size
#define I0 0.75				// I0 value in Meu9I0 rheology     Values are based on  Minatti &  Paris (2015)
#define mm 200.0
#define stress_cal_method 1	// Method 1; viscosity is directly used in momentum equation. Method 2: first the stress tensor is calculated then it is used in momentum equation
#define visc_itr_num 1
#define visc_error 0.0
#define visc_ave 0.0
#define Cd 0.47			// Drag coefficient
#define VF_min 0.25		// Minimum volume fraction
#define VF_max 0.65		// Maximum volume fraction
// Numerical
#define DT 0.00025				// Time step (s) (0.02/0.01) (0.00025)
#define FIN_TIM 1.0		// Time of simulation (s)
#define OPT_FQC 100			// Number of iterations to determine the output interval
#define MPS_TYP	0			// Explicit MPS = 0 ;  Weakly compressible MPS = 1
#define GRD_TYP	0			// Pressure gradient: Pj - Pmin = 0 ;  Pj + Pi = 1 ; Pj + Pi - 2*Pmin = 0
#define RLX_PRS 1.0			// Relaxation factor for pressure correction
#define SND 10.00			// Sound speed (m/s) (10)
#define GAM 7.0				// Gamma weakly compressible MPS
#define ARF 5000000.0			// Wall coefficent repulsive force (10000/100000) (500000)
#define SLP 1				// No-slip = 0 ; Slip = 1 
#define CRT_NUM 0.2			// Courant (CFL) condition number
#define PND_TRS 0.93			// Surface threshold PND
#define NGH_TRS 0.85		// Surface threshold Neighboors (not implemented)
#define COL_RAT 0.2			// Collision ratio
#define DST_LMT_RAT 0.85	// Coefficient of distance which does not allow any further access between particles (0.9)
#define WEI(dst, re) ((re/dst) - 1.0)	// Weight function
//#define WEI(dst, re) ((re/dst+dst/re) - 2.0)	// Weight function
//#define WEI(dst, re) (re/dst-dst/re)	// Weight function
//#define WEI(dst, re) (pow(1-dst/re,3.0))	// Weight function
#define WEI_GRAD(dst, re) ((re/dst) - 1.0)	// Weight function Gradient
//#define WEI_GRAD(dst, re) (re/dst-dst/re)	// Weight function Gradient
//#define WEI_GRAD(dst, re) (pow(1-dst/re,3.0))	// Weight function

FILE* fp;
char filename[256];
int iLP,iF;
double TIM, vMax, Crt, timer_sta, timer_end;
int nP;
double *Acc,*Pos,*Vel,*Prs,*pav;
int *Typ;
double r,r2;
double DB,DB2,DBinv;
int nBx,nBy,nBz,nBxy,nBxyz;
int *bfst,*blst,*nxt;
double n0,lmd,A1,A2,A3,A4,rlim,rlim2,COL;
double Dns[NUM_TYP],invDns[NUM_TYP];
// Polygon
double n0Grad;
double *wallPos, *mirrorPos, *niw, *pndi;
int *numNeigh, *Bc, *Nw;

double *F1, *F2;
//double *Posk, *Velk, *Acv;

// Non-Newtonian
double A1_M, NEU;
double *C, *II, *MEU, *MEU_Y, *Inertia, *pnew, *p_rheo_new, *RHO, *p_smooth, *VF;
int *PTYPE;

// Vector with ID of particles near to mesh
std::vector<int> partNearMesh;

double get_dtime(void){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 0.000001);
}

void ChkPcl(int i){
	if(	Pos[i*3  ]>MAX_X || Pos[i*3  ]<MIN_X ||
		Pos[i*3+1]>MAX_Y || Pos[i*3+1]<MIN_Y ||
		Pos[i*3+2]>MAX_Z || Pos[i*3+2]<MIN_Z)
	{
		Typ[i] = GST;
		Bc[i] = OTR;
		Nw[i] = 0;
		Prs[i]=Vel[i*3]=Vel[i*3+1]=Vel[i*3+2]=0.0;
	}
}

void RdDat(void) {
	fp = fopen(IN_FILE, "r");
	fscanf(fp,"%d",&nP);
	printf("nP: %d\n",nP);
	Acc = (double*)malloc(sizeof(double)*nP*3);	// Particle acceleration
	Pos = (double*)malloc(sizeof(double)*nP*3);	// Particle coordinates
	Vel = (double*)malloc(sizeof(double)*nP*3);	// Particle velocity
	Prs = (double*)malloc(sizeof(double)*nP);	// Particle pressure
	pav = (double*)malloc(sizeof(double)*nP);	// Time averaged particle pressure
	Typ = (int*)malloc(sizeof(int)*nP);			// Particle type

	// Polygons
	wallPos = (double*)malloc(sizeof(double)*nP*3);		// Particle at wall coordinate
	mirrorPos = (double*)malloc(sizeof(double)*nP*3);	// Mirrored particle coordinate
	niw = (double*)malloc(sizeof(double)*nP);			// PND wall
	numNeigh = (int*)malloc(sizeof(int)*nP);			// Number of neighboors
	pndi = (double*)malloc(sizeof(double)*nP);			// PND
	Bc = (int*)malloc(sizeof(int)*nP);					// BC particle type
	F1 = (double*)malloc(sizeof(double)*nP*3);		// Wall-Particle force
	F2 = (double*)malloc(sizeof(double)*nP*3);		// Wall-Particle force
	Nw = (int*)malloc(sizeof(int)*nP);					// Particle near wall

//	Posk = (double*)malloc(sizeof(double)*nP*3);	// Particle coordinates
//	Velk = (double*)malloc(sizeof(double)*nP*3);	// Particle velocity
//	Acv = (double*)malloc(sizeof(double)*nP*3);	// Part

	// Non-Newtonian
	C = (double*)malloc(sizeof(double)*nP);				// Concentration
	II = (double*)malloc(sizeof(double)*nP);			// Invariant
	PTYPE = (int*)malloc(sizeof(int)*nP);				// Type of fluid
	MEU = (double*)malloc(sizeof(double)*nP);			// Dynamic viscosity
	MEU_Y = (double*)malloc(sizeof(double)*nP);			// Dynamic viscosity ??
	Inertia = (double*)malloc(sizeof(double)*nP);		//
	pnew = (double*)malloc(sizeof(double)*nP);			// New pressure
	p_rheo_new = (double*)malloc(sizeof(double)*nP);	//
	RHO = (double*)malloc(sizeof(double)*nP);			// Fluid density
	p_smooth = (double*)malloc(sizeof(double)*nP);	//
	VF = (double*)malloc(sizeof(double)*nP);	//

	for(int i=0;i<nP;i++) {
		int a[2];
		double b[8];
		fscanf(fp," %d %d %lf %lf %lf %lf %lf %lf %lf %lf",&a[0],&a[1],&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7]);
		Typ[i]=a[1];
		Pos[i*3]=b[0];	Pos[i*3+1]=b[1];	Pos[i*3+2]=b[2];
		Vel[i*3]=b[3];	Vel[i*3+1]=b[4];	Vel[i*3+2]=b[5];
		Prs[i]=b[6];		pav[i]=b[7];

//		Posk[i*3]=b[0];	Posk[i*3+1]=b[1];	Posk[i*3+2]=b[2];
//		Velk[i*3]=b[3];	Velk[i*3+1]=b[4];	Velk[i*3+2]=b[5];
	}
	fclose(fp);
	for(int i=0;i<nP;i++) {ChkPcl(i);}
	for(int i=0;i<nP*3;i++) {Acc[i]=0.0;/*Acv[i]=0.0;*/}
	for(int i=0;i<nP*3;i++) {wallPos[i]=0.0;mirrorPos[i]=0.0;F1[i]=0.0;F2[i]=0.0;}
	for(int i=0;i<nP;i++) {
		niw[i]=0.0;numNeigh[i]=0.0;pndi[i]=0.0;Bc[i]=0;Nw[i]=0;
		C[i]=0.0;II[i]=0.0;PTYPE[i]=0;MEU[i]=0.0;MEU_Y[i]=0.0;Inertia[i]=0.0;pnew[i]=0.0;p_rheo_new[i]=0.0;p_smooth[i]=0.0;VF[i]=0.0;

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
	}
}
/*
void RdMes(void) {
	fp = fopen(IN_FILE, "r");
	fscanf(fp,"%d",&nP);
	mesh mesh(IN_FILE, nP);
	mesh.readMeshFile(IN_MESH);
}
*/
void WrtDat(void) {
	char outout_filename[256];
	sprintf(outout_filename, "output%05d.prof",iF);
	fp = fopen(outout_filename, "w");
	fprintf(fp,"%d\n",nP);
	for(int i=0;i<nP;i++) {
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

void WrtVtu(void) {
	char outout_filename[256];
	sprintf(outout_filename, OUT_FOLDER"/output%05d.vtu",iF);
	fp = fopen(outout_filename, "w");
	fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
	fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
	fprintf(fp,"<UnstructuredGrid>\n");
 	fprintf(fp,"<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",nP,nP);

 	fprintf(fp,"<Points>\n");
 	fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
	for(int i=0;i<nP;i++)fprintf(fp,"%lf %lf %lf ",Pos[i*3],Pos[i*3+1],Pos[i*3+2]);
 	fprintf(fp,"\n</DataArray>\n");
  	fprintf(fp,"</Points>\n");

  	fprintf(fp,"<PointData>\n");
 	fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Velocity' format='ascii'>\n");
	for(int i=0;i<nP;i++){
		//double val=sqrt(Vel[i*3]*Vel[i*3]+Vel[i*3+1]*Vel[i*3+1]+Vel[i*3+2]*Vel[i*3+2]);
		fprintf(fp,"%f %f %f ",(float)Vel[i*3],(float)Vel[i*3+1],(float)Vel[i*3+2]);
	}
 	fprintf(fp,"\n</DataArray>\n");
 	//fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='pressave' format='ascii'>\n");
	//for(int i=0;i<nP;i++){	fprintf(fp,"%f ",(float)(pav[i]/OPT_FQC));}
	//fprintf(fp,"\n</DataArray>\n");
	fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='pressure' format='ascii'>\n");
	for(int i=0;i<nP;i++){	fprintf(fp,"%f ",(float)Prs[i]);}
	fprintf(fp,"\n</DataArray>\n");
	if(OUTPND){
 		fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='pnd' format='ascii'>\n");
		for(int i=0;i<nP;i++){	fprintf(fp,"%f\n",(float)pndi[i]);}
 		fprintf(fp,"\n</DataArray>\n");
 	}

 	/*
	fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='F1' format='ascii'>\n");
	for(int i=0;i<nP;i++){
		fprintf(fp,"%f %f %f ",(float)F1[i*3],(float)F1[i*3+1],(float)F1[i*3+2]);
	}
 	fprintf(fp,"\n</DataArray>\n");

 	fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='F2' format='ascii'>\n");
	for(int i=0;i<nP;i++){
		fprintf(fp,"%f %f %f ",(float)F2[i*3],(float)F2[i*3+1],(float)F2[i*3+2]);
	}
 	fprintf(fp,"\n</DataArray>\n");
 	*/
 	fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='RHO' format='ascii'>\n");
		for(int i=0;i<nP;i++){	fprintf(fp,"%f\n",(float)RHO[i]);}
 		fprintf(fp,"\n</DataArray>\n");
 	fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Conc' format='ascii'>\n");
		for(int i=0;i<nP;i++){	fprintf(fp,"%f\n",(float)C[i]);}
 		fprintf(fp,"\n</DataArray>\n");
 	fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='MEU' format='ascii'>\n");
		for(int i=0;i<nP;i++){	fprintf(fp,"%f\n",(float)MEU[i]);}
 		fprintf(fp,"\n</DataArray>\n");
 	fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='MEUy' format='ascii'>\n");
		for(int i=0;i<nP;i++){	fprintf(fp,"%f\n",(float)MEU_Y[i]);}
 		fprintf(fp,"\n</DataArray>\n");
 	fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='II' format='ascii'>\n");
		for(int i=0;i<nP;i++){	fprintf(fp,"%f\n",(float)II[i]);}
 		fprintf(fp,"\n</DataArray>\n");
 	fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='PTYPE' format='ascii'>\n");
		for(int i=0;i<nP;i++){	fprintf(fp,"%d ",PTYPE[i]);}
 		fprintf(fp,"\n</DataArray>\n");
 	fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='p_smooth' format='ascii'>\n");
	for(int i=0;i<nP;i++){	fprintf(fp,"%f ",(float)p_smooth[i]);}
	fprintf(fp,"\n</DataArray>\n");

  	if(OUTAUX){
 		fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
		for(int i=0;i<nP;i++){fprintf(fp,"%d ",Typ[i]);}
 		fprintf(fp,"\n</DataArray>\n");
 		fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='BC' format='ascii'>\n");
		for(int i=0;i<nP;i++){fprintf(fp,"%d ",Bc[i]);}
 		fprintf(fp,"\n</DataArray>\n");
 		//fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='Nwall' format='ascii'>\n");
		//for(int i=0;i<nP;i++){fprintf(fp,"%d ",Nw[i]);}
 		//fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Fwall' format='ascii'>\n");
		//for(int i=0;i<nP;i++){
		//	fprintf(fp,"%f %f %f ",(float)Fwall[i*3],(float)Fwall[i*3+1],(float)Fwall[i*3+2]);
		//}
 		//fprintf(fp,"\n</DataArray>\n");
 	}
  	fprintf(fp,"</PointData>\n");

 	fprintf(fp,"<Cells>\n");
 	fprintf(fp,"<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
	for(int i=0;i<nP;i++)fprintf(fp,"%d ",i);
 	fprintf(fp,"\n</DataArray>\n");
 	fprintf(fp,"<DataArray type='Int32' Name='offsets' format='ascii'>\n");
	for(int i=0;i<nP;i++)fprintf(fp,"%d ",i+1);
 	fprintf(fp,"\n</DataArray>\n");
 	fprintf(fp,"<DataArray type='UInt8' Name='types' format='ascii'>\n");
	for(int i=0;i<nP;i++)fprintf(fp,"1 ");
 	fprintf(fp,"\n</DataArray>\n");
 	fprintf(fp,"</Cells>\n");
 	fprintf(fp,"</Piece>\n");
 	fprintf(fp,"</UnstructuredGrid>\n");
 	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
	iF++;
}

void WrtPvd(void) {
	char pvd_filename[256];
	sprintf(pvd_filename, OUT_FOLDER".pvd");
	fp = fopen(pvd_filename, "w");
	int nIter = ceil(FIN_TIM/DT);

	//fprintf(fp,"<VTKFile type=""Collection"" version=""0.1"" byte_order=""LittleEndian"">\n");
	fprintf(fp,"<VTKFile type='Collection' version='0.1' byte_order='LittleEndian'>\n");
	fprintf(fp,"	<Collection>\n");
	int j = 0;
	for(int i=0;i<nIter;i++){
		if(i % OPT_FQC == 0){
			double timeStep = DT*i;
 			fprintf(fp,"		<DataSet timestep='%.6f' group='A' part='0' file='%s/output%05d.vtu'/>\n",timeStep,OUT_FOLDER,j);
 			j++;
		}
	}
	fprintf(fp,"	</Collection>\n");
 	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
}

void AlcBkt(void) {
	r = PCL_DST*2.1;							// Influence radius
	r2 = r*r;
	DB = r*(1.0+CRT_NUM);						// Length of one bucket side
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

void SetPara(void){
	n0 = n0Grad = lmd = 0.0;
	for(int ix= -4;ix<5;ix++){
	for(int iy= -4;iy<5;iy++){
	for(int iz= -4;iz<5;iz++){
		double x = PCL_DST* (double)ix;
		double y = PCL_DST* (double)iy;
		double z = PCL_DST* (double)iz;
		double dst2 = x*x+y*y+z*z;
		if(dst2 <= r2){
			if(dst2==0.0)continue;
			double dst = sqrt(dst2);
			n0 += WEI(dst, r);			// Initial particle number density
			n0Grad += WEI_GRAD(dst, r);	// Initial particle number density to Gradient
			lmd += dst2 * WEI(dst, r);
		}
	}}}
	lmd = lmd/n0;					// Coefficient λ of Laplacian model
	A1 = 2.0*KNM_VS1*DIM/(n0*lmd);	// Coefficient used to calculate viscosity term
	A1_M = 2.0*DIM/(n0*lmd);		// Coefficient used to calculate viscosity term Multiphase
	A2 = SND*SND/n0;				// Coefficient used to calculate pressure E-MPS
	A3 = -DIM/n0Grad;				// Coefficient used to calculate pressure gradient term
	A4 = SND*SND;					// Coefficient used to calculate pressure WC-MPS
	Dns[FLD]=DNS_FLD;			Dns[WLL]=DNS_WLL;
	invDns[FLD]=1.0/DNS_FLD;	invDns[WLL]=1.0/DNS_WLL;
	rlim = PCL_DST * DST_LMT_RAT;	// A distance that does not allow further access between particles
	rlim2 = rlim*rlim;
	COL = 1.0 + COL_RAT;
	iLP=iF=0;						// Number of iterations // File number
	TIM=0.0;						// Simulation time

	std::cout << "lo: " << PCL_DST << " m, dt: " << DT << " s, PND0: " << n0 << " PND0Grad: " << n0Grad << " lambda: " << lmd << std::endl;
}

void MkBkt(void) {
	for(int i=0;i< nBxyz ;i++){	bfst[i] = -1;	}
	for(int i=0;i< nBxyz ;i++){	blst[i] = -1;	}
	for(int i=0;i< nP ;i++){	nxt[i] = -1;	}
	for(int i=0;i<nP;i++){
		if(Typ[i] == GST)continue;
		int ix = (int)((Pos[i*3  ] - MIN_X)*DBinv) +1;
		int iy = (int)((Pos[i*3+1] - MIN_Y)*DBinv) +1;
		int iz = (int)((Pos[i*3+2] - MIN_Z)*DBinv) +1;
		int ib = iz*nBxy + iy*nBx + ix;
		int j = blst[ib];
		blst[ib] = i;
		if(j == -1){	bfst[ib] = i;	}
		else{			nxt[j] = i;}
	}
}

// Calculation of the volume of fraction if phase II in the mixture
void VolFract_omp()
{
	if (Fraction_method == 1){   //Linear distribution
#pragma omp parallel for schedule(dynamic,64)
		for(int i=0;i<nP;i++){
		if(Typ[i] == FLD){
			double sum1 = 0, sum2 = 0;
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
			int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
			int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
			for(int jz=iz-1;jz<=iz+1;jz++){
			for(int jy=iy-1;jy<=iy+1;jy++){
			for(int jx=ix-1;jx<=ix+1;jx++){
				int jb = jz*nBxy + jy*nBx + jx;
				int j = bfst[jb];
				if(j == -1) continue;
				for(;;){
					double v0 = Pos[j*3  ] - pos_ix;
					double v1 = Pos[j*3+1] - pos_iy;
					double v2 = Pos[j*3+2] - pos_iz;
					double dst2 = v0*v0+v1*v1+v2*v2;
					if(dst2<r2){
					if(j!=i && Typ[j]==FLD){
						sum1 = sum1 + 1;
						if (PTYPE[j] >= 2) sum2 = sum2 + 1;
						}}
					j = nxt[j];
					if(j==-1) break;
				}
			}}}
			if(sum1 == 0)
				C[i] = 0.0;
			else 
				C[i] = sum2 / sum1;
		}}
	}
	else if (Fraction_method == 2){   //Non linear :  Smoothed using the weight funtion
#pragma omp parallel for schedule(dynamic,64)
		for(int i=0;i<nP;i++){
		if(Typ[i] == FLD){
			double sum1 = 0, sum2 = 0;
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
			int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
			int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
			for(int jz=iz-1;jz<=iz+1;jz++){
			for(int jy=iy-1;jy<=iy+1;jy++){
			for(int jx=ix-1;jx<=ix+1;jx++){
				int jb = jz*nBxy + jy*nBx + jx;
				int j = bfst[jb];
				if(j == -1) continue;
				for(;;){
					double v0 = Pos[j*3  ] - pos_ix;
					double v1 = Pos[j*3+1] - pos_iy;
					double v2 = Pos[j*3+2] - pos_iz;
					double dst2 = v0*v0+v1*v1+v2*v2;
					if(dst2<r2){
					if(j!=i && Typ[j]==FLD){
						double dst = sqrt(dst2);
						double w = WEI(dst, r);
						sum1 = sum1 + w;
						if(PTYPE[j] >= 2) sum2 = sum2 + w;
						}}
					j = nxt[j];
					if(j==-1) break;
				}
			}}}
			if(sum1 == 0)
				C[i] = 0.0;
			else 
				C[i] = sum2 / sum1;
		}}
	}
}

//void NonNwtVscTrm_omp(double *x_vel, double *y_vel, double *z_vel){
// Viscosity interaction values for "real" fluid particles
void VscIntVal_omp(){

	double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu0, normal_stress;//,grain_VF, *p_smooth;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
	Xmin = 0.0 - PCL_DST*3; Xmax = 1.65 + PCL_DST*3;
	Ymin = 0.0 - PCL_DST*3; Ymax = 0.15 + PCL_DST*3;
	Zmin = 0.0 - PCL_DST*3; //Zmax = 0.7 + PCL_DST*30;
	// Changed !!!

	// Search free-surface particles for each interval of aa = 2 particles in wall
	int aa = 2, kx, ky;
	int kx_max = int((Xmax - Xmin) / aa / PCL_DST) + 1;
	int ky_max = int((Ymax - Ymin) / aa / PCL_DST) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;

	S11 = new double[nP + 1];
	S22 = new double[nP + 1];
	S33 = new double[nP + 1];
	S12 = new double[nP + 1];
	S13 = new double[nP + 1];
	S23 = new double[nP + 1];
	//p_smooth = new double[nP + 1];
	BL = new double*[kx_max + 1];  // bed level
	WL = new double*[kx_max + 1];  // water level
	PS = new double*[kx_max + 1];  // pressure sediment

#pragma omp parallel for
	for (int m = 1; m <= kx_max; m++)
	{
		BL[m] = new double[ky_max + 1];
		WL[m] = new double[ky_max + 1];
		PS[m] = new double[ky_max + 1];
	}

	// Determining the bed level
#pragma omp parallel for schedule(dynamic,64)
	for (kx = 1; kx <= kx_max; kx++)
	{
		for (ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0;i<nP;i++){
		if(Typ[i] == FLD){
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		
			kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
			ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;


			//if (pos_iz>BL[kx][ky] && C[i]>0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = pnew[i]; }
			if (pos_iz>BL[kx][ky] && C[i]>0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = Prs[i]; }
			if (pos_iz>WL[kx][ky] && PTYPE[i] == 1) { WL[kx][ky] = pos_iz; }
		}
	}

	// Strain rate calculation
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0;i<nP;i++){
	if(Typ[i] == FLD){
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;

//		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2<r2){
				if(j!=i && Typ[j]!=GST){
					double dst = sqrt(dst2);
					double w = WEI(dst, r);
					double vec_ijx = Vel[j*3  ]-vec_ix;	
					double vec_ijy = Vel[j*3+1]-vec_iy;	
					double vec_ijz = Vel[j*3+2]-vec_iz;

					sum1 += vec_ijx*v0*w/dst2;
					sum2 += vec_ijx*v1*w/dst2;
					sum3 += vec_ijx*v2*w/dst2;
					
					sum4 += vec_ijy*v0*w/dst2;
					sum5 += vec_ijy*v1*w/dst2;
					sum6 += vec_ijy*v2*w/dst2;
					
					sum7 += vec_ijz*v0*w/dst2;
					sum8 += vec_ijz*v1*w/dst2;
					sum9 += vec_ijz*v2*w/dst2;
					
					sum10 += Prs[j]*w;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}

		// A3 is a negative cte (-DIM/n0Grad)
		Uxx = -A3*sum1; Uxy = -A3*sum2; Uxz = -A3*sum3;
		Uyx = -A3*sum4; Uyy = -A3*sum5; Uyz = -A3*sum6;
		Uzx = -A3*sum7; Uzy = -A3*sum8; Uzz = -A3*sum9;

		p_smooth[i] = sum10 / n0Grad;
		if (p_smooth[i]<0) p_smooth[i] = 0;

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
		if (II[i]<0 || II[i] * 0 != 0) II[i] = 0;
		//II=fabs(S11[i]*S22[i]-S12[i]*S12[i]);
		
		//std::cout << " II: " << II[i] << std::endl;
	}}

	// Newtonian viscosity
	if (Fluid2_type == 0)
	{
#pragma omp parallel for
		for(int i=0;i<nP;i++){
		if(Typ[i] == FLD){
			if (PTYPE[i] <= 1)MEU[i] = KNM_VS1 * DNS_FL1;
			if (PTYPE[i] != 1)MEU[i] = KNM_VS2 * DNS_FL2;
		}}

//		if (TURB>0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2*sqrt(II[i]);

//			if (NEUt[i] * 0 != 0)  NEUt[i] = 0;
//			if (NEUt[i]>1)     NEUt[i] = 1;
//		}
	}

	// Granular Fluid
	if (Fluid2_type == 1)
	{
	// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
//#pragma omp parallel for
		for(int i=0;i<nP;i++){
		if(Typ[i] == FLD){

//			double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];

			if(PTYPE[i] == 1) {
				MEU[i] = KNM_VS1 * DNS_FL1;
			}
			else if(PTYPE[i] == 2){
				// phi: internal friction angle
				// phi2: ?
				phi = (C[i] - 0.25)*PHI / (1 - 0.25);
				phi2 = (C[i] - 0.25)*PHI_2 / (1 - 0.25);
				if (C[i] <= 0.25) { phi = 0.00001; phi2 = 0.00001; } // phi close to zero
				if (PTYPE[i] <= 0) phi = PHI_BED; // ghost

				// Normal stress calculation (mechanical pressure)
				p_rheo_new[i] = p_smooth[i];

				kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
				ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress=(BL[k]-pos_iy+DL/2)*(DNS_FL2)*9.81;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - pos_iz + PCL_DST / 2)*(DNS_FL2 - DNS_FL1)*9.81 - (vec_ix*vec_ix + vec_iy*vec_iy + vec_iz*vec_iz)*(DNS_FL2 - DNS_FL1) / 2.0;	// normal_stress= Gama.H

				if (p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81<0) p_smooth[i] = (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
				if (TIM <= 1) normal_stress = 1.0*(1 - TIM)*(p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81) + 1.0*(TIM)*normal_stress;

//				normal_stress = p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
//				normal_stress = p_smooth[i];
				//normal_stress=normal_stress*0.61*1500/DNS_FL2;
				if (normal_stress < 1 || C[i] < 0.5) normal_stress = 1;

				p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Inertia[i] = sqrt(II[i])*dg/sqrt(normal_stress/DNS_SDT);		// Free-fall (dry granular material)
				Inertia[i] = sqrt(II[i])*dg/sqrt(normal_stress/(DNS_FL1*Cd));	// Grain inertia (submerged)
				//Inertia[i] = sqrt(II[i])*(KNM_VS1*DNS_FL1)/normal_stress ;	// Viscous regime

//				Inertia[i] = 1.0;
				// VF_max VF_min
				VF[i] = VF_max - (VF_max - VF_min)*Inertia[i];
				if (VF[i] < VF_min) VF[i] = VF_min;
				RHO[i] = DNS_SDT * VF[i] + (1-VF[i])*DNS_FL1;
				phi = phi * VF[i] / VF_max;

				// Mohr-Coulomb
				double yield_stress = cohes * cos(phi) + normal_stress * sin(phi);

				if (yield_stress < 0) yield_stress = 0;

				double visc_max = (yield_stress*mm*0.5 + MEU0);

				if (II[i]>0)
					MEU_Y[i] = yield_stress * (1 - exp(-mm * sqrt(II[i]))) / 2.0 / sqrt(II[i]);
				else
					MEU_Y[i] = visc_max;

				// H-B rheology

				//meu0 = MEU0;

				// Non-linear Meu(I) rheology
				//meu0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*dg/(I0*sqrt(normal_stress/DNS_FL2)+sqrt(II[i])*dg);			//free fall
				meu0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*dg/(I0*sqrt(normal_stress/(DNS_FL1*Cd))+sqrt(II[i])*dg);		//grain inertia
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(KNM_VS1*DNS_FL1)/(I0*normal_stress+sqrt(II[i])*(KNM_VS1*DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*dg*sqrt(normal_stress*DNS_FL2)/I0;		//free fall
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*dg*sqrt(normal_stress*DNS_FL1*Cd)/I0;	//grain inertia
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*(KNM_VS1*DNS_FL1)/I0;					//viscous

				if (II[i] <= 0 || (meu0 * 0) != 0) meu0 = MEU0;

				visc_max = (yield_stress*mm*0.5 + meu0);

				//if(isnan(II[i]) || isinf(II[i])){
					//std::cout << " viscmax: " << II[i] << std::endl;
				//	assert(visc_max >= 0 || visc_max <= 0);
				//}
				
				// Herschel bulkley papanastasiou
				MEU[i] = MEU_Y[i] + MEU0 * pow(4 * II[i], (N - 1) / 2);

				// MEU_Y rheological model
				//MEU[i] = MEU_Y[i] + meu0;
				
				if (II[i] == 0 || MEU[i]>visc_max) MEU[i] = visc_max;
				if (PTYPE[i] <= 0) MEU[i] = MEU[i] * C[i] + DNS_FL1*KNM_VS1*(1 - C[i]);

				//if (MEU[i]/RHO[i] > maxVIS) maxVIS = MEU[i]/RHO[i];
			}
			
			if(PTYPE[i] >= 2) {
				if (C[i] > 0.5) RHO[i] = DNS_FL2;
				else RHO[i] = C[i] * DNS_FL2 + (1 - C[i]) * DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if (stress_cal_method == 2)
//		{
//			for (i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
//				for (l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if (i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2 * MEU[i] * MEU[j] / (MEU[i] + MEU[j]);
//						if ((NEUt[i] + NEUt[j])>0) meuij = meuij + 2 * NEUt[i] * RHO[i] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j]);

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

	} // if (Fluid2_type == 1)

	//---------------------------------------------------------------

	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS; //delete[]p_smooth;
	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL; //p_smooth = NULL;
}

// Slip condition. Viscosity interaction values
void WallSlipVscIntVal_omp(){

	double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu0, normal_stress;//, grain_VF, *p_smooth;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
	Xmin = 0.0 - PCL_DST*3; Xmax = 1.65 + PCL_DST*3;
	Ymin = 0.0 - PCL_DST*3; Ymax = 0.15 + PCL_DST*3;
	Zmin = 0.0 - PCL_DST*3; //Zmax = 0.7 + PCL_DST*30;
	// Changed !!!

	// Search free-surface particles for each interval of aa = 2 particles in wall
	int aa = 2, kx, ky;
	int kx_max = int((Xmax - Xmin) / aa / PCL_DST) + 1;
	int ky_max = int((Ymax - Ymin) / aa / PCL_DST) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
	double aUxx, aUxy, aUxz, aUyx, aUyy, aUyz, aUzx, aUzy, aUzz;

	S11 = new double[nP + 1];
	S22 = new double[nP + 1];
	S33 = new double[nP + 1];
	S12 = new double[nP + 1];
	S13 = new double[nP + 1];
	S23 = new double[nP + 1];
	//p_smooth = new double[nP + 1];
	BL = new double*[kx_max + 1];  // bed level
	WL = new double*[kx_max + 1];  // water level
	PS = new double*[kx_max + 1];  // pressure sediment

#pragma omp parallel for
	for (int m = 1; m <= kx_max; m++)
	{
		BL[m] = new double[ky_max + 1];
		WL[m] = new double[ky_max + 1];
		PS[m] = new double[ky_max + 1];
	}

	// Determining the bed level
#pragma omp parallel for schedule(dynamic,64)
	for (kx = 1; kx <= kx_max; kx++)
	{
		for (ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0;i<nP;i++){
		if(Typ[i] == FLD){
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			
			kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
			ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;

			//if (pos_iz>BL[kx][ky] && C[i]>0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = pnew[i]; }
			if (pos_iz>BL[kx][ky] && C[i]>0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = Prs[i]; }
			if (pos_iz>WL[kx][ky] && PTYPE[i] == 1) { WL[kx][ky] = pos_iz; }
		}
	}

	// Strain rate calculation
	int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	for(int im=0;im<nPartNearMesh;im++){
	//for(int i=0;i<nP;i++){
	int i = partNearMesh[im];
	if(Typ[i] == FLD){
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

	    if (normalMod2 > 0.00000001) {
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
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0 = Pos[j*3  ] - pos_mix;
				double v1 = Pos[j*3+1] - pos_miy;
				double v2 = Pos[j*3+2] - pos_miz;

				double dst2 = v0*v0+v1*v1+v2*v2;
				// If inside neighboor of i and im (intersection)
				if(dstij2<r2 && dst2<r2){
				if(j!=i && Typ[j]!=GST){
					double dst = sqrt(dst2);
					double w = WEI(dst, r);

					double vec_mijx = Vel[j*3  ]-vec_mix;	
					double vec_mijy = Vel[j*3+1]-vec_miy;	
					double vec_mijz = Vel[j*3+2]-vec_miz;

					sum1 += vec_mijx*v0*w/dst2;
					sum2 += vec_mijx*v1*w/dst2;
					sum3 += vec_mijx*v2*w/dst2;
					
					sum4 += vec_mijy*v0*w/dst2;
					sum5 += vec_mijy*v1*w/dst2;
					sum6 += vec_mijy*v2*w/dst2;
					
					sum7 += vec_mijz*v0*w/dst2;
					sum8 += vec_mijz*v1*w/dst2;
					sum9 += vec_mijz*v2*w/dst2;
					
					sum10 += Prs[j]*w;
				}}
				j = nxt[j];
				if(j==-1) break;
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
		if (p_smooth[i]<0) p_smooth[i] = 0;

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
		if (II[i]<0 || II[i] * 0 != 0) II[i] = 0;
		//II=fabs(S11[i]*S22[i]-S12[i]*S12[i]);
	}}
	
	// Newtonian viscosity
	if (Fluid2_type == 0)
	{
	// Loop only for particles near mesh
#pragma omp parallel for
		for(int im=0;im<nPartNearMesh;im++){
		//for(int i=0;i<nP;i++){
		int i = partNearMesh[im];
		if(Typ[i] == FLD){
			if (PTYPE[i] <= 1)MEU[i] = KNM_VS1 * DNS_FL1;
			if (PTYPE[i] != 1)MEU[i] = KNM_VS2 * DNS_FL2;
		}}

//		if (TURB>0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2*sqrt(II[i]);

//			if (NEUt[i] * 0 != 0)  NEUt[i] = 0;
//			if (NEUt[i]>1)     NEUt[i] = 1;
//		}
	}

	// Granular Fluid
	if (Fluid2_type == 1)
	{
		// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
		// Loop only for particles near mesh
//#pragma omp parallel for
		for(int im=0;im<nPartNearMesh;im++){
		//for(int i=0;i<nP;i++){
		int i = partNearMesh[im];
		if(Typ[i] == FLD){

//			double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];

			if(PTYPE[i] == 1){
				MEU[i] = KNM_VS1 * DNS_FL1;
			}
			else if(PTYPE[i] == 2){
				phi = (C[i] - 0.25)*PHI / (1 - 0.25);
				phi2 = (C[i] - 0.25)*PHI_2 / (1 - 0.25);
				if (C[i] <= 0.25) { phi = 0.00001; phi2 = 0.00001; } // phi close to zero
				if (PTYPE[i] <= 0) phi = PHI_BED;

				// Normal stress calculation (mehcanical pressure)
				p_rheo_new[i] = p_smooth[i];

				kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
				ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress=(BL[k]-pos_iy+DL/2)*(DNS_FL2)*9.81;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - pos_iz + PCL_DST / 2)*(DNS_FL2 - DNS_FL1)*9.81 - (vec_ix*vec_ix + vec_iy*vec_iy + vec_iz*vec_iz)*(DNS_FL2 - DNS_FL1) / 2.0;	// normal_stress= Gama.H

				if (p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81<0) p_smooth[i] = (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
				if (TIM <= 1) normal_stress = 1.0*(1 - TIM)*(p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81) + 1.0*(TIM)*normal_stress;

//				normal_stress = p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
//				normal_stress = p_smooth[i];
				//normal_stress=normal_stress*0.61*1500/DNS_FL2;
				if (normal_stress < 1 || C[i] < 0.5) normal_stress = 1;

				p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Inertia[i] = sqrt(II[i])*dg/sqrt(normal_stress/DNS_SDT);		// Free-fall (dry granular material)
				Inertia[i] = sqrt(II[i])*dg/sqrt(normal_stress/(DNS_FL1*Cd));	// Grain inertia (submerged)
				//Inertia[i] = sqrt(II[i])*(KNM_VS1*DNS_FL1)/normal_stress ;	// Viscous regime

				// VF_max VF_min
				VF[i] = VF_max - (VF_max - VF_min)*Inertia[i];
				if (VF[i] < VF_min) VF[i] = VF_min;
				RHO[i] = DNS_SDT * VF[i] + (1-VF[i])*DNS_FL1;
				phi = phi * VF[i] / VF_max;

				double yield_stress = cohes * cos(phi) + normal_stress * sin(phi);

				if (yield_stress < 0) yield_stress = 0;

				double visc_max = (yield_stress*mm*0.5 + MEU0);

				if (II[i]>0)
					MEU_Y[i] = yield_stress * (1 - exp(-mm * sqrt(II[i]))) / 2.0 / sqrt(II[i]);
				else
					MEU_Y[i] = visc_max;

				// H-B rheology

				//meu0 = MEU0;

				// Non-linear Meu(I) rheology
				//meu0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*dg/(I0*sqrt(normal_stress/DNS_FL2)+sqrt(II[i])*dg);			//free fall
				meu0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*dg/(I0*sqrt(normal_stress/(DNS_FL1*Cd))+sqrt(II[i])*dg);		//grain inertia
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(KNM_VS1*DNS_FL1)/(I0*normal_stress+sqrt(II[i])*(KNM_VS1*DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*dg*sqrt(normal_stress*DNS_FL2)/I0;		//free fall
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*dg*sqrt(normal_stress*DNS_FL1*Cd)/I0;	//grain inertia
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*(KNM_VS1*DNS_FL1)/I0;					//viscous

				if (II[i] <= 0 || (meu0 * 0) != 0) meu0 = MEU0;

				visc_max = (yield_stress*mm*0.5 + meu0);

				// Herschel bulkley papanastasiou
				MEU[i] = MEU_Y[i] + MEU0 * pow(4 * II[i], (N - 1) / 2);

				// MEU_Y rheological model
				//MEU[i] = MEU_Y[i] + meu0;
				
				if (II[i] == 0 || MEU[i]>visc_max) MEU[i] = visc_max;
				if (PTYPE[i] <= 0) MEU[i] = MEU[i] * C[i] + DNS_FL1*KNM_VS1*(1 - C[i]);
			}
			
			if(PTYPE[i] >= 2) {
				if (C[i] > 0.5) RHO[i] = DNS_FL2;
				else RHO[i] = C[i] * DNS_FL2 + (1 - C[i]) * DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if (stress_cal_method == 2)
//		{
//			for (i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
//				for (l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if (i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2 * MEU[i] * MEU[j] / (MEU[i] + MEU[j]);
//						if ((NEUt[i] + NEUt[j])>0) meuij = meuij + 2 * NEUt[i] * RHO[i] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j]);

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

	} // if (Fluid2_type == 1)

	//---------------------------------------------------------------

	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS;// delete[]p_smooth;
	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL;// p_smooth = NULL;
}

// No-Slip condition. Viscosity interaction values
void WallNoSlipVscIntVal_omp(){

	double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0.0, phi2 = 0.0, meu0, normal_stress;//, grain_VF, *p_smooth;
	double **BL, **WL, **PS;

	// Changed !!!
	// Be carefull to assign all domain
	double Xmin, Xmax, Ymin, Ymax, Zmin; // Minimum and maximum of searching grid
	// dam1610
	Xmin = 0.0 - PCL_DST*3; Xmax = 1.65 + PCL_DST*3;
	Ymin = 0.0 - PCL_DST*3; Ymax = 0.15 + PCL_DST*3;
	Zmin = 0.0 - PCL_DST*3; //Zmax = 0.7 + PCL_DST*30;
	// Changed !!!

	// Search free-surface particles for each interval of aa = 2 particles in wall
	int aa = 2, kx, ky;
	int kx_max = int((Xmax - Xmin) / aa / PCL_DST) + 1;
	int ky_max = int((Ymax - Ymin) / aa / PCL_DST) + 1;

	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;

	S11 = new double[nP + 1];
	S22 = new double[nP + 1];
	S33 = new double[nP + 1];
	S12 = new double[nP + 1];
	S13 = new double[nP + 1];
	S23 = new double[nP + 1];
	//p_smooth = new double[nP + 1];
	BL = new double*[kx_max + 1];  // bed level
	WL = new double*[kx_max + 1];  // water level
	PS = new double*[kx_max + 1];  // pressure sediment

#pragma omp parallel for
	for (int m = 1; m <= kx_max; m++)
	{
		BL[m] = new double[ky_max + 1];
		WL[m] = new double[ky_max + 1];
		PS[m] = new double[ky_max + 1];
	}

	// Determining the bed level
#pragma omp parallel for schedule(dynamic,64)
	for (kx = 1; kx <= kx_max; kx++)
	{
		for (ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}
	}

#pragma omp parallel for
	for(int i=0;i<nP;i++){
		if(Typ[i] == FLD){
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		
			kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
			ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;

			//if (pos_iz>BL[kx][ky] && C[i]>0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = pnew[i]; }
			if (pos_iz>BL[kx][ky] && C[i]>0.5) { BL[kx][ky] = pos_iz; PS[kx][ky] = Prs[i]; }
			if (pos_iz>WL[kx][ky] && PTYPE[i] == 1) { WL[kx][ky] = pos_iz; }
		}
	}

	// Strain rate calculation
	int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	for(int im=0;im<nPartNearMesh;im++){
	//for(int i=0;i<nP;i++){
	int i = partNearMesh[im];
	if(Typ[i] == FLD){
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

	    if (normalMod2 > 0.00000001) {
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
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0 = Pos[j*3  ] - pos_mix;
				double v1 = Pos[j*3+1] - pos_miy;
				double v2 = Pos[j*3+2] - pos_miz;

				double dst2 = v0*v0+v1*v1+v2*v2;
				// If inside neighboor of i and im (intersection)
				if(dstij2<r2 && dst2<r2){
				if(j!=i && Typ[j]!=GST){
					double dst = sqrt(dst2);
					double w = WEI(dst, r);

					double vec_mijx = Vel[j*3  ]-vec_mix;	
					double vec_mijy = Vel[j*3+1]-vec_miy;	
					double vec_mijz = Vel[j*3+2]-vec_miz;

					sum1 += vec_mijx*v0*w/dst2;
					sum2 += vec_mijx*v1*w/dst2;
					sum3 += vec_mijx*v2*w/dst2;
					
					sum4 += vec_mijy*v0*w/dst2;
					sum5 += vec_mijy*v1*w/dst2;
					sum6 += vec_mijy*v2*w/dst2;
					
					sum7 += vec_mijz*v0*w/dst2;
					sum8 += vec_mijz*v1*w/dst2;
					sum9 += vec_mijz*v2*w/dst2;
					
					sum10 += Prs[j]*w;
				}}
				j = nxt[j];
				if(j==-1) break;
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
		if (p_smooth[i]<0) p_smooth[i] = 0;

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
		if (II[i]<0 || II[i] * 0 != 0) II[i] = 0;
		//II=fabs(S11[i]*S22[i]-S12[i]*S12[i]);
	}}
	
	// Newtonian viscosity
	if (Fluid2_type == 0)
	{
	// Loop only for particles near mesh
#pragma omp parallel for
		for(int im=0;im<nPartNearMesh;im++){
		//for(int i=0;i<nP;i++){
		int i = partNearMesh[im];
		if(Typ[i] == FLD){
			if (PTYPE[i] <= 1)MEU[i] = KNM_VS1 * DNS_FL1;
			if (PTYPE[i] != 1)MEU[i] = KNM_VS2 * DNS_FL2;
		}}

//		if (TURB>0)
//		{
//			NEUt[i] = Cs*DL*Cs*DL*2*sqrt(II[i]);

//			if (NEUt[i] * 0 != 0)  NEUt[i] = 0;
//			if (NEUt[i]>1)     NEUt[i] = 1;
//		}
	}

	// Granular Fluid
	if (Fluid2_type == 1)
	{
		// PROBLEMS TO USE OPENMP HERE. MAYBE THE ACCESS TO BL, WL
		// Loop only for particles near mesh
//#pragma omp parallel for
		for(int im=0;im<nPartNearMesh;im++){
		//for(int i=0;i<nP;i++){
		int i = partNearMesh[im];
		if(Typ[i] == FLD){

//			double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
			double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
			double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];

			if(PTYPE[i] == 1) 
				MEU[i] = KNM_VS1 * DNS_FL1;
			else if(PTYPE[i] == 2){
				phi = (C[i] - 0.25)*PHI / (1 - 0.25);
				phi2 = (C[i] - 0.25)*PHI_2 / (1 - 0.25);
				if (C[i] <= 0.25) { phi = 0.00001; phi2 = 0.00001; } // phi close to zero
				if (PTYPE[i] <= 0) phi = PHI_BED;

				// Normal stress calculation (mechanical pressure)
				p_rheo_new[i] = p_smooth[i];

				kx = int((pos_ix - Xmin) / aa / PCL_DST) + 1;
				ky = int((pos_iy - Ymin) / aa / PCL_DST) + 1;

				// Effective pressure = total pressure (from EOS) - hydrostatic pressure
				//normal_stress=(BL[k]-pos_iy+DL/2)*(DNS_FL2)*9.81;	// normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - pos_iz + PCL_DST / 2)*(DNS_FL2 - DNS_FL1)*9.81 - (vec_ix*vec_ix + vec_iy*vec_iy + vec_iz*vec_iz)*(DNS_FL2 - DNS_FL1) / 2.0;	// normal_stress= Gama.H

				if (p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81<0) p_smooth[i] = (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
				if (TIM <= 1) normal_stress = 1.0*(1 - TIM)*(p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81) + 1.0*(TIM)*normal_stress;

//				normal_stress = p_smooth[i] - (WL[kx][ky] - pos_iz)*DNS_FL1*9.81;
//				normal_stress = p_smooth[i];
				//normal_stress=normal_stress*0.61*1500/DNS_FL2;
				if (normal_stress < 1 || C[i] < 0.5) normal_stress = 1;

				p_rheo_new[i] = normal_stress;

				// Yield stress calculation
				//Inertia[i] = sqrt(II[i])*dg/sqrt(normal_stress/DNS_SDT);		// Free-fall (dry granular material)
				Inertia[i] = sqrt(II[i])*dg/sqrt(normal_stress/(DNS_FL1*Cd));	// Grain inertia (submerged)
				//Inertia[i] = sqrt(II[i])*(KNM_VS1*DNS_FL1)/normal_stress ;	// Viscous regime

				// VF_max VF_min
				VF[i] = VF_max - (VF_max - VF_min)*Inertia[i];
				if (VF[i] < VF_min) VF[i] = VF_min;
				RHO[i] = DNS_SDT * VF[i] + (1-VF[i])*DNS_FL1;
				phi = phi * VF[i] / VF_max;

				double yield_stress = cohes * cos(phi) + normal_stress * sin(phi);

				if (yield_stress < 0) yield_stress = 0;

				double visc_max = (yield_stress*mm*0.5 + MEU0);

				if (II[i]>0)
					MEU_Y[i] = yield_stress * (1 - exp(-mm * sqrt(II[i]))) / 2.0 / sqrt(II[i]);
				else
					MEU_Y[i] = visc_max;

				// H-B rheology

				//meu0 = MEU0;

				// Non-linear Meu(I) rheology
				//meu0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*dg/(I0*sqrt(normal_stress/DNS_FL2)+sqrt(II[i])*dg);			//free fall
				meu0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*dg/(I0*sqrt(normal_stress/(DNS_FL1*Cd))+sqrt(II[i])*dg);		//grain inertia
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*normal_stress*(KNM_VS1*DNS_FL1)/(I0*normal_stress+sqrt(II[i])*(KNM_VS1*DNS_FL1));	//viscous
			   	
			   	// Linear Meu(I) rheology
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*dg*sqrt(normal_stress*DNS_FL2)/I0;		//free fall
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*dg*sqrt(normal_stress*DNS_FL1*Cd)/I0;	//grain inertia
			   	//meu0 = 0.5*(tan(phi2) - tan(phi))*(KNM_VS1*DNS_FL1)/I0;					//viscous

				if (II[i] <= 0 || (meu0 * 0) != 0) meu0 = MEU0;

				visc_max = (yield_stress*mm*0.5 + meu0);

				// Herschel bulkley papanastasiou
				MEU[i] = MEU_Y[i] + MEU0 * pow(4 * II[i], (N - 1) / 2);

				// MEU_Y rheological model
				//MEU[i] = MEU_Y[i] + meu0;
				
				if (II[i] == 0 || MEU[i]>visc_max) MEU[i] = visc_max;
				if (PTYPE[i] <= 0) MEU[i] = MEU[i] * C[i] + DNS_FL1*KNM_VS1*(1 - C[i]);
			}
			
			if(PTYPE[i] >= 2) {
				if (C[i] > 0.5) RHO[i] = DNS_FL2;
				else RHO[i] = C[i] * DNS_FL2 + (1 - C[i]) * DNS_FL1;
			}
		}}

		//---------------------------------- Direct stress calculation method -----------------------------------------
//		if (stress_cal_method == 2)
//		{
//			for (i = 1; i <= NUM; i++)
//			{
//				double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
//				for (l = 2; l <= neighb[i][1]; l++)
//				{
//					j = neighb[i][l];
//					d = DIST(i, j);
//					if (i != j && d <= re)
//					{
//						w = W(d, KTYPE, 2);

//						double meuij = 2 * MEU[i] * MEU[j] / (MEU[i] + MEU[j]);
//						if ((NEUt[i] + NEUt[j])>0) meuij = meuij + 2 * NEUt[i] * RHO[i] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j]);

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

	} // if (Fluid2_type == 1)

	//---------------------------------------------------------------

	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS; //delete[]p_smooth;
	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL;// p_smooth = NULL;
}

void VscTrm_omp(){
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0;i<nP;i++){
	if(Typ[i] == FLD){
		double meu_i = MEU[i];
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double vec_ix = Vel[i*3  ];	double vec_iy = Vel[i*3+1];	double vec_iz = Vel[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2<r2){
				if(j!=i && Typ[j]!=GST){
					double dst = sqrt(dst2);
					double w = WEI(dst, r);

					NEU = 2 * meu_i * MEU[j] / (meu_i + MEU[j]);

//NEU = KNM_VS2 * DNS_FL2;

					if (PTYPE[i] == 1) NEU = NEU/DNS_FL1;
					else NEU = NEU/DNS_FL2;
//					NEU = NEU/RHO[i];
					
					//if ((NEUt[i] + NEUt[j])>0) NEU = NEU + (2 * NEUt[i] * RHO[j] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j])) / RHO[i];

					// Original
//					Acc_x +=(Vel[j*3  ]-vec_ix)*w;
//					Acc_y +=(Vel[j*3+1]-vec_iy)*w;
//					Acc_z +=(Vel[j*3+2]-vec_iz)*w;
					// Modified
					Acc_x +=(Vel[j*3  ]-vec_ix)*w*NEU;
					Acc_y +=(Vel[j*3+1]-vec_iy)*w*NEU;
					Acc_z +=(Vel[j*3+2]-vec_iz)*w*NEU;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		// Original
//		Acc[i*3  ]=Acc_x*A1 + G_X;
//		Acc[i*3+1]=Acc_y*A1 + G_Y;
//		Acc[i*3+2]=Acc_z*A1 + G_Z;
		// Modified
		Acc[i*3  ]=Acc_x*A1_M + G_X;
		Acc[i*3+1]=Acc_y*A1_M + G_Y;
		Acc[i*3+2]=Acc_z*A1_M + G_Z;
	}}
}

void UpPcl1_omp(){
#pragma omp parallel for
	for(int i=0;i<nP;i++){
		if(Typ[i] == FLD){
			Vel[i*3  ] +=Acc[i*3  ]*DT;	Vel[i*3+1] +=Acc[i*3+1]*DT;	Vel[i*3+2] +=Acc[i*3+2]*DT;
			Pos[i*3  ] +=Vel[i*3  ]*DT;		Pos[i*3+1] +=Vel[i*3+1]*DT;		Pos[i*3+2] +=Vel[i*3+2]*DT;
			Acc[i*3]=Acc[i*3+1]=Acc[i*3+2]=0.0;
			F1[i*3]=F1[i*3+1]=F1[i*3+2]=0.0;
			F2[i*3]=F2[i*3+1]=F2[i*3+2]=0.0;
			ChkPcl(i);
		}
	}
}

void ChkCol_omp(){
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0;i<nP;i++){
	if(Typ[i] == FLD){
//		double mi = Dns[Typ[i]];
		double mi;
		if (PTYPE[i] == 1) mi = DNS_FL1;
		else mi = DNS_FL2;

		double pos_ix =  Pos[i*3  ];double pos_iy =  Pos[i*3+1];double pos_iz =  Pos[i*3+2];
		double vec_ix =  Vel[i*3  ];double vec_iy =  Vel[i*3+1];double vec_iz =  Vel[i*3+2];
		double vec_ix2 = Vel[i*3  ];double vec_iy2 = Vel[i*3+1];double vec_iz2 = Vel[i*3+2];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2<rlim2){
				if(j!=i && Typ[j]!=GST){
					double fDT = (vec_ix-Vel[j*3  ])*v0+(vec_iy-Vel[j*3+1])*v1+(vec_iz-Vel[j*3+2])*v2;
					if(fDT > 0.0){
//						double mj = Dns[Typ[j]];
						double mj;
						if (PTYPE[j] == 1) mj = DNS_FL1;
						else mj = DNS_FL2;

						fDT *= COL*mj/(mi+mj)/dst2;
						vec_ix2 -= v0*fDT;		vec_iy2 -= v1*fDT;		vec_iz2 -= v2*fDT;
					}
					/*
					double fDT = (Vel[j*3  ]-vec_ix)*v0+(Vel[j*3+1]-vec_iy)*v1+(Vel[j*3+2]-vec_iz)*v2;
					double mj = Dns[Typ[j]];
					fDT *= COL*mj/(mi+mj)/dst2;
					vec_ix2 += v0*fDT;		vec_iy2 += v1*fDT;		vec_iz2 += v2*fDT;
					*/
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		Acc[i*3  ]=vec_ix2;	Acc[i*3+1]=vec_iy2;	Acc[i*3+2]=vec_iz2;
	}}
#pragma omp parallel for
	for(int i=0;i<nP;i++){
		// CHANGED !!!
		//Pos[i*3  ]+=(Acc[i*3  ]-Vel[i*3  ])*DT; Pos[i*3+1]+=(Acc[i*3+1]-Vel[i*3+1])*DT; Pos[i*3+2]+=(Acc[i*3+2]-Vel[i*3+2])*DT;
		Vel[i*3  ]=Acc[i*3  ];	Vel[i*3+1]=Acc[i*3+1];	Vel[i*3+2]=Acc[i*3+2];

		//Velk[i*3  ]=Vel[i*3  ];	Velk[i*3+1]=Vel[i*3+1];	Velk[i*3+2]=Vel[i*3+2];
		//Pos[i*3  ]=Posk[i*3  ]+Vel[i*3  ]*DT; Pos[i*3+1]=Posk[i*3+1]+Vel[i*3+1]*DT; Pos[i*3+2]=Posk[i*3+2]+Vel[i*3+2]*DT;
		Acc[i*3  ]=0.0;	Acc[i*3+1]=0.0;	Acc[i*3+2]=0.0;	
	}
}

void MkPrs_omp(){
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0;i<nP;i++){
	if(Typ[i] != GST){
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double ni = 0.0;
		// Add PND due Wall polygon
		ni += niw[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2<r2){
				if(j!=i && Typ[j]!=GST){
					double dst = sqrt(dst2);
					double w = WEI(dst, r);
					ni += w;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		pndi[i] = ni;
//		double mi = Dns[Typ[i]];
		double mi;
		if (PTYPE[i] == 1) mi = DNS_FL1;
		else mi = DNS_FL2;

		double pressure = 0.0;
		Bc[i] = SRF;
		if (ni > PND_TRS*n0){
			Bc[i] = INR;
			pressure = (ni - n0) * A2 * mi;
		}
		if (pressure < 0.0)
			pressure = 0.0;
		Prs[i] = pressure;
	}}
}

void MkPrsWc_omp(){
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0;i<nP;i++){
	if(Typ[i] != GST){
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double ni = 0.0;
		// Add PND due Wall polygon
		ni += niw[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2<r2){
				if(j!=i && Typ[j]!=GST){
					double dst = sqrt(dst2);
					double w = WEI(dst, r);
					ni += w;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		pndi[i] = ni;
//		double mi = Dns[Typ[i]];
		double mi;
		if (PTYPE[i] == 1) mi = DNS_FL1;
		else mi = DNS_FL2;

		double pressure = 0.0;
		Bc[i] = SRF;
		if (ni > PND_TRS*n0){
			Bc[i] = INR;
			pressure = (mi*A4/GAM)*(pow(ni/n0,GAM)-1);
			//pressure = (ni - n0) * A2 * mi;
		}
		if (pressure < 0.0)
			pressure = 0.0;
		Prs[i] = pressure;
	}}
}

void PrdPrsGrdTrm_omp(){
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0;i<nP;i++){
	if(Typ[i] == FLD){
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pre_min = Prs[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2<r2){
				if(j!=i && Typ[j]!=GST){
					if(pre_min > Prs[j])pre_min = Prs[j];
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2<r2){
				if(j!=i && Typ[j]!=GST){
					double dst = sqrt(dst2);
					double w = WEI_GRAD(dst, r);
					if(GRD_TYP == 0)
						w *= (Prs[j] - pre_min)/dst2;
					else if(GRD_TYP == 1)
						w *= (Prs[j] + Prs[i])/dst2;
					else if(GRD_TYP == 2)
						w *= (Prs[j] + Prs[i] - 2*pre_min)/dst2;
					Acc_x += v0*w;	Acc_y += v1*w;	Acc_z += v2*w;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
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

void PrsGrdTrm_omp(){
#pragma omp parallel for schedule(dynamic,64)
	for(int i=0;i<nP;i++){
	if(Typ[i] == FLD){
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pre_min = Prs[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2<r2){
				if(j!=i && Typ[j]!=GST){
					if(pre_min > Prs[j])pre_min = Prs[j];
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				double v0 = Pos[j*3  ] - pos_ix;
				double v1 = Pos[j*3+1] - pos_iy;
				double v2 = Pos[j*3+2] - pos_iz;
				double dst2 = v0*v0+v1*v1+v2*v2;
				if(dst2<r2){
				if(j!=i && Typ[j]!=GST){
					double dst = sqrt(dst2);
					double w = WEI_GRAD(dst, r);
					if(GRD_TYP == 0)
						w *= (Prs[j] - pre_min)/dst2;
					else if(GRD_TYP == 1)
						w *= (Prs[j] + Prs[i])/dst2;
					else if(GRD_TYP == 2)
						w *= (Prs[j] + Prs[i] - 2*pre_min)/dst2;
					
					Acc_x += v0*w;	Acc_y += v1*w;	Acc_z += v2*w;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		// Original
//		Acc[i*3  ]=RLX_PRS*Acc_x*invDns[FLD]*A3;
//		Acc[i*3+1]=RLX_PRS*Acc_y*invDns[FLD]*A3;
//		Acc[i*3+2]=RLX_PRS*Acc_z*invDns[FLD]*A3;
		// Modified
		Acc[i*3  ]=RLX_PRS*Acc_x*A3/RHO[i];
		Acc[i*3+1]=RLX_PRS*Acc_y*A3/RHO[i];
		Acc[i*3+2]=RLX_PRS*Acc_z*A3/RHO[i];
	}}
}

void WallPrdPrsGrdTrm_omp(){
	int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	for(int im=0;im<nPartNearMesh;im++){
	//for(int i=0;i<nP;i++){
	int i = partNearMesh[im];
	if(Typ[i] == FLD){
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double pre_min = Prs[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;
				// If inside neighboor of i
				if(dstij2<r2){
				if(j!=i && Typ[j]!=GST){
					if(pre_min > Prs[j])pre_min = Prs[j];
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0 = Pos[j*3  ] - pos_mix;
				double v1 = Pos[j*3+1] - pos_miy;
				double v2 = Pos[j*3+2] - pos_miz;

				double dst2 = v0*v0+v1*v1+v2*v2;
				// If inside neighboor of i and im (intersection)
				if(dstij2<r2 && dst2<r2){
				if(j!=i && Typ[j]!=GST){

					Nw[i]=1; // Only to show particles near polygon

					double dst = sqrt(dst2);
					double w = WEI_GRAD(dst, r);
					if(GRD_TYP == 0)
						w *= (Prs[j] - pre_min)/dst2;
					else if(GRD_TYP == 1)
						w *= (Prs[j] + Prs[i])/dst2;
					else if(GRD_TYP == 2)
						w *= (Prs[j] + Prs[i] - 2*pre_min)/dst2;
					Acc_x += v0*w;	Acc_y += v1*w;	Acc_z += v2*w;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		// Add "i" contribution ("i" is a neighboor of "mirror i")
	  	double v0 = pos_ix - pos_mix;
		double v1 = pos_iy - pos_miy;
		double v2 = pos_iz - pos_miz;
		double dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2<r2){

			Nw[i]=1; // Only to show particles near polygon

			double dst = sqrt(dst2);
			double w = WEI_GRAD(dst, r);
			if(GRD_TYP == 0)
				w *= (Prs[i] - pre_min)/dst2;
			else if(GRD_TYP == 1)
				w *= (Prs[i] + Prs[i])/dst2;
			else if(GRD_TYP == 2)
				w *= (Prs[i] + Prs[i] - 2*pre_min)/dst2;
			Acc_x += v0*w;	Acc_y += v1*w;	Acc_z += v2*w;
	  	}
		// Wall gradient Mitsume`s model
	    double Rref_i[9], normaliw[3], normaliwSqrt;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);

	    if (normaliwSqrt > 0.00000001) {
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

		// Repulsive force
	  	double rpsForce[3];

	  	rpsForce[0]=rpsForce[1]=rpsForce[2] = 0.0;
	  	if (normaliwSqrt < 0.5*PCL_DST && normaliwSqrt != 0){
    		rpsForce[0] = - ARF*(0.5*PCL_DST/(normaliwSqrt) - 1)*normaliw[0];
    		rpsForce[1] = - ARF*(0.5*PCL_DST/(normaliwSqrt) - 1)*normaliw[1];
    		rpsForce[2] = - ARF*(0.5*PCL_DST/(normaliwSqrt) - 1)*normaliw[2];

    		//if (i==6){
    			//printf("x:%lf y:%lf z:%lf x:%lf y:%lf z:%lf\n", Pos[i*3],Pos[i*3+1],Pos[i*3+2],mirrorPos[i*3],mirrorPos[i*3+1],mirrorPos[i*3+2]);
    			//printf("nx:%lf ny: %lf nz: %lf nN:%lf\n", normaliw[0],normaliw[1],normaliw[2],normaliwSqrt);
    			//printf("Fx:%lf Fy:%lf Fz:%lf\n", rpsForce[0],rpsForce[1],rpsForce[2]);
    		//}
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

void WallPrsGrdTrm_omp(){
	int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	for(int im=0;im<nPartNearMesh;im++){
	//for(int i=0;i<nP;i++){
	int i = partNearMesh[im];
	if(Typ[i] == FLD){
		double Acc_x = 0.0;			double Acc_y = 0.0;			double Acc_z = 0.0;
		double pos_ix = Pos[i*3  ];	double pos_iy = Pos[i*3+1];	double pos_iz = Pos[i*3+2];
		double pos_mix = mirrorPos[i*3  ];	double pos_miy = mirrorPos[i*3+1];	double pos_miz = mirrorPos[i*3+2];
		double pre_min = Prs[i];
		int ix = (int)((pos_ix - MIN_X)*DBinv) +1;
		int iy = (int)((pos_iy - MIN_Y)*DBinv) +1;
		int iz = (int)((pos_iz - MIN_Z)*DBinv) +1;
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;
				// If inside neighboor of i
				if(dstij2<r2){
				if(j!=i && Typ[j]!=GST){
					if(pre_min > Prs[j])pre_min = Prs[j];
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0 = Pos[j*3  ] - pos_mix;
				double v1 = Pos[j*3+1] - pos_miy;
				double v2 = Pos[j*3+2] - pos_miz;

				double dst2 = v0*v0+v1*v1+v2*v2;

				// If inside neighboor of i and im (intersection)
				if(dstij2<r2 && dst2<r2){
				if(j!=i && Typ[j]!=GST){

					Nw[i]=1; // Only to show particles near polygon

					double dst = sqrt(dst2);
					double w = WEI_GRAD(dst, r);
					if(GRD_TYP == 0)
						w *= (Prs[j] - pre_min)/dst2;
					else if(GRD_TYP == 1)
						w *= (Prs[j] + Prs[i])/dst2;
					else if(GRD_TYP == 2)
						w *= (Prs[j] + Prs[i] - 2*pre_min)/dst2;
					Acc_x += v0*w;	Acc_y += v1*w;	Acc_z += v2*w;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}
		// Add "i" contribution ("i" is a neighboor of "mirror i")
	  	double v0 = pos_ix - pos_mix;
		double v1 = pos_iy - pos_miy;
		double v2 = pos_iz - pos_miz;
		double dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2<r2){

			Nw[i]=1; // Only to show particles near polygon

			double dst = sqrt(dst2);
			double w = WEI_GRAD(dst, r);
			if(GRD_TYP == 0)
				w *= (Prs[i] - pre_min)/dst2;
			else if(GRD_TYP == 1)
				w *= (Prs[i] + Prs[i])/dst2;
			else if(GRD_TYP == 2)
				w *= (Prs[i] + Prs[i] - 2*pre_min)/dst2;
			Acc_x += v0*w;	Acc_y += v1*w;	Acc_z += v2*w;
	  	}
		// Wall gradient Mitsume`s model
	    double Rref_i[9], normaliw[3], normaliwSqrt;
	    // Normal fluid-wall particle = 0.5*(Normal fluid-mirror particle)
	    normaliw[0] = 0.5*(pos_ix - pos_mix); normaliw[1] = 0.5*(pos_iy - pos_miy); normaliw[2] = 0.5*(pos_iz - pos_miz);
	    normaliwSqrt = sqrt(normaliw[0]*normaliw[0] + normaliw[1]*normaliw[1] + normaliw[2]*normaliw[2]);

	    if (normaliwSqrt > 0.00000001) {
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

		// Repulsive force
	  	double rpsForce[3];

	  	rpsForce[0]=rpsForce[1]=rpsForce[2] = 0.0;
	  	if (normaliwSqrt < 0.5*PCL_DST && normaliwSqrt != 0){
    		rpsForce[0] = - ARF*(0.5*PCL_DST/(normaliwSqrt) - 1)*normaliw[0];
    		rpsForce[1] = - ARF*(0.5*PCL_DST/(normaliwSqrt) - 1)*normaliw[1];
    		rpsForce[2] = - ARF*(0.5*PCL_DST/(normaliwSqrt) - 1)*normaliw[2];
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

void WallSlipVscTrm_omp(){
	int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", nPartNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	for(int im=0;im<nPartNearMesh;im++){
	//for(int i=0;i<nP;i++){
	int i = partNearMesh[im];
	if(Typ[i] == FLD){
		
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

	    if (normalMod2 > 0.00000001) {
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
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0 = Pos[j*3  ] - pos_mix;
				double v1 = Pos[j*3+1] - pos_miy;
				double v2 = Pos[j*3+2] - pos_miz;

				double dst2 = v0*v0+v1*v1+v2*v2;
				// If inside neighboor of i and im (intersection)
				if(dstij2<r2 && dst2<r2){
				if(j!=i && Typ[j]!=GST){
					double dst = sqrt(dst2);
					double w = WEI(dst, r);

					NEU = 2 * meu_i * MEU[j] / (meu_i + MEU[j]);

//NEU = KNM_VS2 * DNS_FL2;

					if (PTYPE[i] == 1) NEU = NEU/DNS_FL1;
					else NEU = NEU/DNS_FL2;

					//if ((NEUt[i] + NEUt[j])>0) NEU = NEU + (2 * NEUt[i] * RHO[j] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j])) / RHO[i];

					// Original
//					Acc_x +=(Vel[j*3  ]-vec_mix)*w;
//					Acc_y +=(Vel[j*3+1]-vec_miy)*w;
//					Acc_z +=(Vel[j*3+2]-vec_miz)*w;
					// Modified
					Acc_x +=(Vel[j*3  ]-vec_mix)*w*NEU;
					Acc_y +=(Vel[j*3+1]-vec_miy)*w*NEU;
					Acc_z +=(Vel[j*3+2]-vec_miz)*w*NEU;

					//Acc_x +=(Velk[j*3  ]-vec_mix)*w;
					//Acc_y +=(Velk[j*3+1]-vec_miy)*w;
					//Acc_z +=(Velk[j*3+2]-vec_miz)*w;
				}}
				j = nxt[j];
				if(j==-1) break;
			}
		}}}

		// Add "i" contribution ("i" is a neighboor of "mirror i")
		double v0 = pos_ix - pos_mix;
		double v1 = pos_iy - pos_miy;
		double v2 = pos_iz - pos_miz;
		double dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2<r2){
			double dst = sqrt(dst2);
			double w = WEI(dst, r);

			NEU = 2 * meu_i * meu_i / (meu_i + meu_i);

//NEU = KNM_VS2 * DNS_FL2;

			if (PTYPE[i] == 1) NEU = NEU/DNS_FL1;
			else NEU = NEU/DNS_FL2;

			// Original
//			Acc_x +=(vec_ix-vec_mix)*w;
//			Acc_y +=(vec_iy-vec_miy)*w;
//			Acc_z +=(vec_iz-vec_miz)*w;

			// Modified
			Acc_x +=(vec_ix-vec_mix)*w*NEU;
			Acc_y +=(vec_iy-vec_miy)*w*NEU;
			Acc_z +=(vec_iz-vec_miz)*w*NEU;
	  	}

	  	// Wall laplacian Mitsume`s model
      	// Correction of velocity
      	// Original
//      Acc[i*3  ] += (Rref_i[0]*Acc_x + Rref_i[1]*Acc_y + Rref_i[2]*Acc_z)*A1;
//		Acc[i*3+1] += (Rref_i[3]*Acc_x + Rref_i[4]*Acc_y + Rref_i[5]*Acc_z)*A1;
//		Acc[i*3+2] += (Rref_i[6]*Acc_x + Rref_i[7]*Acc_y + Rref_i[8]*Acc_z)*A1;

		// Modified
		Acc[i*3  ] += (Rref_i[0]*Acc_x + Rref_i[1]*Acc_y + Rref_i[2]*Acc_z)*A1_M;
		Acc[i*3+1] += (Rref_i[3]*Acc_x + Rref_i[4]*Acc_y + Rref_i[5]*Acc_z)*A1_M;
		Acc[i*3+2] += (Rref_i[6]*Acc_x + Rref_i[7]*Acc_y + Rref_i[8]*Acc_z)*A1_M;
	}}
}


void WallNoSlipVscTrm_omp(){
	int nPartNearMesh = partNearMesh.size();
	//printf(" Mesh %d \n", partNearMesh);
	// Loop only for particles near mesh
#pragma omp parallel for schedule(dynamic,64)
	for(int im=0;im<nPartNearMesh;im++){
	//for(int i=0;i<nP;i++){
	int i = partNearMesh[im];
	if(Typ[i] == FLD){

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

	    if (normalMod2 > 0.00000001) {
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
		for(int jz=iz-1;jz<=iz+1;jz++){
		for(int jy=iy-1;jy<=iy+1;jy++){
		for(int jx=ix-1;jx<=ix+1;jx++){
			int jb = jz*nBxy + jy*nBx + jx;
			int j = bfst[jb];
			if(j == -1) continue;
			for(;;){
				// Particle distance r_ij = Xj - Xi_temporary_position
				double v0ij = Pos[j*3  ] - pos_ix;
				double v1ij = Pos[j*3+1] - pos_iy;
				double v2ij = Pos[j*3+2] - pos_iz;

				double dstij2 = v0ij*v0ij+v1ij*v1ij+v2ij*v2ij;

				// Mirror particle distance r_imj = Xj - Xim_temporary_position
				double v0 = Pos[j*3  ] - pos_mix;
				double v1 = Pos[j*3+1] - pos_miy;
				double v2 = Pos[j*3+2] - pos_miz;

				double dst2 = v0*v0+v1*v1+v2*v2;
				// If inside neighboor of i and im (intersection)
				if(dstij2<r2 && dst2<r2){
				if(j!=i && Typ[j]!=GST){
					double dst = sqrt(dst2);
					double w = WEI(dst, r);

					NEU = 2 * meu_i * MEU[j] / (meu_i + MEU[j]);

//NEU = KNM_VS2 * DNS_FL2;


					if (PTYPE[i] == 1) NEU = NEU/DNS_FL1;
					else NEU = NEU/DNS_FL2;
					
					//if ((NEUt[i] + NEUt[j])>0) NEU = NEU + (2 * NEUt[i] * RHO[j] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j])) / RHO[i];

					// Original
//					Acc_x +=(Vel[j*3  ]-vec_mix)*w;
//					Acc_y +=(Vel[j*3+1]-vec_miy)*w;
//					Acc_z +=(Vel[j*3+2]-vec_miz)*w;
					// Modified
					Acc_x +=(Vel[j*3  ]-vec_mix)*w*NEU;
					Acc_y +=(Vel[j*3+1]-vec_miy)*w*NEU;
					Acc_z +=(Vel[j*3+2]-vec_miz)*w*NEU;
					
					//Acc_x +=(Velk[j*3  ]-vec_mix)*w;
					//Acc_y +=(Velk[j*3+1]-vec_miy)*w;
					//Acc_z +=(Velk[j*3+2]-vec_miz)*w;

					//if(i==2817){
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
				if(j==-1) break;
			}
		}}}

		// Add "i" contribution ("i" is a neighboor of "mirror i")
	  	double v0 = pos_ix - pos_mix;
		double v1 = pos_iy - pos_miy;
		double v2 = pos_iz - pos_miz;
		double dst2 = v0*v0+v1*v1+v2*v2;
		if(dst2<r2){
			double dst = sqrt(dst2);
			double w = WEI(dst, r);

			NEU = 2 * meu_i * meu_i / (meu_i + meu_i);

//NEU = KNM_VS2 * DNS_FL2;

			if (PTYPE[i] == 1) NEU = NEU/DNS_FL1;
			else NEU = NEU/DNS_FL2;

			// Original
//			Acc_x +=(vec_ix-vec_mix)*w;
//			Acc_y +=(vec_iy-vec_miy)*w;
//			Acc_z +=(vec_iz-vec_miz)*w;

			// Modified
			Acc_x +=(vec_ix-vec_mix)*w*NEU;
			Acc_y +=(vec_iy-vec_miy)*w*NEU;
			Acc_z +=(vec_iz-vec_miz)*w*NEU;

			//Acc_x += (Rinv_i[0]*(vec_ix-vec_mix)+ Rinv_i[1]*(vec_iy-vec_miy) + Rinv_i[2]*(vec_iz-vec_miz))*w;
			//Acc_y += (Rinv_i[3]*(vec_ix-vec_mix)+ Rinv_i[4]*(vec_iy-vec_miy) + Rinv_i[5]*(vec_iz-vec_miz))*w;
			//Acc_z += (Rinv_i[6]*(vec_ix-vec_mix)+ Rinv_i[7]*(vec_iy-vec_miy) + Rinv_i[8]*(vec_iz-vec_miz))*w;
	  	}

	  	//if (i==6){
	  	//	printf("Vx:%lf Vy:%lf Vz:%lf Vx:%lf Vy:%lf Vz:%lf\n",vec_ix,vec_mix,vec_iy,vec_miy,vec_iz,vec_miz);
		//	printf("Accx:%lf Bccy:%lf Bccz:%lf \n", Acc[i*3],Acc[i*3+1],Acc[i*3+2]);
			//printf("Bccx:%lf Bccy:%lf Bccz:%lf \n", Acc[i*3],Acc[i*3+1],Acc[i*3+2]);
		//}
		// Wall laplacian Mitsume`s model
      	// Correction of velocity
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

		//if(i==2817){
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
		//if(i==6){
			//printf("Accx:%lf Accy:%lf Accz:%lf \n", Acc[i*3],Acc[i*3+1],Acc[i*3+2]);
		//	printf("Time:%e\n", TIM);
		//	printf("Xi:%e %e %e Xm:%e %e %e\n", pos_ix,pos_iy,pos_iz,pos_mix,pos_miy,pos_miz);
		//	printf("Vi:%e %e %e Vm:%e %e %e\n", vec_ix,vec_iy,vec_iz,vec_mix,vec_miy,vec_miz);
		//	printf("Acc:%e %e %e\n", AA[0],AA[1],AA[2]);
		//}
	}}
}

void UpPcl2_omp(void){
	vMax = 0.0;						// Maximum velocity
#pragma omp parallel for
	for(int i=0;i<nP;i++){
		if(Typ[i] == FLD){
			Vel[i*3  ]+=Acc[i*3  ]*DT;	Vel[i*3+1]+=Acc[i*3+1]*DT;	Vel[i*3+2]+=Acc[i*3+2]*DT;
			Pos[i*3  ]+=Acc[i*3  ]*DT*DT;	Pos[i*3+1]+=Acc[i*3+1]*DT*DT;	Pos[i*3+2]+=Acc[i*3+2]*DT*DT;
			Acc[i*3]=Acc[i*3+1]=Acc[i*3+2]=0.0;

			//Pos[i*3  ]=Posk[i*3 ]+Vel[i*3  ]*DT;	Pos[i*3+1]=Posk[i*3+1]+Vel[i*3+1]*DT;	Pos[i*3+2]=Posk[i*3+2]+Vel[i*3+2]*DT;
			//Posk[i*3  ]=Pos[i*3  ];	Posk[i*3+1]=Pos[i*3+1];	Posk[i*3+2]=Pos[i*3+2];
			//Velk[i*3  ]=Vel[i*3  ];	Velk[i*3+1]=Vel[i*3+1];	Velk[i*3+2]=Vel[i*3+2];

			//F1[i*3  ]=Acc[i*3  ];	F1[i*3+1]=Acc[i*3+1];	F1[i*3+2]=Acc[i*3+2];
			//F2[i*3  ]=Acv[i*3  ];	F2[i*3+1]=Acv[i*3+1];	F2[i*3+2]=Acv[i*3+2];
			//Acc[i*3]=Acc[i*3+1]=Acc[i*3+2]=0.0;
			//Acv[i*3]=Acv[i*3+1]=Acv[i*3+2]=0.0;

			double vMod2 = Vel[i*3  ]*Vel[i*3  ] + Vel[i*3+1]*Vel[i*3+1] + Vel[i*3+2]*Vel[i*3+2];
			if (vMod2 > vMax*vMax)
				vMax = sqrt(vMod2);
			ChkPcl(i);
		}
	}
	Crt = DT*vMax/PCL_DST;
}

void ClcEMPS(mesh mesh){
	while(1){
		if(iLP%100==0){
			int p_num=0;
			for(int i=0;i<nP;i++){if(Typ[i] != GST)p_num++;}
			timer_end = get_dtime();
			printf("%5d th TIM: %lf / p_num: %d / Vmax: %lf / Courant: %lf", iLP,TIM,p_num, vMax, Crt);
			printf(" / RunTime: %13.6lf sec\n",timer_end -timer_sta);
		}
		if(iLP%OPT_FQC == 0 ){
			//WrtDat();
			WrtVtu();
			if(TIM >= FIN_TIM ){	break;}
		}
		// Small subdomains (buckets)
		MkBkt();
		
//		if(Fluid2_type==1){
//			VolFract_omp(); // Calculation of the volume of fraction if phase II in the mixture
//			VscIntVal_omp();// Calculation of viscosity
//		}
		// Calculation of acceleration due laplacian of viscosity and gravity
		VscTrm_omp();
		// Add acceleration due pressure gradient (Prediction)
		if(RLX_PRS!=1.0){
			PrdPrsGrdTrm_omp();
			// Pressure gradient on polygon wall
			WallPrdPrsGrdTrm_omp();
		}
		// Update velocity and positions
		UpPcl1_omp();
		// Check collision
		ChkCol_omp();
		// Fluid particles: Calculation of PND due wall and number of neighboors
		// Positions of wall and mirror particles
		mesh.closestPointPNDBoundaryAABB(DIM, r2, nP, Typ, FLD, Pos, wallPos, mirrorPos, niw, numNeigh, partNearMesh);

		// Pressure calculation
		if(MPS_TYP==0)
			MkPrs_omp();
		else if(MPS_TYP==1)
			MkPrsWc_omp();
		// Calculation of acceleration due pressure gradient
		PrsGrdTrm_omp();
		// Add acceleration due pressure gradient on polygon wall
		WallPrsGrdTrm_omp();
		// Add acceleration due laplacian of viscosity on wall
		if(SLP){
			if(Fluid2_type==1){
				VolFract_omp(); // Calculation of the volume of fraction if phase II in the mixture

				VscIntVal_omp();// Calculation of viscosity due fluid particles

				WallSlipVscIntVal_omp();// Calculation of viscosity due polygon
			}
			WallSlipVscTrm_omp(); // Slip
		}
		else{
			if(Fluid2_type==1){
				VolFract_omp(); // Calculation of the volume of fraction if phase II in the mixture

				VscIntVal_omp();// Calculation of viscosity due fluid particles

				WallNoSlipVscIntVal_omp();// Calculation of viscosity due polygon
			}
			WallNoSlipVscTrm_omp(); // No-slip
		}
		// Update velocity and positions
		UpPcl2_omp();
		// Pressure calculation
		if(MPS_TYP==0)
			MkPrs_omp();
		else if(MPS_TYP==1)
			MkPrsWc_omp();

		for(int i=0;i<nP;i++){pav[i] += Prs[i];}
		iLP++;
		TIM += DT;
	}
}

int main( int argc, char** argv) {
	printf("start emps_omp.\n");
    // Creating a directory 
    if (mkdir(OUT_FOLDER, 0777) == -1){
    	printf("Unable to create directory.\n");
    }
    else{
    	printf("Directory created.\n");
    }
	// Read data
	RdDat();
	// Read mesh
	//RdMes();
	fp = fopen(IN_FILE, "r");
	fscanf(fp,"%d",&nP);
	mesh mesh(IN_FILE, nP);
	mesh.readMeshFile(IN_MESH);
	// Allocation
	AlcBkt();
	// Setting parameters
	SetPara();
	// PVD file
	WrtPvd();

	timer_sta = get_dtime();

	ClcEMPS(mesh);

	free(Acc);	free(Pos);	free(Vel);
	free(Prs);	free(pav);	free(Typ);
	free(bfst);	free(blst);	free(nxt);

	free(wallPos);	free(mirrorPos);	free(niw);
	free(pndi);	free(numNeigh); free(Bc);
	free(F1);	free(F2);	free(Nw);
//	free(Posk); free(Velk); free(Acv);

	printf("end emps_omp.\n");
	return 0;
}
