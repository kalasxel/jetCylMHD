#include <iostream> 
#include <cmath>
#include <iomanip>
#include <omp.h>

#include "Param.h"
#include "CFullVec.h"
#include "CHlld.h"
#include "outputPark.h"
#include "setInital.h"

using namespace std;

FullVec TurnIn(FullVec &AAA);
FullVec TurnOut(FullVec &AAA);
void step(FullVec **U);


int main(int argc, char **argv)
{{
	cout << "Space: r0=" << RI << ", rf=" << RF << "; z0=" << ZI << ", zf=" << ZF << ';' << endl;
	cout << "Spatial steps: dr=" << DR << ", dz=" << DZ << ';' << endl;
	cout << "Grid: Nr=" << NR << ", Nz=" << NZ << ';' << endl; 
	cout << "Time: Tfin=" << TFIN << ", dT=" << DT << "; Steps' number: " << ST << ';' << endl;
	cout << "Mn=" << fixed << scientific << Mnuclon << ';' << endl;
	cout << "GO!" << endl;

	omp_set_dynamic(0);      // запретить библиотеке openmp менять число потоков во время исполнения
	omp_set_num_threads(4); // установить число потоков

	double gR[NR+2], gZ[NZ+2];
	setGridR(gR); setGridZ(gZ);
	outputGrid(gR,gZ);

	FullVec **U;
	U = new FullVec* [NR+2];
	for(int i=0; i<NR+2; i++) U[i] = new FullVec [NZ+2];


	setInital(U); 
	setBound(U);

	int st(0), nbF(0); double t(0);
	cout << "time: " << t << "; step: " << st << ';' << endl;

	outputLayer(U,t,nbF);


	while(t<TFIN)
	{
		t += DT; st++;
		if(st%CVAR==0) cout << "time: " << t << "; step: " << st << ';' << endl;

		step(U);

		if(st%10==0)
		{
			nbF++;
			outputLayer(U,t,nbF);
		}
	}


	for(int i=0; i<NR+2; i++) delete [] U[i];
	delete [] U;

	cout << "Calculations has done!" << endl;
	return 0;
}}



FullVec TurnIn(FullVec &AAA)
{{
	FullVec TMP;

	TMP.rho = AAA.rho;

	TMP.u = AAA.v;
	TMP.v = AAA.w;
	TMP.w = AAA.u;

	TMP.B = AAA.H;
	TMP.H = AAA.D;
	TMP.D = AAA.B;

	TMP.e = AAA.e;

	return TMP;
}}

FullVec TurnOut(FullVec &AAA)
{{
	FullVec TMP;

	TMP.rho = AAA.rho;

	TMP.u = AAA.w;
	TMP.v = AAA.u;
	TMP.w = AAA.v;

	TMP.B = AAA.D;
	TMP.H = AAA.B;
	TMP.D = AAA.H;

	TMP.e = AAA.e;

	return TMP;
}}


//#define CARTESIAN  // FOR CARTESIAN COORDINATES
#define VISCOSITY

void step(FullVec **U)
{{
	FullVec **Fr = new FullVec * [NR+1];
	for(int n=0; n<NR+1; n++) Fr[n] = new FullVec [NZ+1];
	FullVec **Fz = new FullVec * [NR+1];
	for(int n=0; n<NR+1; n++) Fz[n] = new FullVec [NZ+1];


			// Stephan-Bolzman radiation
			double Lamda[NR+2][NZ+1];
			for(int n=1; n<=NR; n++)
			{
				for(int m=1; m<=NZ; m++)
				{
					Lamda[n][m] = SteBolz*(U[n][m].getT()*U[n][m].getT())*2*(1/DZ + 1/DR);
					//cout << "SB " << Lamda[n][m] << endl;
					//U[n][m].e += Lamda*DT;
				}
			}


	int n(0), m(0);
#pragma omp parallel for shared(U,Fr,Fz) private(n,m)
	for(m=1; m<=NZ; m++)
	{
		for(n=1; n<=NR; n++)
		{
			FullVec L, R, TMP;
			Hlld flow;

			L=U[n-1][m]; R=U[n][m];
			Fr[n][m] = flow.getFlowHLLD(L,R);

			L=U[n][m-1]; R=U[n][m];
			L=TurnIn(L); R=TurnIn(R);
			L=TurnIn(L); R=TurnIn(R);
			TMP = flow.getFlowHLLD(L,R);
			TMP = TurnOut(TMP);
			Fz[n][m] = TurnOut(TMP);
		}
	}

	#ifdef VISCOSITY
		for(int n=2; n<NR+1; n++)
		{
			for(int m=1; m<NZ+1; m++)
			{
				Fr[n][m].w -= VISCISITY1*( U[n][m].w/U[n][m].rho - U[n][m-1].w/U[n][m-1].rho )/DR;
				Fz[n][m].w -= 2*VISCISITY1*( U[n][m].w/U[n][m].rho - U[n][m-1].w/U[n][m-1].rho )/DZ;
			}
		}

	#endif



	double gR[NR+2]; setGridR(gR);

	for(int n=1; n<NR; n++)
	{
		for(int m=1; m<NZ; m++)
		{
			FullVec RHS;


		#ifdef CARTESIAN
			RHS = Fr[n+1][m] - Fr[n][m];
			RHS = RHS*(DT/DR);
			U[n][m] = U[n][m] - RHS;
		#else // CYLINDRICAL COORDINATES 
			if(true)
			{
				FullVec RADl, RADr;
				RADl = Fr[n][m]*gR[n]; RADr = Fr[n+1][m]*gR[n+1];

				RHS = ( RADr-RADl )*( 2/(gR[n+1]+gR[n]) )*(DT/DR);
				U[n][m] = U[n][m] - RHS;

				FullVec GeomSource;
				double Ptotal = (GAMMA-1)*( U[n][m].e -0.5*U[n][m].Vsq()/U[n][m].rho -0.5*U[n][m].Bsq() ) +0.5*U[n][m].Bsq();
				double Vx, Vy, By; // they will different 
				Vx = U[n][m].v*U[n][m].v/U[n][m].rho + Ptotal - U[n][m].H*U[n][m].H;
				Vy = -U[n][m].u*U[n][m].v/U[n][m].rho + U[n][m].B*U[n][m].H;
				By = U[n][m].u*U[n][m].H/U[n][m].rho - U[n][m].B*U[n][m].v/U[n][m].rho;
				GeomSource.setVec(0, Vx,Vy,0, 0,By,0, 0);
				GeomSource = GeomSource*( 2/(gR[n+1]+gR[n]) )*DT;
				U[n][m] = U[n][m] + GeomSource;
			}
		#endif


			RHS = ( Fz[n][m+1]-Fz[n][m] )*(DT/DZ);
			U[n][m] = U[n][m] - RHS; 

			U[n][m].e += Lamda[n][m]*DT;

		}
	}	


	for(int n=0; n<NR+1; n++) delete [] Fr[n];
	delete [] Fr;
	for(int n=0; n<NR+1; n++) delete [] Fz[n];
	delete [] Fz;

}}
