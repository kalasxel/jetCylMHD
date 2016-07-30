#include <iostream> 
#include <cmath>
#include <iomanip>

#include "Param.h"
#include "CFullVec.h"
#include "CHlld.h"
#include "outputPark.h"
#include "setInital.h"

using namespace std;

FullVec TurnIn(FullVec &AAA);
FullVec TurnOut(FullVec &AAA);
void step(FullVec **U);

double **divCalc(FullVec **U);
void DIV(FullVec **U, int st);


int main(int argc, char **argv)
{{
	cout << "Space: r0=" << RI << ", rf=" << RF << "; z0=" << ZI << ", zf=" << ZF << ';' << endl;
	cout << "Spatial steps: dr=" << DR << ", dz=" << DZ << ';' << endl;
	cout << "Grid: Nr=" << NR << ", Nz=" << NZ << ';' << endl; 
	cout << "Time: Tfin=" << TFIN << ", dT=" << DT << "; Steps' number: " << ST << ';' << endl;
	cout << "GO!" << endl;

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
	//DIV(U,st);


	while(t<TFIN)
	{
		t += DT; st++;
		if(st%CVAR==0) cout << "time: " << t << "; step: " << st << ';' << endl;

		step(U);

		if(st%10==0)
		{
			nbF++;
			outputLayer(U,t,nbF);
			//DIV(U,nbF);
		}
	}

	//step(U);

/*
	FullVec J,I, RRR, tJ, tI;
	Hlld flow;

	J.setVec(0.1,-0.1,2,4,5,2,3,0.4); I.setVec(2,6,7,8,1,2,3,0.4);
	tJ=J; tI=I;
	//J.slimU(); printVec(J);
	printVec(J,I);
	RRR = flow.getFlowHLLD(J,I);
	printVec(RRR,tJ);
*/




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

void step(FullVec **U)
{{
	FullVec **Fr = new FullVec * [NR+1];
	for(int n=0; n<NR+1; n++) Fr[n] = new FullVec [NZ+1];
	FullVec **Fz = new FullVec * [NR+1];
	for(int n=0; n<NR+1; n++) Fz[n] = new FullVec [NZ+1];


	for(int m=1; m<=NZ; m++)
	{
		for(int n=1; n<=NR; n++)
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
				//cout << n << ' ' << RADl.u << ' ' << RADr.u << endl;

				RHS = ( RADr-RADl )*( 2/(gR[n+1]+gR[n]) )*(DT/DR);
				U[n][m] = U[n][m] - RHS;
				//cout << fixed << setprecision(wide) << RHS.u << endl;

				FullVec GeomSource;
				double Ptotal = (GAMMA-1)*( U[n][m].e -0.5*U[n][m].Vsq()/U[n][m].rho -0.5*U[n][m].Bsq() ) +0.5*U[n][m].Bsq();
				double Vx, Vy, By; // they will different 
				Vx = U[n][m].v*U[n][m].v/U[n][m].rho + Ptotal - U[n][m].H*U[n][m].H;
				Vy = -U[n][m].u*U[n][m].v/U[n][m].rho + U[n][m].B*U[n][m].H;
				By = U[n][m].u*U[n][m].H/U[n][m].rho - U[n][m].B*U[n][m].v/U[n][m].rho;
				GeomSource.setVec(0, Vx,Vy,0, 0,By,0, 0);
				GeomSource = GeomSource*( 2/(gR[n+1]+gR[n]) )*DT;
				U[n][m] = U[n][m] + GeomSource;
				//cout << fixed << setprecision(wide) << GeomSource.u << endl;
			}
		#endif

			RHS = ( Fz[n][m+1]-Fz[n][m] )*(DT/DZ);
			U[n][m] = U[n][m] - RHS; 

		}
	}	


	for(int n=0; n<NR+1; n++) delete [] Fr[n];
	delete [] Fr;
	for(int n=0; n<NR+1; n++) delete [] Fz[n];
	delete [] Fz;

}}



double **divCalc(FullVec **U)
{{
	double **div = new double * [NR+2];
	for(int n=0; n<NR+2; n++) div[n] = new double [NZ+2];

	double gR[NR+2]; setGridR(gR);


	for(int m=1; m<NZ+1; m++)
	{
		#ifdef CARTESIAN
			for(int n=1; n<NR+1; n++) div[n][m] = ( U[n+1][m].B - U[n-1][m].B )/(2*DR);
		#else
			for(int n=1; n<NR+1; n++) div[n][m] = ( U[n+1][m].B*gR[n+1] - U[n-1][m].B*gR[n-1] )/(2*DR)/gR[n];
		#endif
	}

	for(int n=1; n<NR+1; n++)
	{
		for(int m=1; m<NZ+1; m++) div[n][m] += ( U[n][m+1].D - U[n][m-1].D )/(2*DZ);
	}

	for(int m=1; m<NZ+1; m++) div[1][m] = 0;

	return div;

	for(int n=0; n<NR+2; n++) delete [] div[n];
	delete [] div;
}}


void DIV(FullVec **U, int st)
{
	double **divs = new double * [NR+2];
	for(int n=0; n<NR+2; n++) divs[n] = new double [NZ+2];

	divs = divCalc(U); 
	outputLayer(divs,st);

	for(int n=0; n<NR+2; n++) delete [] divs[n];
	delete [] divs;	
}