#include <iostream> 
#include <cmath>

#include "Param.h"
#include "CFullVec.h"
#include "CHlld.h"



double setGridR(double gR[NR+2])
{{
	for(int n=0; n<NR+2; n++)
	{
		gR[n] = (n-1)*DR;
		//std::cout << n << ' ' << gR[n] << std::endl; 
	}
}}
double setGridZ(double gZ[NZ+2])
{{
	for(int m=0; m<NZ+2; m++)
	{
		gZ[m] = (m-1)*DZ;
		//cout << gZ[m] << endl;
	}
}}

int NfromR(double R)
{{
	return static_cast<int>((R-RI)/DR)+1;
}}
int NfromZ(double Z)
{{
	return static_cast<int>((Z-ZI)/DZ)+1;
}}



	double differ(10);
	double Rho_in(7e-6), Vx_in(0), Vy_in(0), Vz_in(5e6), Bx_in(0), By_in(15e2), Bz_in(0), P_in(2e6);
	double Rho_out(7e-7), Vx_out(0), Vy_out(0), Vz_out(0), Bx_out(0), By_out(0), Bz_out(0), P_out(2e6);
	double thick=0.3*(RF-RI);

	int njm(0), njp(njm + NfromR(thick));
	int mjm(0), mjp( NfromZ((RF-RI)/7e-1) ); 



void setInital(FullVec **A)
{{

	int n(3);
	int m(3);
	for( m=mjm; m<mjp; m++)
	{
		for( n=0; n<njm; n++)
		{
			A[n][m].setVecP(Rho_out, Vx_out,Vy_out,Vy_out, Bx_out,By_out,Bz_out, P_out);
		}


		for( n; n<njp; n++)
		{
			double Bp, gR[NR+2];
			setGridR(gR);
			Bp = (By_in/thick)*gR[n];

			A[n][m].setVecP(Rho_in, Vx_in,Vy_in,Vz_in, Bx_in,Bp,Bz_in, P_in);
		}


		for( n; n<NR+1; n++)
		{
			double Bp, gR[NR+2];
			setGridR(gR);
			Bp = (By_in*thick)/gR[n];

			A[n][m].setVecP(Rho_out, Vx_out,Vy_out,Vy_out, Bx_out,Bp,Bz_out, P_out);
		}


	}
	for( m; m<NZ+1; m++)
	{
		for( n=0; n<njm; n++)
		{
			A[n][m].setVecP(Rho_out, Vx_out,Vy_out,Vy_out, Bx_out,By_out,Bz_out, P_out);
		}
		for( n; n<njp; n++)
		{
			A[n][m].setVecP(Rho_out, Vx_out,Vy_out,Vy_out, Bx_out,By_out,Bz_out, P_out);
		}
		for( n; n<NR+1; n++)
		{
			A[n][m].setVecP(Rho_out, Vx_out,Vy_out,Vy_out, Bx_out,By_out,Bz_out, P_out);
		}
	}

	for(int m=0; m<mjm; m++)
	{
		for( n=0; n<njm; n++)
		{
			A[n][m].setVecP(Rho_out, Vx_out,Vy_out,Vy_out, Bx_out,By_out,Bz_out, P_out);
		}
		for( n; n<njp; n++)
		{
			A[n][m].setVecP(Rho_out, Vx_out,Vy_out,Vy_out, Bx_out,By_out,Bz_out, P_out);
		}
		for( n; n<NR+1; n++)
		{
			A[n][m].setVecP(Rho_out, Vx_out,Vy_out,Vy_out, Bx_out,By_out,Bz_out, P_out);
		}
	}


/*
	double gR[NR+2]; setGridR(gR);
	double gZ[NZ+2]; setGridZ(gZ);
	for(int n=0	; n<NR+2; n++)
	{
		for(int m=0; m<NZ+2; m++)
		{
			A[n][m].setVecP(Rho_out, Vx_out,Vy_out,Vy_out, Bx_out,By_out,Bz_out, P_out);

			//A[n][m].B = gR[n]*gZ[m];
			//A[n][m].D = -gZ[m]*gZ[m];
			//A[n][m].B = gR[n]*gR[n]*gZ[m];
			//A[n][m].D = -gZ[m]*gZ[m]*gR[n];

			int alpha = 0;
			A[n][m].B = pow(gR[n],alpha+1)*pow(gZ[m],alpha+1);
			A[n][m].D = -pow(gR[n],alpha)*pow(gZ[m],alpha+2);
		}
	}
*/

}}

void setBound(FullVec **U)
{{
	FullVec R1, R2, Z1, Z2;

	R1.setVecP(Rho_out, Vx_out,Vy_out,Vy_out, Bx_out,By_out,Bz_out, P_out);
	R2=R1; Z2=R1; Z1=R1;

	for(int n=0; n<NR+2; n++)
	{
		U[n][0] = Z1;
		U[n][NZ+1] = Z2;
	}
	for(int m=0; m<NZ+2; m++)
	{
		U[0][m] = R1;
		U[NR+1][m] = R2;
	}
}}


