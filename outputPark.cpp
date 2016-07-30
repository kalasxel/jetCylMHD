#include <iostream> 
#include <fstream>
#include <iomanip>

#include "Param.h"
#include "CFullVec.h"

using namespace std;


#define GRIDR "results/GRIDR.dat"
#define GRIDZ "results/GRIDZ.dat"
const char *NAMES[] = 
	{"results/0_Ro/Ro_%d.dat",
	"results/1_Vx/Vx_%d.dat",
	"results/2_Vy/Vy_%d.dat",
	"results/3_Vz/Vz_%d.dat",
	"results/4_Bx/Bx_%d.dat",
	"results/5_By/By_%d.dat",
	"results/6_Bz/Bz_%d.dat",
	"results/7_E/E_%d.dat",
	"results/8_P/P_%d.dat",
	"results/9_Ptot/Ptot_%d.dat"};
#define TIME "results/TIME.dat"
#define DIV "results/10_div/Div_%d"



void outputLayer(FullVec **U, float t, int st)
{{
	ofstream times;
	times.open(TIME, ios::app);
	times << fixed << setprecision(4) << t << endl;
	times.close();

	for(int i=0; i<10; i++)
	{
		char fullName[30];
		sprintf(fullName,NAMES[i], st);

		ofstream flow;
		flow.open(fullName, ios_base::out | ios_base::trunc);

		for(int m=1; m<NZ+1; m++)
		{
			for(int n=1; n<NR+1; n++)
			{
				flow << fixed << setprecision(wide) << showpos << U[n][m].returnVec	(i) << "  ";
			}
			flow << endl;
		}

		flow.close();		
	}
}}


void outputLayer(double **U, int st) // for div
{{
		char fullName[30];
		sprintf(fullName,DIV, st);

		ofstream flow;
		flow.open(fullName, ios_base::out | ios_base::trunc);

		for(int m=1; m<NZ+1; m++)
		{
			for(int n=1; n<NR+1; n++)
			{
				flow << fixed << setprecision(wide) << showpos << U[n][m] << "  ";
			}
			flow << endl;
		}

		flow.close();		
}}

void outputGrid(double gR[NR+2], double gZ[NZ+2])
{{
	ofstream gridR;
	gridR.open(GRIDR, ios_base::out | ios_base::trunc);
	for(int n=1; n<NR+1; n++) gridR << gR[n] << endl;
	gridR.close();

	ofstream gridZ;
	gridZ.open(GRIDZ, ios_base::out | ios_base::trunc);
	for(int m=1; m<NZ+1; m++) gridZ << gZ[m] << endl;
	gridZ.close();
}}


void printVec(FullVec &A)
{{
	cout << A.rho << endl;
	cout << A.u << endl;
	cout << A.v << endl;
	cout << A.w << endl;
	cout << A.B << endl;
	cout << A.H << endl;
	cout << A.D << endl;
	cout << A.e << endl;
}}

void printVec(FullVec &A, FullVec &B)
{{
	cout << fixed << setprecision(wide) << A.rho << ' ' << B.rho << endl;
	cout << fixed << setprecision(wide) << A.u << ' ' << B.u << endl;
	cout << fixed << setprecision(wide) << A.v << ' ' << B.v << endl;
	cout << fixed << setprecision(wide) << A.w << ' ' << B.w << endl;
	cout << fixed << setprecision(wide) << A.B << ' ' << B.B << endl;
	cout << fixed << setprecision(wide) << A.H << ' ' << B.H << endl;
	cout << fixed << setprecision(wide) << A.D << ' ' << B.D << endl;
	cout << fixed << setprecision(wide) << A.e << ' ' << B.e << endl;
}}

