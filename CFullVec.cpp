#include <iostream> 

#include "Param.h"
#include "CFullVec.h"

using namespace std;


FullVec::FullVec()
{{
	rho = u = v = w = B = H = D = e = 0;
}}

FullVec::FullVec(double t_rho, double t_u, double t_v, double t_w, double t_B, 
				 double t_H, double t_D, double t_e)
{{
	rho = t_rho;
	u = t_u; v = t_v; w = t_w;
	B = t_B; H = t_H; D = t_D;
	e = t_e;
}}



double FullVec::Vsq()
{{
	return u*u + v*v + w*w;
}} 

double FullVec::Bsq()
{{
	return B*B + H*H + D*D;
}} 

double FullVec::VBsq()
{{
	return u*B + v*H + w*D;
}}

// 2>3
double FullVec::P()
{{
	return (GAMMA-1)*( e - rho*Vsq()/2 - Bsq()/2 );
}}

// 2>3
double FullVec::Ptot()
{{
	return P() + Bsq()/2;
}}



void FullVec::setVec(double t_rho, double t_u, double t_v, double t_w, 
			double t_B, double t_H, double t_D, double t_e)
{{
	rho = t_rho;
	u = t_u; v = t_v; w = t_w;
	B = t_B; H = t_H; D = t_D;
	e = t_e;

}}

double FullVec::EfromP(double P)
{{
	return P/(GAMMA-1) + rho*Vsq()/2 + Bsq()/2;
}}

void FullVec::setVecP(double t_rho, double t_u, double t_v, double t_w, 
			double t_B, double t_H, double t_D, double t_P)
{{
	rho = t_rho;
	u = t_u; v = t_v; w = t_w;
	B = t_B; H = t_H; D = t_D;
	e = t_P/(GAMMA-1) + rho*( t_u*t_u+t_v*t_v+t_w*t_w )/2 + ( t_B*t_B+t_H*t_H+t_D*t_D )/2;
	u = t_rho*t_u; v = t_rho*t_v; w = t_rho*t_w; // to make U

}}





FullVec FullVec::operator+(FullVec &AAA)
{{
	FullVec TMP;

	TMP.rho = rho + AAA.rho;
	TMP.u = u + AAA.u;
	TMP.v = v + AAA.v;
	TMP.w = w + AAA.w;
	TMP.B = B + AAA.B;
	TMP.H = H + AAA.H;
	TMP.D = D + AAA.D;
	TMP.e = e + AAA.e;

	return TMP;
}}


FullVec FullVec::operator-(FullVec &AAA)
{{
	FullVec TMP;

	TMP.rho = rho - AAA.rho;
	TMP.u = u - AAA.u;
	TMP.v = v - AAA.v;
	TMP.w = w - AAA.w;
	TMP.B = B - AAA.B;
	TMP.H = H - AAA.H;
	TMP.D = D - AAA.D;
	TMP.e = e - AAA.e;

	return TMP;
}}

FullVec FullVec::operator*(double t)
{{
	FullVec TMP;

	TMP.rho = t*rho;
	TMP.u = t*u;
	TMP.v = t*v;
	TMP.w = t*w;
	TMP.B = t*B;
	TMP.H = t*H;
	TMP.D = t*D;
	TMP.e = t*e;

	return TMP;
}}



double FullVec::returnVec(int type)
{{
	switch(type)
	{
		case 0:
			return rho; break;
		case 1:
			return u/rho; break;
		case 2:
			return v/rho; break;
		case 3:
			return w/rho; break;
		case 4:
			return B; break;
		case 5:
			return H; break;
		case 6:
			return D; break;
		case 7:
			return e; break;
		case 8:
			return PFat(); break;
		case 9:
			return PtotFat(); break;
		default:
			cout << "AN ERROR IN PRINT" << endl; 
			return 0;
	}
}}


// to get the physical vector instead the U ( U contains rho*u etc. )
void FullVec::slimU()
{{
	u = u/rho;
	v = v/rho;
	w = w/rho;
}}

void FullVec::fatU()
{{
	u = u*rho;
	v = v*rho;
	w = w*rho;
}}

//2-3:
double FullVec::PFat()
{{
	return (GAMMA-1)*( e - Vsq()/(2*rho) - Bsq()/2 );
}}
double FullVec::PtotFat()
{{
	return PFat() + Bsq()/2;
}}