#include <iostream> 
#include <cmath>

#include "Param.h"
#include "CFullVec.h"
#include "CHlld.h"

//using namespace std;


FullVec Hlld::getU(FullVec &AAA)
{{
	FullVec TMP;

	TMP.rho = AAA.rho;
	TMP.u = AAA.rho*AAA.u;
	TMP.v = AAA.rho*AAA.v;
	TMP.w = AAA.rho*AAA.w;
	TMP.B = AAA.B;
	TMP.H = AAA.H;
	TMP.D = AAA.D;
	TMP.e = AAA.e;

	return TMP;
}}


FullVec Hlld::getF(FullVec &AAA)
{{
	FullVec TMP;

	TMP.rho = AAA.rho*AAA.u;
	TMP.u = AAA.rho*AAA.u*AAA.u + AAA.Ptot() - AAA.B*AAA.B;
	TMP.v = AAA.rho*AAA.v*AAA.u - AAA.B*AAA.H;
	TMP.w = AAA.rho*AAA.w*AAA.u - AAA.B*AAA.D;
	TMP.B = 0;
	TMP.H = AAA.H*AAA.u - AAA.B*AAA.v;
	TMP.D = AAA.D*AAA.u - AAA.B*AAA.w;
	TMP.e = ( AAA.e + AAA.Ptot() )*AAA.u - AAA.B*AAA.VBsq();

	return TMP;
}}


// 3
double Hlld::Cf(FullVec &A)
{{
	double sq = GAMMA*A.P() + A.Bsq();
	
	return sqrt( sq + sqrt(sq*sq-4*GAMMA*A.P()*A.B*A.B) )/sqrt(2*A.rho);
}}



// 67
void Hlld::setSlr(FullVec &L, FullVec &R)
{{
	//Sl = min(L.u,R.u) - 1.5*max(Cf(L),Cf(R));
	//Sr = max(L.u,R.u) + 1.5*max(Cf(L),Cf(R));

	//Sl = min(L.u,R.u) - 3*max(Cf(L),Cf(R));
	//Sr = max(L.u,R.u) + 3*max(Cf(L),Cf(R));

	Sl = std::min(L.u,R.u) - std::max(Cf(L),Cf(R)); //Sl *= 0.5;
	Sr = std::max(L.u,R.u) + std::max(Cf(L),Cf(R)); //Sr *= 2;

	//Sl = min( (L.u-Cf(L)), (R.u-Cf(R)) );
	//Sr = max( (L.u+Cf(L)), (R.u+Cf(R)) );
}}

// 38
void Hlld::setSm(FullVec &L, FullVec &R)
{{
	setSlr(L,R);

	double lU = (Sl-L.u)*L.rho;
	double rU = (Sr-R.u)*R.rho;

	Sm = ( rU*R.u - lU*L.u - R.Ptot() + L.Ptot() )/( rU - lU );
}}

// 41
double Hlld::PtotStar(FullVec &L, FullVec &R)
{{
	double lU = (Sl-L.u)*L.rho;
	double rU = (Sr-R.u)*R.rho;
	return ( rU*L.Ptot() - lU*R.Ptot() + rU*lU*(R.u-L.u) )/(rU-lU);
}}


// 39-48
void Hlld::setLStar(FullVec &L, FullVec &R)
{{
	FullVec LR;
	double Sa;

	setSm(L,R);

	LR=L; 
	Sa=Sl;

	LS.rho = LR.rho*(Sa-LR.u)/(Sa-Sm);
	LS.u = Sm;

	double down = LR.rho*(Sa-LR.u)*(Sa-Sm) - LR.B*LR.B;

	LS.v = LR.v - LR.B*LR.H*(Sm-LR.u)/down;
	LS.w = LR.w - LR.B*LR.D*(Sm-LR.u)/down;

	LS.B = LR.B; // it does not matter
	LS.H = LR.H*( LR.rho*(Sa-LR.u)*(Sa-LR.u) - LR.B*LR.B)/down;
	LS.D = LR.D*( LR.rho*(Sa-LR.u)*(Sa-LR.u) - LR.B*LR.B)/down;
/*
	double lU = (Sl-L.u)*L.rho;
	double rU = (Sr-R.u)*R.rho;
	tPtotLS = ( rU*L.Ptot() - lU*R.Ptot() + rU*lU*(R.u-L.u) )/(rU-lU);
*/
	LS.e = ( (Sa-LR.u)*LR.e - LR.Ptot()*LR.u + PtotStar(L,R)*Sm + LR.B*(LR.VBsq()-LS.VBsq()) )/(Sa-Sm);
}}
void Hlld::setRStar(FullVec &L, FullVec &R)
{{
	FullVec LR;
	double Sa;

	setSm(L,R);

	LR=R; 
	Sa=Sr;

	RS.rho = LR.rho*(Sa-LR.u)/(Sa-Sm);
	RS.u = Sm;

	double down = LR.rho*(Sa-LR.u)*(Sa-Sm) - LR.B*LR.B;

	RS.v = LR.v - LR.B*LR.H*(Sm-LR.u)/down;
	RS.w = LR.w - LR.B*LR.D*(Sm-LR.u)/down;

	RS.B = LR.B; // it does not matter
	RS.H = LR.H*( LR.rho*(Sa-LR.u)*(Sa-LR.u) - LR.B*LR.B)/down;
	RS.D = LR.D*( LR.rho*(Sa-LR.u)*(Sa-LR.u) - LR.B*LR.B)/down;
/*
	double lU = (Sl-L.u)*L.rho;
	double rU = (Sr-R.u)*R.rho;
	tPtotRS = ( rU*L.Ptot() - lU*R.Ptot() + rU*lU*(R.u-L.u) )/(rU-lU);
*/
	RS.e = ( (Sa-LR.u)*LR.e - LR.Ptot()*LR.u + PtotStar(L,R)*Sm + LR.B*(LR.VBsq()-RS.VBsq()) )/(Sa-Sm);
}}

// 51
void Hlld::setSlrStar(FullVec &L, FullVec &R)
{{
	setLStar(L,R); setRStar(L,R);

	SlStar = Sm - fabs(LS.B)/sqrt(LS.rho);
	SrStar = Sm + fabs(RS.B)/sqrt(RS.rho);
}}


// 49-63
void Hlld::setLRStarStar(FullVec &L, FullVec &R) 
{{
	setSlrStar(L,R);

	LSS.rho=LS.rho;
	RSS.rho=RS.rho;

	LSS.u = RSS.u = Sm;
	LSS.v = RSS.v = ( sqrt(LS.rho)*LS.v + sqrt(RS.rho)*RS.v + (RS.H-LS.H)*SIGNUM(LS.B) )/( sqrt(LS.rho)+sqrt(RS.rho) );
	LSS.w = RSS.w = ( sqrt(LS.rho)*LS.w + sqrt(RS.rho)*RS.w + (RS.D-LS.D)*SIGNUM(LS.B) )/( sqrt(LS.rho)+sqrt(RS.rho) );

	LSS.B = RSS.B = LS.B; // No Matter
	LSS.H = RSS.H = ( sqrt(LS.rho)*RS.H + sqrt(RS.rho)*LS.H + sqrt(LS.rho*RS.rho)*(RS.v-LS.v)*SIGNUM(LS.B) )/( sqrt(LS.rho)+sqrt(RS.rho) );
	LSS.D = RSS.D = ( sqrt(LS.rho)*RS.D + sqrt(RS.rho)*LS.D + sqrt(LS.rho*RS.rho)*(RS.w-LS.w)*SIGNUM(LS.B) )/( sqrt(LS.rho)+sqrt(RS.rho) );

	LSS.e = LS.e - sqrt(LS.rho)*SIGNUM(LS.B)*(LS.VBsq() - LSS.VBsq());
	RSS.e = RS.e + sqrt(RS.rho)*SIGNUM(RS.B)*(RS.VBsq() - RSS.VBsq());
}}

////////////////////////////////////////////////////////////////////////////////////////////

FullVec Hlld::getFlStar(FullVec &L, FullVec &R)
{{
	FullVec TMP;
	FullVec Ul, UlStar, Fl;

	setLStar(L,R); setRStar(L,R);

	Fl=getF(L);

	Ul = getU(L); UlStar = getU(LS);
	Ul = Ul*Sl; UlStar = UlStar*Sl;

	TMP = Fl - Ul + UlStar;

	return TMP;
}}


FullVec Hlld::getFlStarStar(FullVec &L, FullVec &R)
{{
	FullVec TMP;
	FullVec UlStar, UlStarStar, FlStar;

	FlStar = getFlStar(L,R);
	setLRStarStar(L,R);

	UlStar = getU(LS); UlStarStar = getU(LSS);
	UlStar = UlStar*SlStar; UlStarStar = UlStarStar*SlStar;

	TMP = FlStar - UlStar + UlStarStar;

	return TMP;
}}


FullVec Hlld::getFrStarStarDBG(FullVec &L, FullVec &R)
{{
	FullVec TMP;
	FullVec UlStarStar, UrStarStar, FlStarStar;

	FlStarStar = getFlStarStar(L,R);

	UlStarStar = getU(LSS); UrStarStar = getU(RSS);
	UlStarStar = UlStarStar*Sm; UrStarStar = UrStarStar*Sm;

	TMP = FlStarStar - UlStarStar + UrStarStar;

	return TMP;
}}


FullVec Hlld::getFrStarDBG(FullVec &L, FullVec &R)
{{
	FullVec TMP;
	FullVec UrStarStar, UrStar, FrStarStar;

	FrStarStar = getFrStarStarDBG(L,R);

	UrStarStar = getU(RSS); UrStar = getU(RS);
	UrStarStar = UrStarStar*SrStar; UrStar = UrStar*SrStar;

	TMP = FrStarStar - UrStarStar + UrStar;

	return TMP;
}}


FullVec Hlld::getFrDBG(FullVec &L, FullVec &R)
{{
	FullVec TMP;
	FullVec UrStar, Ur, FrStar;

	FrStar = getFrStarDBG(L,R);

	UrStar = getU(RS); Ur = getU(R);
	UrStar = UrStar*Sr; Ur = Ur*Sr;

	TMP = FrStar - UrStar + Ur;

	return TMP;
}}


FullVec Hlld::getFrStar(FullVec &L, FullVec &R)
{{
	FullVec TMP;
	FullVec Ur, UrStar, Fr;

	setLStar(L,R); setRStar(L,R);

	Fr = getF(R);

	Ur = getU(R); UrStar = getU(RS);
	Ur = Ur*Sr; UrStar = UrStar*Sr;

	TMP = Fr + UrStar - Ur;

	return TMP;
}}


FullVec Hlld::getFrStarStar(FullVec &L, FullVec &R)
{{
	FullVec TMP;
	FullVec UrStar, UrStarStar, FrStar;

	FrStar = getFrStar(L,R);
	setLRStarStar(L,R);

	UrStar = getU(RS); UrStarStar = getU(RSS);
	UrStar = UrStar*SrStar; UrStarStar = UrStarStar*SrStar;

	TMP = FrStar + UrStarStar - UrStar;

	return TMP;
}}


#include <cstdlib>
FullVec Hlld::getFlowHLLD(FullVec &L, FullVec &R)
{{
	L.slimU(); R.slimU();

	FullVec TMP;

	setSlr(L,R); setSlrStar(L,R); setSm(L,R);

	//Sl = fabs(Sl); // DBG

	if(Sl>0){
		TMP = getF(L);
		//TMP = getU(L); // DBG
		//std::cout << "Choice 1" << std::endl;
	}
	else if(SlStar>=0){
		TMP = getFlStar(L,R);
		//std::cout << "Choice 2" << std::endl;
	}
	else if(Sm>=0){
		TMP = getFlStarStar(L,R);
		//std::cout << "Choice 3" << std::endl;
	}
	else if(SrStar>=0){ 
		TMP = getFrStarStar(L,R);
		//std::cout << "Choice 4" << std::endl;
	}
	else if(Sr>=0){ 
		TMP = getFrStar(L,R);
		//std::cout << "Choice 5" << std::endl;
	}
	else if(Sr<0){
		TMP = getF(R);
		//std::cout << "Choice 6" << std::endl;
	}
	else
	{
		std::cout << "AN ERROR IN THE SOLVER" << std::endl;
		std::cout << Sl << ' ' << SlStar << ' ' << Sm << ' ' << SrStar << ' ' << Sr << std::endl;
		exit(4);
	}

	return TMP;
}}

