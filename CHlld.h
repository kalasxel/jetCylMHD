class Hlld
{
private:

	double Cf(FullVec &A);

	double Sl, Sr, Sm;
	double SlStar, SrStar;

	void setSlr(FullVec &L, FullVec &R);
	void setSm(FullVec &L, FullVec &R);
	void setSlrStar(FullVec &L, FullVec &R);

	FullVec LS, RS;
	double PtotStar(FullVec &L, FullVec &R);
	void setLStar(FullVec &L, FullVec &R);
	void setRStar(FullVec &L, FullVec &R);
	FullVec LSS, RSS;
	void setLRStarStar(FullVec &L, FullVec &R);

	FullVec getFlStar(FullVec &L, FullVec &R);
	FullVec getFlStarStar(FullVec &L, FullVec &R);
	FullVec getFrStar(FullVec &L, FullVec &R);
	FullVec getFrStarStar(FullVec &L, FullVec &R);
	FullVec getFrStarStarDBG(FullVec &L, FullVec &R);
	FullVec getFrStarDBG(FullVec &L, FullVec &R);

public:

	FullVec getU(FullVec &AAA);
	FullVec getF(FullVec &AAA);

	FullVec getFrDBG(FullVec &L, FullVec &R);

	FullVec getFlowHLLD(FullVec &L, FullVec &R);
};
