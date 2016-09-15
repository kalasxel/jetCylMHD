class FullVec
{
public:

	FullVec();
	FullVec(double t_rho, double t_u, double t_v, double t_w, 
			double t_B, double t_H, double t_D, double t_e);

	double rho;
	double u, v, w;  // Vx, Vy, Vz respectively
	double B, H, D;  // Bx, By, Bz respectively
	double e;

	double Vsq();
	double Bsq();
	double VBsq();
	double getT();

	double P();
	double Ptot();

	void setVec(double t_rho, double t_u, double t_v, double t_w, 
			double t_B, double t_H, double t_D, double t_e);
	void setVecT(double nn, double t_u, double t_v, double t_w,
				double t_B, double t_H, double t_D, double T);



	FullVec operator+(FullVec &AAA);
	FullVec operator-(FullVec &AAA);
	FullVec operator*(double t); // multiply from the right!!!!!!!


	double returnVec(int type);
	double returnVecT(int type);

	void slimU(); // to get the physical vector
	void fatU();

};