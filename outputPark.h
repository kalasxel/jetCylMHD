void outputGrid(double gR[NR+2], double gZ[NZ+2]);

void outputLayer(FullVec **U, float t, int st);
void outputLayer(double **U, int st); // for div

void printVec(FullVec &A);
void printVec(FullVec &A, FullVec &B);