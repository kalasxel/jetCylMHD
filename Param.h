#define PI 3.1415926535897932384626433832795
//#define AvNum 6.022140857e23 // 1/mol
//#define Boltz 8.6173324e-5 // eV/K
#define eVtoErg 1.6021766208e-12
#define eVtoK 1.16045221e4

#define GAMMA 1.66666666666666666666666667
#define Mnuclon 3.3211e-23

#define VISCISITY1 0.001

template <class XXX> int SIGNUM(XXX val){ return val>0 ? 1 : 0; }
//#define SIGNUM(Value) ((Value) < 0 ? (-1) : !!(Value))

#define RI 0.0
#define RF 5
#define DR 1e-1 // "make clear" after each changind of step
const int NR = static_cast<int>( (RF-RI)/DR ) + 1; 


#define ZI 0.0
#define ZF 40
#define DZ 1e-1
const int NZ = static_cast<int>( (ZF-ZI)/DZ ) + 1; 

#define DT 1e-8 // "make clear" after each changind of step
#define TFIN 10e-6
const int ST = static_cast<int>( TFIN/DT );
const int CVAR = static_cast<int>( ST/10 );



#define wide 3