#define PI 3.1415926535897932384626433832795
#define GAMMA 1.66666666666666666666666667

template <class XXX> int SIGNUM(XXX val){ return val>0 ? 1 : 0; }
//#define SIGNUM(Value) ((Value) < 0 ? (-1) : !!(Value))

#define RI 0.0
#define RF 5
#define DR 5e-2 // "make clear" after each changind of step
const int NR = static_cast<int>( (RF-RI)/DR ) + 1; 


#define ZI 0.0
#define ZF 40
#define DZ 5e-2
const int NZ = static_cast<int>( (ZF-ZI)/DZ ) + 1; 

#define DT 2e-9 // "make clear" after each changind of step
#define TFIN 15e-6
const int ST = static_cast<int>( TFIN/DT );
const int CVAR = static_cast<int>( ST/10 );


#define wide 12