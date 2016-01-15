#include "common.h"
class dataproc
{
public:
	//the values of the "M" distribution at the grid points are stored in m_dis(I)(I=1,2...NVAR)
	//and the cluster concentrations are stored in m_dis(I+NVAR)(I=1,2...NDISCR)
	float m_dis[COUNT_AX];
	//the particle sizes at the grid points and for the cluster sizes are stored in p_size
	float p_size[COUNT_AX];
	//W(I) contains the "W" value at the grid points.
	float w[COUNT_W];
	//
	float heading[INPUT_CASE_COUNT];
	//store ten user-defined switches which are read in as input data
	int iparm[INPUT_CASE_COUNT];
	//case number
	int cscnt;
	//size of smallest particle
	float xa;
	//size of largest particle
	float xb;
	//reletive error(see ode code write up)
	float relerr;
	//absolute integration error (see ode code write up)
	float abserr;
	//
	int nvar;
	//maximum number of grid points
	float nmxvar;
	//number of discrete sizes
	float ndiscr;
	//number of discrete multiplets given in the continuous regime
	float nmulgv;
	//number of outputs
	float nprnt;
	//the output number for which the program is to switch interpolation methods
	int iswtch;
	//output number at which the steady state approximation is to be used
	int isstat;
	//number of grid points given
	int igrid;
	//number of interpolated points between grid points 
	int nptbet;
	//number of continuous sources
	int ncsorc;
	//number of continuous removals
	int ncremv;
	//number of discrete sources
	int ndsorc;
	//number of discrete removals
	int ndremv;
	//0: no condensation / 1: with condensation
	int icondn;
	//flag to indicate 
	int isplin;
	//flag to call setup before interpolating (0:no/1:yes)
	int isetup;
	//flag to call setup for switched interpolation method
	int ispswt;
	//output number to switch interpolation method
	int isetsw;
	//flag for time dependent coagulation and evaporation coefficients
	int itimdp;
	//the i-th quadrature weight factor for integration over the entire continuous domain
	int factor[24];
	
	int issflg;
	
	//number of quadrature points ( must be an even number )
	int nquadp;
	//half the number of quadrature points
	int nquad2;	

	float tout[16];

	int ndcls[20];

	int nmulp1;
	int nmulp2;
	int igdmp1;
	int igdmp2;

	int nvarm1;
	int nvarp1;
	int nvarp2;
	int nvart;

	float xboxa;
	float algxba;

	float a[60];
	float x[60];
public:
	int readCaseCnt(ifstream* in);
	void readCase(ifstream* in);
};
