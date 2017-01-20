#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include "dco.hpp"

using namespace dco;
typedef ga1s<double> DCO_MODE;
typedef DCO_MODE::type DCO_TYPE;
typedef DCO_MODE::tape_t DCO_TAPE_TYPE;

using namespace std;
int Nsize;
#define Nsize 16
/*Function for assigning values to 2D Array (Preferably for Initial Conditions)*/
template <typename TYPE>
void assign(TYPE **a,int row,int col,double val)
{
    int i,j;
    for (i=0;i<row;i++){
      for (j=0;j<col;j++){
        a[i][j]=val;}}
}

/*Fuction for assigning boundary Conditions*/
template <typename TYPE>
void boundary(TYPE **a,int row,int col,double val,int bound){
  int i,num;
  if(bound==1 || bound ==3){
    if(bound==1){num=0;}
    else{num=col-1;}
  for (i=0;i<row;i++){
    a[i][num]=val;
  }}
  if(bound==2 || bound ==4){
    if(bound==2){num=0;}
    else{num=row-1;}
  for (i=0;i<col;i++){
    a[num][i]=val;
  }}
}
/*Function for equating two 2D Arrays*/
template <typename TYPE>
void equate(TYPE **a,TYPE **b,int row,int col){
  int i,j;
  for(i=0;i<row;i++){
    for(j=0;j<col;j++){
      a[i][j]=b[i][j];}}
}
void Cvalue(c){
	int i,j;
	double amp=10;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			c[i][j]=amp+amp*pow(sin(2*M_PI*i/N),2)*pow(cos(2*M_PI*j/N),2);
		}
	}
}
template <typename TYPE>
inline void simulation(const TYPE& k, int N, double endT, TYPE *const U)
{
  int i,j;
  //Initialization of a 2D Array
  TYPE **U_old= new TYPE*[N];
  for(int i = 0; i < N; i++){
		U_old[i] = new TYPE[N];}
  TYPE **U_new= new TYPE*[N];
  for(int i = 0; i < N; i++){
  	U_new[i] = new TYPE[N];}
  /*Initial Condition*/
  assign(U_old,N,N,25);

  /*  2D - Boundary numbering
      _________B3________
     |                   |
     |                   |
     |                   |
     |                   |
     B2       Grid       B4
     |                   |
     |                   |
     |                   |
     |_________B1________|
  */
  /*Boundary Condition*/
  boundary(U_old,N,N,100,1);
  boundary(U_old,N,N,100,2);
  boundary(U_old,N,N,100,4);
  boundary(U_old,N,N,100,3);
  equate(U_new,U_old,N,N);
  //double k=20;//Thermal Conductivity
  double dt=0.005;//(Half) Time step
  double dx=1.0/N-1;//stepping in x direction
  double dy=1.0/N-1;//stepping in y direction
  double DT;//Total time
  /* R_x and R_y */
  TYPE Rx=k*dt/ (dx * dx);
  TYPE Ry=k*dt/ (dy * dy);
  TYPE a[N],b[N],c[N],d[N],m;

  while(DT<endT){//end time
  for (j=1;j<N-1;j++){
    std::fill_n(a,N,-Rx/2.0);
    std::fill_n(b,N,(1.0+Rx));
    std::fill_n(c,N,-Rx/2.0);
    d[1]=(Ry/2.0)*(U_old[1][j-1]+U_old[1][j+1]) + (1.0-Ry)*(U_old[1][j]) - a[1]*U_old[0][j];
    for (i=2;i<=N-2;i++){/*Thomas Algorithm*/
        m=a[i]/b[i-1];
        b[i]=b[i] - (m*c[i-1]);
        if (i==N-2){d[i]= (Ry /2.0)*(U_old[i][j-1]+U_old[i][j+1]) + (1.0-Ry)*(U_old[i][j]) - c[N-1]*U_old[N-1][j];}
        else{d[i]= (Ry /2.0)*(U_old[i][j-1]+U_old[i][j+1]) + (1.0-Ry)*(U_old[i][j]);}
        d[i]=d[i]-(m*d[i-1]);
    }
    U_new[N-2][j]=d[N-2]/b[N-2];
    for (i=N-3;i>=1;i--){ /*back Substitution*/
      U_new[i][j]=(d[i]-(c[i]*U_new[i+1][j]))/b[i];
    }
  }
equate(U_old,U_new,N,N);
for (i=1;i<N-1;i++){
  std::fill_n(a,N,-Ry/2.0);
  std::fill_n(b,N,(1.0+Ry));
  std::fill_n(c,N,-Ry/2.0);
  d[1]=(Rx/2.0)*(U_old[i-1][1]+U_old[i+1][1]) + (1-Rx)*(U_old[i][1]) - a[1]*U_old[i][0];
  for (j=2;j<=N-2;j++){//Thomas Algorithm//
      m=a[j]/b[j-1];
      b[j]=b[j] - (m*c[j-1]);
      if (j==N-2){d[j] = (Rx/2.0)*(U_old[i-1][j]+U_old[i+1][j])+(1-Rx)*(U_old[i][j]) - c[N-1]*U_old[i][N-2];}
      else{d[j] = (Rx/2.0)*(U_old[i-1][j]+U_old[i+1][j])+(1-Rx)*(U_old[i][j]);}
      d[j]=d[j] - (m*d[j-1]);

  }
  U_new[i][N-2]=d[N-2]/b[N-2];
  for (j=N-3;j>=1;j--){ //back Substitution
    U_new[i][j]=(d[j]-(c[j]*U_new[i][j+1]))/b[j];
  }
}
equate(U_old,U_new,N,N);
DT=DT+dt*2;
}
int point,count=-1;
  for(i=0;i<N*N;i++){
    point=i%(N);
    if(point==0){count=count+1;}
    U[i]=U_old[count][point];
  }

  //for printing a 2D Array
   /* for (i=0;i<N;i++){
      for (j=0;j<N;j++){
        cout<<U_old[i][j];cout<<"\t";}
        cout<<endl;}*/
  delete[] U_old,U_new;
}

template <typename TYPE>
inline TYPE objective_function(const TYPE& k, int N, double endT, TYPE *const U,const double *const O)
{
	simulation(k,N,endT,U);
	
	TYPE objective_function = 0;
	for(int i=0;i<N*N-1;++i)
		objective_function += (U[i]-O[i])*(U[i]-O[i]);
	objective_function/=((N*N)-1);

	return objective_function;
}

void derivative_a1s(const double& k, int N, double endT, const double *const U, const double *const O, double &obj, double &grad_obj)
{
DCO_MODE::global_tape=DCO_TAPE_TYPE::create();
DCO_TYPE *a1s_U = new DCO_TYPE[N*N];
for (int i=0;i<N*N;++i) a1s_U[i] = U[i];
DCO_TYPE a1s_k = k;
DCO_MODE::global_tape->register_variable(a1s_k);
DCO_TYPE a1s_obj = objective_function(a1s_k,N,endT,a1s_U,O);
DCO_MODE::global_tape->register_output_variable(a1s_obj);

derivative(a1s_obj) = 1;

DCO_MODE::global_tape->interpret_adjoint();
grad_obj = derivative(a1s_k);
obj = value(a1s_obj);
DCO_MODE::global_tape->interpret_adjoint();
DCO_TAPE_TYPE::remove(DCO_MODE::global_tape);

delete [] a1s_U;

}

		
int main(int argc, char* argv[])
{
  double *const U = new double [Nsize*Nsize];
  double *const O = new double [Nsize*Nsize];
  double *const local_U = new double [Nsize*Nsize];
  double endT=30, k=1.0;
  simulation(k,Nsize,endT,O);
  k=2.0;
  //steepest descent
  double epsilon = 1e-2;
  double gradient_obj=epsilon;
  int it = 0;
  
  /*simulation(k,Nsize,endT,U);
  double obj = objective_function(k,Nsize,endT,U,O);*/
  while (fabs(gradient_obj)>=epsilon){
	double obj;
	derivative_a1s(k,Nsize,endT,U,O,obj,gradient_obj);
	cout<<"Gradient_obj="<<gradient_obj<<endl;

	double local_obj = obj, alfa = 0.01, local_k;
	while (local_obj>=obj){
		local_k = k-alfa*gradient_obj;
		cout<<"k="<<local_k<<endl;
		for(int i=0;i<Nsize*Nsize;++i) local_U[i] = U[i];
		local_obj = objective_function(local_k,Nsize,endT,local_U,O);
		cout<<"local_obj="<<local_obj<<endl;
		alfa/=2.0;
	}
  	k = local_k;
	cout<<"Thermal diffusivity = "<<k<<endl;
  }

  delete[] U;
  delete[] O;
  delete[] local_U;

  return 0;
}
  

