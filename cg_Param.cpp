#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include "dco.hpp"
using namespace std;
using namespace dco;
typedef ga1s<double> DCO_MODE;
typedef DCO_MODE::type DCO_TYPE;
typedef DCO_MODE::tape_t DCO_TAPE_TYPE;
int Nsize;
#define Nsize 32
int nc;
#define nc 9

/*Function for assigning values to 2D Array (Preferably for Initial Conditions)*/
template <typename TYPE>
inline void assign(TYPE **a,int row,int col,double val)
{
    int i,j;
    for (i=0;i<row;i++){
      for (j=0;j<col;j++){
        a[i][j]=val;}}
}

template <typename TYPE>
inline void cdist(TYPE *k, int cvals, int size, TYPE *c)
{
	for(int i=0;i<cvals-1;i++){
		for(int j=i*ceil(size/cvals);j<(i+1)*ceil(size/cvals);j++){
			c[j] = k[i];
		}
	}
	for(int i=(cvals-1)*ceil(size/cvals);i<size*size;i++)c[i] = k[cvals-1];
}
	
	
/*Function for assigning boundary Conditions*/
template <typename TYPE>
inline void boundary(TYPE **a,int row,int col,double val,int bound){
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
inline void equate(TYPE **a,TYPE **b,int row,int col){
  int i,j;
  for(i=0;i<row;i++){
    for(j=0;j<col;j++){
      a[i][j]=b[i][j];}}
}
template <typename TYPE>
inline void multiply(TYPE *u, TYPE  *ans,TYPE *c, TYPE  rx, TYPE  ry, int size){
int i;
for(i=0;i<size*size;i++){
  if(i<size+1 || (i+1)%size<2|| i>(size*size - size-2)){ans[i]=u[i]*(1);}
  else{ans[i]=u[i]*(1 + 2*rx*c[i] + 2*ry*c[i]) + u[i+1]*(-rx*c[i+1]) + u[i-1]*(-rx*c[i-1]) + u[i-size]*(-ry*c[i-size]) + u[i+size]*(-ry*c[i+size]);}
}}
template <typename TYPE>
inline TYPE dotpr(TYPE *u,TYPE *v,int size){
  int i;
  TYPE x;
  x=0;
  for(i=0;i<size;i++){
    x=x+(u[i]*v[i]);
  }
  return x;
}
template <typename TYPE>
inline void simulation(TYPE *k, int cvals, int size, double endT, TYPE *U)
{
  int i,j;
  TYPE  *c;
  c = new TYPE [size*size];
  cdist(k,cvals,size,c);
  /*for(i=0;i<ceil(size/16)*(size);i++)c[i] = k[0];
  for(i=ceil(size/10)*(size);i<2*ceil(size/4)*(size);i++)c[i] = k[1];
  for(i=2*ceil(size/4)*(size);i<3*ceil(size/4)*(size);i++)c[i] = k[2];
  for(i=3*ceil(size/4)*(size);i<size*size;i++)c[i] = k[3];*/
  
  //for(i=0;i<size*size;i++)cout<<c[i]<<endl;

  //Initialization of a 2D Array
  TYPE **U_old= new TYPE*[size];
  for(i = 0; i < size; i++){
		U_old[i] = new TYPE[size];}
  /*Initial Condition*/
  assign(U_old,size,size,25.0);
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
  boundary(U_old,size,size,100.0,1);
  boundary(U_old,size,size,100.0,2);
  boundary(U_old,size,size,100.0,4);
  boundary(U_old,size,size,100.0,3);
  //equate(U_new,U_old,size,size);
  double dt=0.05;//(Half) Time step
  double dx=1.0/size-1;//stepping in x direction
  double dy=1.0/size-1;//stepping in y direction
  double DT;//Total time
  /* R_x and R_y */
  TYPE  Rx=dt/ (dx * dx);
  TYPE  Ry=dt/ (dy * dy);
  //converting to AX=B form
  TYPE  *X,*B,*C;
  X = new TYPE [size*size];
  B = new TYPE [size*size];
  C = new TYPE [size*size];
  int point,count;
  count=-1;
  for(i=0;i<size*size;i++){
    point=i%(size);
    if(point==0){count=count+1;}
    B[i]=U_old[count][point];
    X[i]=0.0;
  }
  TYPE  *res,*Ax,*ADk,*d,alfa,beta,rtr,rtr_temp,tol;
  tol=1e-9;
  res = new TYPE [size*size];
  ADk = new TYPE [size*size];
  Ax =new TYPE [size*size];
  d = new TYPE [size*size];
  while(DT<endT){//end time
  //Conjugate gradient
  //multiply(X,Ax,c,Rx,Ry,size);
  for(i=0;i<size*size;i++){
    res[i]=d[i]=B[i];}
    rtr=dotpr(res,res,size*size);
  for(j=0;j<size*size;j++){
    multiply(d,ADk,c,Rx,Ry,size);
    alfa=rtr/dotpr(ADk,d,size*size);
    for(i=0;i<size*size;i++){
      res[i]=res[i] - alfa*ADk[i];
      X[i]=X[i] + alfa*d[i];
    }
    rtr_temp=rtr;
    rtr=dotpr(res,res,size*size);
    if(sqrt(rtr)<tol){break;}
    beta=rtr/rtr_temp;
    for(i=0;i<size*size;i++){
    d[i]=res[i] + beta*d[i];}
  }
  for(i=0;i<size*size;i++){
      B[i]=X[i];
      X[i]=0.0;
  }
  DT=DT+dt*2;
  }

  for(i=0;i<size*size;i++){
      U[i]=B[i];
  }
  delete[] U_old;
}

template <typename TYPE>
inline TYPE objective_function(TYPE *k, int cvals, int size, double endT, TYPE *U,const double *const O)
{
  simulation(k,cvals,size,endT,U);
  TYPE objective_function = 0;
	for(int i=0;i<size*size-1;++i)
    objective_function += (U[i]-O[i])*(U[i]-O[i]);
  objective_function/=((size*size)-1);


	return objective_function;
}

void derivative_a1s(const double *const k, int cvals, int N, double endT,const double *const U,const double *const O, double &obj, double *grad_obj)
{
DCO_MODE::global_tape=DCO_TAPE_TYPE::create();
DCO_TYPE *a1s_U = new DCO_TYPE[N*N];
DCO_TYPE *a1s_k = new DCO_TYPE[cvals];
//DCO_MODE::global_tape->reset();
for (int i=0;i<N*N;++i){a1s_U[i] = U[i];}
for(int i=0;i<cvals;i++){a1s_k[i] = k[i];}
//DCO_MODE::global_tape->reset();
for(int i=0;i<cvals;i++){
DCO_MODE::global_tape->register_variable(a1s_k[i]);}
DCO_TYPE a1s_obj=objective_function(a1s_k,cvals,N,endT,a1s_U,O);
//cout<<a1s_obj<<endl;
derivative(a1s_obj) = 1.0;
DCO_MODE::global_tape->interpret_adjoint();
for(int i=0;i<cvals;i++){
//DCO_MODE::global_tape->register_output_variable(a1s_obj);
grad_obj[i] = derivative(a1s_k[i]);}
//cout<<grad_obj[10]<<endl;
//for (int i=0;i<N*N;++i){ if(i<10)cout<<grad_obj[i]<<endl;}
obj = value(a1s_obj);
//DCO_MODE::global_tape->interpret_adjoint();
DCO_TAPE_TYPE::remove(DCO_MODE::global_tape);

delete [] a1s_U;

}
int main(int argc, char* argv[]){
  int i;
  double *const U = new double [Nsize*Nsize];
  double *const O = new double [Nsize*Nsize];
  double *const c = new double [nc];
  double local_obj;
  double *const local_U = new double [Nsize*Nsize];

  double *const local_c = new double [nc];
  double endT=5,gradot;
  for (i=0;i<nc;i++)c[i] = 0.5*(i+1);
  double st = clock();
  simulation(c,nc,Nsize,endT,O);
  //for(i=0;i<Nsize*Nsize;i++)cout<<O[i]<<endl;
  std::fill_n(c,Nsize*Nsize,1.0);
  //steepest descent
  double epsilon = 1e-5;
  gradot=epsilon;
  double *gradient_obj = new double [nc];
  std::fill_n(gradient_obj,nc,epsilon);
  //for(i=0;i<10;i++)cout<<gradient_obj[i]<<endl;
  int it = 0;
  while (gradot>=epsilon){
    double obj;
  	derivative_a1s(c,nc,Nsize,endT,U,O,obj,gradient_obj);
  	//cout<<gradient_obj[3]<<endl;
    //cout<<obj<<endl;
    gradot=sqrt(dotpr(gradient_obj,gradient_obj,nc));
    //cout<<"Gradient_obj="<<sqrt(dotpr(gradient_obj,gradient_obj,Nsize*Nsize))<<"\t"<<obj<<endl;
    //for(i=0;i<10;i++)cout<<gradient_obj[i]<<endl;

    double alfa = 0.01;
    local_obj=obj;
    while (local_obj>=obj){
      for(i=0;i<Nsize*Nsize;++i) {
        local_U[i] = U[i];}
      for(i=0;i<nc;i++){
        local_c[i] = c[i]-alfa*gradient_obj[i];
      }
      local_obj = objective_function(local_c,nc,Nsize,endT,local_U,O);
      cout<<local_c[2]<<endl;
      //cout<<"local_obj="<<local_obj<<endl;
      alfa/=2.0;
    }
      for(i=0;i<Nsize*Nsize;i++){c[i] = local_c[i];}
      //cout<<"Thermal diffusivity = "<<c[10]<<endl;
  }
  double et = clock();
  double t = (et-st)/CLOCKS_PER_SEC;
  cout<<"Time="<<t<<endl;
  for(i=0;i<nc;i++)
	{
	cout<<"c["<<i<<"]="<<c[i]<<endl;
	}	

}
