#include "func.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
using namespace std;
int N;
#define N 10

int main(int argc, char* argv[])
{
  int i,j,k;
  double t1,t2;
  //Initialization of a 2D Array
  double **U_old= new double*[N];
  for(int i = 0; i < N; i++){
		U_old[i] = new double[N];}
  double **U_new= new double*[N];
  for(int i = 0; i < N; i++){
  	U_new[i] = new double[N];}
  double **C= new double*[N];
  for(int i = 0; i < N; i++){
  	C[i] = new double[N];}
  /*Initial Condition*/
  assign(U_old,N,N,100);
  assign(C,N,N,0.01);

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
  boundary(U_old,N,N,0,1);
  boundary(U_old,N,N,0,2);
  boundary(U_old,N,N,0,4);
  boundary(U_old,N,N,0,3);
  //equate(U_new,U_old,N,N);
  double dt=0.005;//(Half) Time step
  double dx=1.0/N-1;//stepping in x direction
  double dy=1.0/N-1;//stepping in y direction
  double DT;//Total time
  /* R_x and R_y */
  double Rx=dt/ (dx * dx);
  double Ry=dt/ (dy * dy);
  //converting to AX=B form
  double *X,*B,*c;
  X = new double[N*N];
  B = new double[N*N];
  c = new double[N*N];
  int point,count;
  count=-1;
  for(i=0;i<N*N;i++){
    point=i%(N);
    if(point==0){count=count+1;}
    B[i]=U_old[count][point];
    c[i]=C[count][point];
    X[i]=0.0;
  }
  double *res,*Ax,*ADk,*d,alfa,beta,rtr,rtr_temp,tol;
  tol=1e-9;
  res = new double[N*N];
  ADk = new double[N*N];
  Ax =new double[N*N];
  d = new double[N*N];
  t1=clock();
  while(DT<10){//end time
  //Conjugate gradient
  multiply(X,Ax,c,Rx,Ry,N);
  for(i=0;i<N*N;i++){
    res[i]=d[i]=B[i];}
    rtr=dotpr(res,res,N*N);
  for(j=0;j<N*N;j++){
    multiply(d,ADk,c,Rx,Ry,N);
    alfa=rtr/dotpr(ADk,d,N*N);
    for(i=0;i<N*N;i++){
      res[i]=res[i] - alfa*ADk[i];
      X[i]=X[i] + alfa*d[i];
    }
    rtr_temp=rtr;
    rtr=dotpr(res,res,N*N);
    if(sqrt(rtr)<tol){break;}
    beta=rtr/rtr_temp;
    for(i=0;i<N*N;i++){
    d[i]=res[i] + beta*d[i];}
  }
  for(i=0;i<N*N;i++){
      B[i]=X[i];
      X[i]=0.0;
  }
  DT=DT+dt*2;
  }
t2=clock();
  count=-1;
  for(i=0;i<N*N;i++){
    point=i%(N);
    if(point==0){count=count+1;}
    U_old[count][point]=B[i];
  }

  //for printing a 2D Array
    for (i=0;i<N;i++){
      for (j=0;j<N;j++){
        cout<<U_old[i][j];cout<<"\t";}
        cout<<endl;}
        cout<<"time="<<(t2-t1)/CLOCKS_PER_SEC<<endl;

  delete[] U_old;
  delete[] U_new;
}
