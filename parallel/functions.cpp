#include "functions.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <mpi.h>
using namespace std;
/*Function for assigning values to 2D Array (Preferably for Initial Conditions)*/
void assign(double **a,int row,int col,double val)
{
    int i,j;
    for (i=0;i<row;i++){
      for (j=0;j<col;j++){
        a[i][j]=val;}}
}
/*Fuction for assigning boundary Conditions*/
void boundary(double **a,int row,int col,double val,int bound){
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
double dotpr(double *u,double *v,int N){
  int i;
  double x;
  x=0;
  for(i=0;i<N;i++){
    x=x+(u[i]*v[i]);
  }
  return x;
}
void equate_vector(double *u,double* v,int N){
  int i;
  for(i=0;i<N;i++){
    u[i]=v[i];
  }
}
void divide(int *chunk,int *lcount,int N,int size){
  int i,temp,rows;
  rows=N;
  lcount[0]=0;
  for(i=0;i<size;i++){
    temp=ceil(1.0*rows/(size-i));
    chunk[i]=N*temp;
    rows=rows-temp;
    if(i>0){lcount[i]=lcount[i-1]+chunk[i-1];}
  }
}
void multiply(double *u, double *ans,double *c, double rx, double ry, int N,int start,int len,int rank){
  int i,temp;
  temp=start;
  for(i=0;i<len;i++){
    if(temp<N+1 || (temp+1)%N<2|| temp>(N*N - N-2)){if(rank>0)ans[i]=u[i+N]*(1);else ans[i]=u[i]*(1);}
    else{
      if(rank>0)ans[i]=u[i+N]*(1 + 2*rx*c[i+N] + 2*ry*c[i+N]) + u[i+1+N]*(-rx*c[i+1+N]) + u[i-1+N]*(-rx*c[i-1+N]) + u[i]*(-ry*c[i]) + u[i+2*N]*(-ry*c[i+2*N]);
      else ans[i]=u[i]*(1 + 2*rx*c[i] + 2*ry*c[i]) + u[i+1]*(-rx*c[i+1]) + u[i-1]*(-rx*c[i-1]) + u[i-N]*(-ry*c[i-N]) + u[i+N]*(-ry*c[i+N]);}
  temp=temp+1;}}
  void distribute(double *X,int *lcount,int *chunk,int len,int N,int size,int rank){
    int i;
    if(rank==0){
    for(i=1;i<size;i++){
      if(i==size-1){MPI_Send(&X[lcount[i]-N],chunk[i]+N, MPI_DOUBLE,i,0, MPI_COMM_WORLD);}
      else{MPI_Send(&X[lcount[i]-N],chunk[i]+2*N, MPI_DOUBLE,i,0, MPI_COMM_WORLD);}
    }}
    if(rank==size-1){
      MPI_Recv(X,len+N, MPI_DOUBLE, 0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
      if(rank>0 && rank<size-1){
      MPI_Recv(X,len+2*N, MPI_DOUBLE, 0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
