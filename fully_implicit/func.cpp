#include "func.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
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
/*Function for equating two 2D Arrays*/
void equate(double **a,double **b,int row,int col){
  int i,j;
  for(i=0;i<row;i++){
    for(j=0;j<col;j++){
      a[i][j]=b[i][j];}}
}
void multiply(double *u, double *ans,double *c, double rx, double ry, int N){
int i;
for(i=0;i<N*N;i++){
  if(i<N+1 || (i+1)%N<2|| i>(N*N - N-2)){ans[i]=u[i]*(1);}
  else{ans[i]=u[i]*(1 + 2*rx*c[i] + 2*ry*c[i]) + u[i+1]*(-rx*c[i+1]) + u[i-1]*(-rx*c[i-1]) + u[i-N]*(-ry*c[i-N]) + u[i+N]*(-ry*c[i+N]);}
}}
double dotpr(double *u,double *v,int N){
  int i;
  double x;
  x=0;
  for(i=0;i<N;i++){
    x=x+(u[i]*v[i]);
  }
  return x;
}
