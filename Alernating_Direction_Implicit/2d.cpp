#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
using namespace std;
int N;
#define N 100
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
int main(int argc, char* argv[])
{
  int i,j;
  //Initialization of a 2D Array
  double **U_old= new double*[N];
  for(int i = 0; i < N; i++){
		U_old[i] = new double[N];}
  double **U_new= new double*[N];
  for(int i = 0; i < N; i++){
  	U_new[i] = new double[N];}
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
  double k=20;//Thermal Conductivity
  double dt=0.005;//(Half) Time step
  double dx=1.0/N-1;//stepping in x direction
  double dy=1.0/N-1;//stepping in y direction
  double DT;//Total time
  /* R_x and R_y */
  double Rx=k*dt/ (dx * dx);
  double Ry=k*dt/ (dy * dy);
  double a[N],b[N],c[N],d[N],m;

  while(DT<10){//end time
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

  //for printing a 2D Array
    for (i=0;i<N;i++){
      for (j=0;j<N;j++){
        cout<<U_old[i][j];cout<<"\t";}
        cout<<endl;}
  delete[] U_old,U_new;
}
