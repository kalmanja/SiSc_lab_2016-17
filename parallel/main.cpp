#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include "functions.h"
#include <mpi.h>
using namespace std;

int N;
#define N 100

int main(int argc,char* argv[]){
  int rank, size;//N = mesh size
  int len,start,end;
  double *U,*c,*X,*ans,**U_old,**C,*res,*b,*d,*Adk,rtr,rtr_temp,dAd,alfa,beta,tol, t1,t2;
  int temp,i,j,k,count,point,*chunk=NULL,*lcount=NULL,*rcount=NULL;//collective communication
  double dt=0.05;//(Half) Time step
  double dx=1.0/N-1;//stepping in x direction
  double dy=1.0/N-1;//stepping in y direction
  double DT;//Total time
  tol=1e-9;//tolerence
  /* R_x and R_y */
  double Rx=dt/ (dx * dx);
  double Ry=dt/ (dy * dy);
  MPI_Init (NULL,NULL);      /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */
  if(rank==0){
    U_old= new double*[N];
    for(i = 0; i < N; i++){
      U_old[i] = new double[N];}
    C= new double*[N];
    for(i = 0; i < N; i++){
      C[i] = new double[N];}
    assign(U_old,N,N,100);//initial condition
    assign(C,N,N,0.01);
    boundary(U_old,N,N,0,1);//boundary conditions
    boundary(U_old,N,N,0,2);
    boundary(U_old,N,N,0,4);
    boundary(U_old,N,N,0,3);
    chunk =new int[size];
    lcount =new int[size];
    U= new double[N*N];//converting to AX=B form
    c= new double[N*N];
    X= new double[N*N];
    res= new double[N*N];
    b= new double[N*N];
    count=-1;
    for(i=0;i<N*N;i++){
      point=i%(N);
      if(point==0){count=count+1;}
      b[i]=U_old[count][point];//matrix to vector
      c[i]=C[count][point];//matrix to vector
      X[i]=0.0;//intial guess for CG
    }
    divide(chunk,lcount,N,size);

  }//end of if(rank == 0)

  MPI::COMM_WORLD.Scatter(chunk,1, MPI::INT,&len,1,MPI::INT,0);
  MPI::COMM_WORLD.Scatter(lcount,1, MPI::INT,&start,1,MPI::INT,0);
  if (rank==0){d=new double[N*N];}
  if(rank>0 && rank<size-1){
    b=new double[len+2*N];
    c=new double[len+2*N];
    d=new double[len+2*N];
    res=new double[len];
    X= new double[len];
  }
  if(rank==size-1){
    b=new double[len+N];
    c=new double[len+N];
    d=new double[len+N];
    res=new double[len];
    X= new double[len];
  }
  ans=new double[len];
  Adk=new double[len];
  distribute(c,lcount,chunk,len,N,size,rank);
  if(rank==0)t1=MPI_Wtime();
  while(DT<10){  /*Conjugate Gradient*/
  distribute(b,lcount,chunk,len,N,size,rank);
  if(rank==0){
    for(i=0;i<len;i++){
      res[i]=b[i];}
    equate_vector(d,b,len+N);
  }
  else{
    for(i=N;i<N+len;i++){
      res[i-N]=b[i];
      }
    if(rank<size-1)equate_vector(d,b,len+2*N);
    else{equate_vector(d,b,len+N);}
  }
  rtr=dotpr(res,res,len);
  MPI_Allreduce(MPI_IN_PLACE,&rtr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for (k=0;k<N*N;k++){
    multiply(d,Adk,c,Rx,Ry,N,start,len,rank);
    if(rank==0)dAd=dotpr(d,Adk,len);
    else dAd=dotpr(&d[N],Adk,len);
    MPI_Allreduce(MPI_IN_PLACE,&dAd,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    alfa=rtr/dAd;
    for(i=0;i<len;i++){
      if(rank==0)X[i]=X[i] + alfa*d[i];
      else X[i]=X[i] + alfa*d[i+N];
      res[i]=res[i] -alfa*Adk[i];
    }
    rtr_temp=rtr;
    rtr=dotpr(res,res,len);
    MPI_Allreduce(MPI_IN_PLACE,&rtr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(sqrt(rtr)<tol){break;}
    beta=rtr/rtr_temp;
    for(i=0;i<len;i++){
      if(rank==0)d[i]=res[i] + beta*d[i];
      else d[i+N]=res[i] + beta*d[i+N];
    }
    if(rank==0)MPI_Gatherv(d,len,MPI_DOUBLE,d,chunk,lcount,MPI_DOUBLE,0,MPI_COMM_WORLD);
    else MPI_Gatherv(&d[N],len,MPI_DOUBLE,d,chunk,lcount,MPI_DOUBLE,0,MPI_COMM_WORLD);
    distribute(d,lcount,chunk,len,N,size,rank);
  }//end of k iteration
  MPI_Gatherv(X,len,MPI_DOUBLE,X,chunk,lcount,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if(rank==0)for(i=0;i<N*N;i++){b[i]=X[i];X[i]=0;}
  else{for(i=0;i<len;i++)X[i]=0;}
  DT=DT+dt*2;
  }
  if(rank==0)t2=MPI_Wtime();
  /*if(rank==0){
    count=-1;
    for(i=0;i<N*N;i++){
      point=i%(N);
      if(point==0){count=count+1;}
      U_old[count][point]=b[i];}
    for (i=0;i<N;i++){
      for (j=0;j<N;j++){
        cout<<U_old[i][j];cout<<"\t";}
        cout<<endl;}
  }*/
  if(rank==0)cout<<"time="<<t2-t1<<endl;
  //cout<<rank<<"\t"<<len<<"\t"<<start<<"\t"<<end<<endl;
  MPI_Finalize();
  return 0;

}
