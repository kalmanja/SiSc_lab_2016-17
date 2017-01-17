void assign(double **a,int row,int col,double val);
void boundary(double **a,int row,int col,double val,int bound);
double dotpr(double *u,double *v,int N);
void equate_vector(double *u,double* v,int N);
void divide(int *chunk,int *lcount,int N,int size);
void multiply(double *u, double *ans,double *c, double rx, double ry, int N,int start,int len,int rank);
void distribute(double *X,int *lcount,int *chunk,int len,int N,int size,int rank);
