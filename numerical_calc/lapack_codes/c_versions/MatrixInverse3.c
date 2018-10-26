#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vecLib/vecLib.h>

//#undef LAPACK
#define LAPACK
#undef PRINT

void printMatrix(double *a, double *b, long n);
int MatrixInverse(double *a, double *b, long n);

int main(void){
  clock_t start, end;
  int i,j;
  long n;
  double *a,*b;

  scanf("%ld",&n);
  printf("%d ",n);

  a=(double *)malloc(n*n*sizeof(double));
  b=(double *)malloc(n*sizeof(double));

  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      a[i*n+j]= 2*(double) random() / RAND_MAX - 1.0;
    }
  }

  for (i=0;i<n;i++){
    b[i]= 2*(double) random() / RAND_MAX - 1.0;
  }
  printMatrix(a,b,n);
#ifndef LAPACK
  printf("BOB");
  start = clock();
  MatrixInverse(a,b,n);
  end = clock();
  printf("%5d [dim] %10.4f [sec] #BOB\n", n,(double)(end-start)/CLOCKS_PER_SEC);
#else
  printf("LAPACK");
  long nrhs=1, lda,ldb, info;
  long *ipiv;
  lda=ldb=n;
  ipiv=(long *)malloc(n*sizeof(long));
  double tmp;
  for (i=0;i<n;i++){
    for (j=i+1;j<n;j++){
      tmp=a[i*n+j];
      a[i*n+j]=a[j*n+i];
      a[j*n+i]=tmp;
    }
  }
  printMatrix(a,b,n);
  start = clock();
  dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  end = clock();
  printf("%5d [dim] %10.4f [sec] #LAPACK\n", n, (double)(end-start)/CLOCKS_PER_SEC);
  free(ipiv);
#endif

  printMatrix(a,b,n);

  free(a);
  free(b);
  return 0;
}

int MatrixInverse(double *a, double *b, long n){
  double *x;
  double pvt=0.00005,am;
  int i,j,k;

  x=(double *)malloc(n*sizeof(double));

  for(i=0;i<n-1;i++){
    if(fabs(a[i*n+i])<pvt){
      printf("Pivot %3d=%10.5f is too small.\n",i,a[i*n+i]);
      return 1;
    }
    for(j=i+1;j<n;j++){
      am=a[j*n+i]/a[i*n+i];
      for(k=0;k<n;k++) a[j*n+k]-=am*a[i*n+k];
      b[j]-=am*b[i];
    }
  }
  //Backward substitution
  for(j=n-1;j>=0;j--){
    x[j]=b[j];
    for(k=j+1;k<n;k++){
      x[j]-=a[j*n+k]*x[k];
    }
    b[j]=x[j]/=a[j*n+j];
  }
  free(x);
  return 0;
}

void printMatrix(double *a, double *b, long n){
  int i,j;
#ifdef PRINT
  printf("\n");
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      printf("%10.5f",a[i*n+j]);
    }
    printf(":%10.5f",b[i]);
    printf("\n");
  }
  printf("\n");
#endif
  return;
}
