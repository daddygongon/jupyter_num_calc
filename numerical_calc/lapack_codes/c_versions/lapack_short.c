#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vecLib/vecLib.h>

void printMatrix(double *a, double *b, long n);

int main(void){
  clock_t start, end;
  int i,j;
  double *a,*b;
  long n,nrhs=1, lda,ldb, info, *ipiv;

  scanf("%ld",&n);

  a=(double *)malloc(n*n*sizeof(double));
  b=(double *)malloc(n*sizeof(double));
  lda=ldb=n;
  ipiv=(long *)malloc(n*sizeof(long));
  
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      a[j*n+i]= 2*(double) random() / RAND_MAX - 1.0;
    }
  }

  for (i=0;i<n;i++){
    b[i]= 2*(double) random() / RAND_MAX - 1.0;
  }
  printMatrix(a,b,n);

  start = clock();
  dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  end = clock();
  printf("%5d [dim] %10.4f [sec] #LAPACK\n", 
	 n, (double)(end-start)/CLOCKS_PER_SEC);
  printMatrix(a,b,n);

  free(a);
  free(b);
  free(ipiv);

  return 0;
}
