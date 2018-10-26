#include <stdio.h>
#include <math.h>

#define M 10
void PrintMatrix(double a[M][M], int n);

int main(void){
  double a[M][M],v[M][M];
  double eps=0.0001,div,r,t,s,c,apj,aqj,aip,aiq,vip,viq;
  int i,j,n,iter,count,iterMax=1000000,p,q;
  
  scanf("%d",&n);
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++) scanf("%lf",&a[i][j]);
  }
  PrintMatrix(a,n);
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++) v[i][j]=0.;    
    v[i][i]=1.;
  }
  
  for(iter=1;iter<=iterMax;iter++){
    count=0;
    for(p=1;p<=n-1;p++){
      for(q=p+1;q<=n;q++){
	if(fabs(a[p][q])<eps) continue;
	count++;
	div=a[p][p]-a[q][q];
	if (div != 0.0){
	  r=2.0*a[p][q]/div;
	  t=0.5*atan(r);
	} else {
	  t=0.78539818;
	}
	s=sin(t);
	c=cos(t);
	for(j=1;j<=n;j++){
	  apj=a[p][j];
	  aqj=a[q][j];
	  a[p][j]=apj*c+aqj*s;
	  a[q][j]=-apj*s+aqj*c;
	}
	for(i=1;i<=n;i++){
	  aip=a[i][p];
	  aiq=a[i][q];
	  a[i][p]=aip*c+aiq*s;
	  a[i][q]=-aip*s+aiq*c;
	  vip=v[i][p];
	  viq=v[i][q];
	  v[i][p]=vip*c+viq*s;
	  v[i][q]=-vip*s+viq*c;
	}
	printf("p,q=%3d,%3d\n",p,q);
	PrintMatrix(a,n);
      }
    }
    if (count==0) break;
  }
  printf("Eigen values:\n");
  for(i=1;i<=n;i++) printf("%6.2f",a[i][i]);
  printf("\nEigen vectors:\n");
  PrintMatrix(v,n);
  
  return 0;
}

void PrintMatrix(double a[M][M], int n){
  int i,j;
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++) printf("%6.2f",a[i][j]);
    printf("\n");
  }
  printf("\n");
}

