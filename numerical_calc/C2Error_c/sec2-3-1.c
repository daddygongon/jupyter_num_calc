//C23.c
#include <stdio.h>

int main(void){
  float a=1.23456789;
  printf("  a= %17.10f\n",a);

  float b=100.0;
  float c=a+b;
  printf("%20.10f %20.10f %20.10f\n",a,b,c);
  
  double x=(float)1.23456789;
  double y=(double)100;
  double z=x+y;
    printf("%20.10f %20.10f %20.10f\n",x,y,z);
    printf("%20.12e %20.12e %20.12e\n",x,y,z);

  x=(double)1.23456789;
  y=(double)100;
  z=x+y;
    printf("%20.12e %20.12e %20.12e\n",x,y,z);
    //    printf("%20.16e %20.16e %20.16e\n",x,y,z);
  return 0;
}
