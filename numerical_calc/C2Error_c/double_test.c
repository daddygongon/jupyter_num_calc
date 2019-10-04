#include <stdio.h>

int main(void){
  float a,b,d;
  double x,y,z;

  a=1.23456789;
  printf("   a=%17.10f\n",a);

  x=(float)1.23456789;
  y=(double)100;
  z=x+y;
  printf("%20.12e, %20.12e, %20.12e\n",x,y,z);

  x=(double)1.23456789;
  y=(double)100;
  z=x+y;
  printf("%20.12e, %20.12e, %20.12e\n",x,y,z);
}
