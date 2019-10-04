#include <stdio.h>
#include <math.h>
#define EPSILON 0.01
int main(void){
  float x=77777,y=7,y1,z,z1;
  
  y1=1/y;
  z=x/y;
  z1=x*y1;
  printf("%10.2f %10.2f\n", z,z1);
  // if (z!=z1){
  if (fabs(z-z1)>EPSILON){
    printf("z is not equal to z1.\n");
  } else {
    printf("z is equal to z1.\n");
  }
  printf("%20.15f %20.15f\n", z,z1);
}
