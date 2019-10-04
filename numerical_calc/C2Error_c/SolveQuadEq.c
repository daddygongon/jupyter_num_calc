#include <stdio.h>
#include <math.h>

int main(void){
  printf("Hello\n");
  //  float a,b,c,x1,x2;
  double a,b,c,x1,x2;
  a=1.0;
  b=-10000000.0;
  c=1.0;

  x1=(-b-sqrt(b*b-4.0*a*c))/(2.0*a);
  x2=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
  printf("%20.10f, %20.10f\n",x1,x2);
  printf("%20.10e, %20.10e\n",x1,x2);
  if (b>0){
    x1=(-b-sqrt(b*b-4.0*a*c))/(2.0*a);
  } else {
    x1=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
  }
    x2=c/(a*x1);
  printf("%20.10f, %20.10f\n",x1,x2);
  printf("%20.10e, %20.10e\n",x1,x2);
}

/*in case:   float a,b,c,x1,x2;
[bob-no-MacBook-Pro:NumericalRecipe/NumRecipe12/C2Error_c] bob% !.
./a.out
Hello
  -1.000000000000000e+07, 9.420700371265411e-03
  -1.000000000000000e+07, -1.000000011686097e-07
*/

/* double
0.0000000997,   9999999.9999998994
  9999999.9999998994,         0.0000001000
*/
