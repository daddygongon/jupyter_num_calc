#include <stdio.h>


int main(void){
  double e=1.0,w;
  //  float e=1.0,w;
  int i;
  w=1.0+e;
  printf("Hello world!!\n");
  while (w>1.0) {
    printf("%-15.10e %-15.10e %-15.10e\n",e,w,w-1.0);
    e=e/2.0;
    w=1.0+e;
  }
}
