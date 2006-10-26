#include <stdio.h>

int main()
{
  double mach_eps;
  
  mach_eps = 1;
  do
    mach_eps /= 2;
  while((mach_eps + 1) > 1);
  mach_eps *= 2;

  printf("Машинное эпсилон = %e\n", mach_eps);

  return 0;
}
