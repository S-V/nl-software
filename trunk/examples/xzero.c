#include <stdio.h>
#include <math.h>
#include "nl.h"

/*
  Пример использования функций из модуля @roots.h@ 
  Корень уравнения
  $x^3 - 2x - 5 = 0$
*/

double fun(double x)
{
  return x*(x*x - 2) - 5;
}


int main()
{
  double x;

  x = roots_zero(&fun, 2, 3, 1e-6); 

  printf("Корень:  %f  \n", x);
  printf("Значение функции:  %e  \n", (*fun)(x));

  return 0;
}