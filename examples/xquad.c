#include <stdio.h>
#include <math.h>
#include "nl.h"

/*
  Использование функций из модуля @quad.h@
  Вычисление определенного интеграла
  $\int_0^1 \sqrt{x} dx$
  адаптивными методами Симпсона, Лобатто-Кронрода, Гаусса--Кронрода
*/

double fun(double x)
{
  return sqrt(x);
}

void print_results(double q, int rc, int nfunevals)
{
   switch (rc)
   {
      case 0:
         printf("Успешное завершение\n");
         break;
      case -1:
         printf("Очередной отрезок не содержит машинных чисел\n");
         printf("Требуемая точность не может быть достигнута\n");
         break;
      case -2:
         printf("Количество вычислений значений подынтегральной функции\n");
         printf("превысило максимально допустимое\n");
         break;
   }

   printf("Значение интеграла: %14.12f\n", q);
   printf("Количество обращений к функции: %d\n", nfunevals);
}


int main()
{
   double q;
   int rc, nfunevals;

   printf ("\nМетод Симпсона:\n");
   q = quad_simpson(&fun, 0, 1, 1e-6, 10000, &rc, &nfunevals);
   print_results(q, rc, nfunevals);

   printf ("\nМетод Лобатто-Кронрода 4,7:\n");
   q = quad_lk47(&fun, 0, 1, 1e-6, 10000, &rc, &nfunevals);
   print_results(q, rc, nfunevals);

   printf ("\nМетод Гаусса-Кронрода 7,15:\n");
   q = quad_gk715(&fun, 0, 1, 1e-6, 10000, &rc, &nfunevals);
   print_results(q, rc, nfunevals);

   return 0;
}

