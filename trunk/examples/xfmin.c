/*
  Пример использования функций из модуля @.h@ 
  Корень уравнения
  $x^3 - 2x - 5 = 0$
*/

#include <stdio.h>
#include <math.h>
#include "nl.h"

/**
  Минимизация функции одного переменного.
  Используется комбинированный метод золотого сечения
  и квадратичной интерполяции.

  - Вход:
	- \f$func\f$ - указатель на функцию
	- \f$a\f$ - левый конец исходного интервала
	- \f$b\f$ - правый конец исходного интервала
	- \f$tol\f$ - желаемая длина интервала неопределенности конечного результата

  - Выход:
	- Функция возвращает найденное решение.

  Если \f$f\f$ имеет непрерывную вторую производную, положительную в точке
  минимума (не совпадающей ни с \f$a\f$, ни с \f$b\f$), то сходимость
  сверхлинейная порядка приблизительно \f$1.324\dots\f$
  В любом случае сходимость не может быть хуже, чем для метода золотого сечения.

  Функция является переводом фортрановской программы fmin из книги [FMM].
*/

double fmin(double (*func)(double), double a, double b, double tol) 
{

  double eps, tol1, tol2, c, d, e, x, u, v, w, fx, fu, fv, fw, xm, p, q, r;


  c = 0.5*(3 - sqrt(5.0));


  /* eps равно корню квадратному из машинного эпсилон */

  eps = 1.0;

  do
  {
    eps = eps/2;
    tol1 = 1.0 + eps;
  }
  while (tol1 > 1.0);

  eps = sqrt(eps);

  v = a + c*(b - a);
  w = v;
  x = v;
  e = 0;
  fx = (*func)(x);
  fv = fx;
  fw = fx;

  /* Основной цикл */

  while (1)
  {
      xm = 0.5*(a + b);
      tol1 = eps*fabs(x) + tol/3;
      tol2 = 2*tol1;

      if (fabs(x - xm) <= tol2 - 0.5*(b - a)) 
        break;

      if (fabs(e) <= tol1)
      {
         /* Шаг метода золотого сечения */

         if (x >= xm) 
           e = a - x;
         else 
           e = b - x;

         d = c*e;
      }
      else
      {
        /* Построение параболической интерполянты */

        r = (x - w)*(fx - fv);
        q = (x - v)*(fx - fw);
        p = (x - v)*q - (x - w)*r;
        q = 2*(q - r);

        if (q > 0) 
          p = -p;

        q = fabs(q);
        r = e;
        e = d;

        /* Приемлема ли парабола? */

        if (fabs(p) >= fabs(0.5*q*r) || p <= q*(a - x) || p >= q*(b - x)) 
        {
           /* Все-таки метод золотого сечения */

           if (x >= xm) 
             e = a - x;
           else 
             e = b - x;

           d = c*e;
        }
        else
        {
          /* Шаг параболической интерполяции */

          d = p/q;
          u = x + d;

          /* f не должна вычисляться слишком близко к a или b */

          if ((u - a) < tol2)
          {
            if (xm >= x)
              d = tol1;
            else
              d = -tol1;
          }
          else if ((b - u) < tol2)
          { 
            if (xm >= x)
              d = tol1;
            else
              d = -tol1;
          }
        }
      }

      /* f не должна вычисляться слишком близко к x */

      if (fabs(d) >= tol1) 
         u = x + d;
      else
      {
         if (d >= 0) 
           u = x + tol1;
         else
           u = x - tol1;
      }

      fu = (*func)(u);

      /* Пересчет a, b, v, w, x */

      if (fu <= fx) 
      {
        if (u >= x) 
          a = x;
        else 
          b = x;

        v = w;
        fv = fw;
        w = x;
        fw = fx;
        x = u;
        fx = fu;
      }
      else
      {
        if (u < x) 
          a = u;
        else
          b = u;

        if (fu <= fw || w == x)
        {
          v = w;  
          fv = fw;
          w = u;  
          fw = fu;
        }
        else if (fu <= fv || v == x || v == w) 
        {
          v = u;  
          fv = fu;
        }
      }
   }

   return x;

}

double func(double x)
{
  return x*(x*x - 2) - 5;
}


int main()
{
  double x, fx;

  x = fmin(&func, 0, 1, 1e-6);
  fx = (*func)(x); 

  printf("Точка минимума:    %f\n", x);
  printf("Значение минимума: %f\n", fx);


  return 0;
}