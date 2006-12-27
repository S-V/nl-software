/*
  Пример использования функций из модуля @.h@ 
  Корень уравнения
  $x^3 - 2x - 5 = 0$
*/

#include <stdio.h>
#include <math.h>
#include "nl.h"

/**
  Решение уравнения с одним неизвестным.
  Используется комбинированный метод деления пополам (бисекции),
  секущих и обратной квадратичной интерполяции.

  - Вход:
	- @func@ - указатель на функцию вычисляющую правую часть уравнения \f$f(x)=0\f$
	- $a$ - левый конец исходного интервала
	- $b$ - правый конец исходного интервала
	- $tol$ - желаемая длина интервала неопределенности конечного результата

  - Выход:
	- Функция возвращает найденное решение.

  Без проверки предполагается, что на концах отрезка \f$[a, b]\f$
  функция \f$f(x)\f$ имеет разные знаки. Корень уравнения ищется в пределах допуска
  на ошибку \f$4\cdot \varepsilon_{\rm M} \cdot |x| + tol\f$,
  где \f$\varepsilon_{\rm M}\f$ - машинная точность.

  Функция является переводом фортрановской программы zeroin из книги [FMM].
*/

double fzero(double (*func)(double), double a, double b, double tol) 
{

  double eps, tol1, c, d, e, fa, fb, fc, p, q, r, s, xm;

  /* Вычисление машинного эпсилон */

  eps = 1.0;

  do
  {
    eps = eps/2;
    tol1 = 1.0 + eps;
  }
  while (tol1 > 1.0);

  /* Начальные приготовления */

  fa = (*func)(a);
  fb = (*func)(b);
  c = a;
  fc = fa;
  d = b - a;
  e = d;

  /* Основной цикл */

  while(1)
  {

    if (fabs(fc) < fabs(fb))
    {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

     /* Проверка сходимости */

     tol1 = 2.0*eps*fabs(b) + 0.5*tol;
     xm = .5*(c - b);
     if (fabs(xm) <= tol1 || fb == 0.0) 
  	break;
   
     if (fabs(e) < tol1 || fabs(fa) < fabs(fb))
     {
          /* Бисекция */

          d = xm;
          e = d;
     }
     else 
     {
        if (a != c)
        {
          /* Обратная квадратичная интерполяция */

          q = fa/fc;
          r = fb/fc;
          s = fb/fa;
          p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
          q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        }
        else
        {
          /* Линейная интерполяция (метод секущих) */

          s = fb/fa;
          p = 2.0*xm*s;
          q = 1.0 - s;
        }
        
        /* Выбрать знаки */

        if (p > 0.0) 
          q = -q;
        p = fabs(p);

        /* Приемлема ли интерполяция? */

        if (2*p < 3*xm*q - fabs(tol1*q) && p < fabs(0.5*e*q))
        {
          e = d;
          d = p/q;
        }
        else
        {
          /* Все-таки бисекция */

          d = xm;
          e = d;
        }
     }

     /* Завершение итерации */

     a = b;
     fa = fb;

     if (fabs(d) > tol1) 
       b = b + d;
     else
       if (xm >= 0)
         b = b + tol1;
       else
         b = b - tol1;

     fb = (*func)(b);

     if (fb*(fc/fabs(fc)) > 0.0)
     {
       c = a;
       fc = fa;
       d = b - a;
       e = d;
     }

  }
    
  return b;

}

double func(double x)
{
  return x*(x*x - 2) - 5;
}


int main()
{
  double x;

  x = fzero(&func, 2, 3, 1e-6); 

  printf("Корень:  %f  \n", x);

  return 0;
}