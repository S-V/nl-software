#ifndef __INTERP_H__
#define __INTERP_H__

/**
  \file
  ������������ ������ �������-��������������� ���������

  � ������ ������ ���������� �������, �������� �������-����������
  ������������. ����� 
  \f$x_0 < x_1 < \dots < x_{n-1}\f$
  \f$y_0 < y_1 < \dots < y_{n-1}\f$.
  ���������� �������-���������� �������: 
  \f$P(x)\f$, ������ ����������� �������� 
  \f$p_j(x) = a_{j0}x^3 + a_{j1}x^2 + a_{j2}x + a{j1}\f$ ��� \f$x_j \le x < x_{j+1}\f$
  \f$(j=0,1,\dots,n - 2)\f$. 
  
  ������� \f$p_j(x)\f$ ������ �������� � ����
  \f$p_j(x)= y_j + d_j(x-x_j) +c_j(x-x_j)^2+b_j(x-x_j)\f$.
  ������� ��������, ��� \f$y_j = p_j(x_j)\f$, \f$d_j = p'_j(x_j)\f$,
  \f$c_j = p''_j(x_j)/2\f$, \f$b_j = p'''_j(x_j)/6\f$.

  \f$P(x)\f$ ���������� </i>���������� ��������� �������������</i>,
  ���� \f$P(x)\f$ ����������, ����� ����������� ����������� � \f$P(x_j)= y_j\f$ \f$(j=0,1,\dots,n - 1)\f$,
  ����� �������, \f$p_j(x_j)=y_j\f$, \f$p_j(x_{j+1})=y_{j+1}\f$ \f$(j=0,1,\dots,n - 1)\f$, 
  \f$p'_j(x_j)=p'_{j+1}(x_j)\f$ \f$(j=0,1,\dots,n - 2)\f$. 
  ���������� ������� ����������� �� �����������.

  \f$P(x)\f$ ���������� </i>���������� ��������</i>,
  ���� \f$P(x)\f$ ����������, ����� ����������� ������ � ������ ����������� 
  � \f$P(x_j)= y_j\f$ \f$(j=0,1,\dots,n - 1)\f$,
  ����� �������, \f$p_j(x_j)=y_j\f$, \f$p_j(x_{j+1})=y_{j+1}\f$ \f$(j=0,1,\dots,n - 1)\f$, 
  \f$p'_j(x_j)=p'_{j+1}(x_j)\f$ \f$(j=0,1,\dots,n - 2)\f$,
  \f$p''_j(x_j)=p''_{j+1}(x_j)\f$ \f$(j=0,1,\dots,n - 2)\f$. 
  ���������� ������ �� �����������.
*/

/**
  \example xinterp.c
*/

/**
  \example xiquad.c
*/

/**
  ���������� ������� �����������

  ������� ������ ���������� ������� �����������. ����������� �������������
  �������������� ��������: �� ��� ��������, ��� ������ ���������
  ����������� ����� ���������; ���� \f$x_j\f$, ���������� ����������
  ����������� (����������), �������� ���������� ����������� (����������) �
  ��� ������������.

  - ����:
      - \f$x\f$ (������ ����� \f$n\f$) - �������� �����
      - \f$y\f$ (������ ����� \f$n\f$) - �������� �����
      - \f$n\f$ ���������� �����

  - �����:
      - \f$b\f$ (������ ����� \f$n - 1\f$) 
        \f$c\f$ (������ ����� \f$n - 1\f$) 
        \f$d\f$ (������ ����� \f$n\f$) - ������������ ���������� ���������.
        \f$d[n-1]\f$ ����� �������� ������������ � ����� \f$x[n-1]\f$

  ������ ������ ������������� ��������: 
  \f$n \ge 3\f$, \f$x[j] < x[j + 1]\f$ \f$(j=0,1,\dots,n - 2)\f$.

  ������������ \f$25n + O(1)\f$
*/

extern void interp_pchip(
   double *x, 
   double *y, 
   size_t n, 
   double *b, 
   double *c, 
   double *d);

/**
  ���������� ������

  ������� ������ ���������� ������. ������ �������������
  ��������������� �������: � ������ \f$x_1\f$ � \f$x_{n-2}\f$ 
  �� ����� ����������� <i>������</i> �����������
  (������� ������� ����� - ������ � ��� �������������� �������� �����������).

  - ����:
      - \f$x\f$ (������ ����� \f$n\f$) - �������� �����
      - \f$y\f$ (������ ����� \f$n\f$) - �������� �����
      - \f$n\f$ ���������� �����

  - �����:
      - \f$b\f$ (������ ����� \f$n - 1\f$) 
        \f$c\f$ (������ ����� \f$n - 1\f$) 
        \f$d\f$ (������ ����� \f$n\f$) - ������������ ���������� ���������.
        \f$d[n-1]\f$ ����� �������� ������� � ����� \f$x[n-1]\f$

  - ������� ������:
      - \f$work\f$ - ������ ����� \f$2n\f$

  ������ ������ ������������� ��������: 
  \f$n \ge 4\f$, \f$x[j] < x[j + 1]\f$ \f$(j=0,1,\dots,n - 2)\f$.

  ������������ \f$32n + O(1)\f$
*/

extern void interp_spline(
   double *x, 
   double *y, 
   size_t n, 
   double *b, 
   double *c, 
   double *d, 
   double *work);


/**
  ���������� �������� ����������� ������������

  ������� ��������� �������� �������� ����������� ������������
  ��� ����������� ������� � �������� �����.

  - ����:
      - \f$x\f$ (������ ����� \f$n\f$) - �������� ����� (����� ������������)
      - \f$y\f$ (������ ����� \f$n\f$) - �������� ����� (����� ������������)
      - \f$b\f$ (������ ����� \f$n - 1\f$) 
        \f$c\f$ (������ ����� \f$n - 1\f$) 
        \f$d\f$ (������ ����� \f$n\f$) - ������������ ���������� ���������.
        (����� ���� ��������� ��������� #interp_pchip ��� #interp_spline)
      - \f$n\f$ ���������� ����� ������������
      - \f$xx\f$ - ������������ �����

  - �����:
      - ������� ���������� �������� ������������ � ����� \f$xx\f$

  ������ ������ ������������� ��������: 
  \f$n \ge 2\f$, \f$x[j] < x[j + 1]\f$ \f$(j=0,1,\dots,n - 2)\f$.
*/

double interp_eval(
   double *x, 
   double *y, 
   double *b, 
   double *c, 
   double *d, 
   size_t n, 
   double xx);

/**
  ���������� �������� ������������� ��������� �� ����������� ������������

  ������� ������� ������������ �������� �� ����������� ������������
  \f$\int_{x_1}^{x_2} P(x) dx\f$. ��������� �������� �������� ������������
  � ������������� ��������� �� �������, �������� �������� ����������
  \f$x_j\f$, \f$y_j\f$

  - ����:
      - \f$x\f$ (������ ����� \f$n\f$) - �������� ����� (����� ������������)
      - \f$y\f$ (������ ����� \f$n\f$) - �������� ����� (����� ������������)
      - \f$b\f$ (������ ����� \f$n - 1\f$) 
        \f$c\f$ (������ ����� \f$n - 1\f$) 
        \f$d\f$ (������ ����� \f$n\f$) - ������������ ���������� ���������.
        (����� ���� ��������� ��������� #interp_pchip ��� #interp_spline)
      - \f$n\f$ ���������� ����� ������������
      - \f$x1\f$ - ������ ������ ��������������
      - \f$x2\f$ - ������� ������ ��������������

  - �����:
      - ������� ���������� �������� ������������� ���������

  ����� \f$x1\f$, \f$x2\f$ ����� ������ ��� ������� \f$x[0], x[n - 1]\f$.
  � ���� ������ ��������� � ���������� ��������� ��� � �����������
  ��������� �� �������� �������� ������� ����� ������ �����������������.

  ������ ������ ������������� ��������: 
  \f$n \ge 2\f$, \f$x[j] < x[j + 1]\f$ \f$(j=0,1,\dots,n - 2)\f$.
*/

double interp_quad(
   double *x, 
   double *y, 
   double *b, 
   double *c, 
   double *d, 
   size_t n, 
   double x1, 
   double x2);

#endif
