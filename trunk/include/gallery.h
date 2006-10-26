#ifndef __GALLERY_H__
#define __GALLERY_H__

/**
  \file
  Различные матрицы.

  Файл содержит функции, генерирующие часто встречающиеся в приложениях матрицы.
*/

/**
  \example xlaplace.c
*/

/**	
  Случайная разреженная симметричная матрица S.
    - Вход:
      - \f$n\f$ - размер матрицы
      - \f$nz\f$ - число ненулевых элементов
      -	\f$M\f$ - вспомогательный массив длины \f$n\f$
      -	\f$IS\f$ - массив длины \f$n + 1\f$
      -	\f$JS\f$ - массив длины \f$nz\f$
    - Выход:
      -	\f$IS\f$, \f$JS\f$ - разреженный портрет матрицы в формате RR(U)U
*/
extern size_t gall_rand_sym(
	size_t n,
	size_t nz,
	size_t *IS,
	size_t *JS,
	size_t *M);

/**
  Матрица сеточной аппроксимации уравнения Пуассона.

  - Вход:  
     - \f$m\f$, \f$n\f$  - размер сетки
  - Выход: 
     - \f$IA\f$, \f$JA\f$, \f$AN\f$ - разреженное представление матрицы в RR(O)U-формате
*/
extern void gall_laplace(size_t m, size_t n, size_t *IA, size_t *JA, double *AN);

/**
  Размеры представления матрицы сеточной аппроксимации уравнения Пуассона.

  - Вход:  
     - \f$m\f$, \f$n\f$  - размер сетки
  - Выход: 
     - \f$N\f$ - порядок матрицы
     - \f$nz\f$ - число ненулевых элементов в RR(O)U-формате
*/
extern void gall_laplace_size(size_t m, size_t n, size_t *N,	size_t *nz);

/**
  Матрица сеточной аппроксимации уравнения Пуассона.

  - Вход:  
     - \f$m\f$, \f$n\f$  - размер сетки
  - Выход: 
     - \f$IA\f$, \f$JA\f$, \f$AN\f$, \f$AD\f$ - разреженное представление матрицы в RR(U)U-формате
*/
extern void gall_laplace_sym(size_t m, size_t n, size_t *IA, size_t *JA, double *AN, double *AD);

/**
  Размеры представления матрицы сеточной аппроксимации уравнения Пуассона.

  - Вход:  
     - \f$m\f$, \f$n\f$  - размер сетки
  - Выход: 
     - \f$N\f$ - порядок матрицы
     - \f$nz\f$ - число ненулевых элементов в RR(U)U-формате
*/
extern void gall_laplace_sym_size(size_t m, size_t n, size_t *N,	size_t *nz);
#endif
