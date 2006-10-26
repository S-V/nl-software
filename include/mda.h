#ifndef __MDA_H
#define __MDA_H

#include "util.h"

/**
  \file
  ����� ����������� ������� ��� ���������������� ������������ ����������� �������.

  ����� ����������� ������� ������������ �� ��������������� �����
  ��������� ������� ����������� ������������ �������. �� ������� ���������������� �����
  � �������� �������� ������������ ������������ ������������ �������, �����, ��� 
  ���������� ��������� ���������� �������, ��� �������, ����� ��������� ����������.

  ��� ������� ������� ��������� \f$Ax=b\f$ � ������������ ������������ ������������ 
  ����������� �������� \f$A\f$ �������:

    - � ������� ������� #mda_init, #mda_convert, #mda_order ����� ���������������� 
      \f$P\f$ ����� � �������� ������� \f$A\f$;
    - � ������� ������� #sp_permute_sym ����������� ������ � ������� � ������� \f$A\f$
      ��������� �������: \f$PAP^{T}\f$;
    - � ������� ������� #sp_chol_symb, #sp_chol_num ����� ���������� ���������
      ������� \f$PAP^{T}\f$;
    - � ������� ������� #nl_dvector_permute ����������� �������� ������� \f$b\f$
      ��������� �������: \f$Pb\f$;
    - � ������� ������� #sp_chol_solve ������ ������� \f$PAP^{T}y=Pb\f$;
    - � ������� ������� #nl_dvector_permute ����������� �������� ������� \f$y\f$
      �, ��� ����� ����� ������� ������� \f$x=P^{T}y\f$.

  ����� ����������, ���������, ��������, ������� #mda_create, #nl_xvector_create, #nl_dvector_create,
  ������� �������� ����������� ������ � ����� ��������� ��������� �������� ���������� �� 
  � ������� ������� #mda_free, #nl_xvector_free, #nl_dvector_free..


  \todo
   (����� ��� ����� � ������� �������� ��������������)
   - ��������� ������ � #sp_order. ��. xmda
   - ���������� #mda_order, ����� ��� �������� � �������� RR(U)O, � �� RR(C)O
   - ����� ����������� ������ ������ #mda_init � #mda_convert. ��� ���������������
     ���� (�� ������ ����� �� ��� �����) ������ �������� � #mda_order
   - ��������� ������ � #mda_chol_symb. ��. xmda__ ��� ������� � ������� #mda_chol_num
     ������ ������ ������ sparse
   - ��������� ������ � #mda_chol_solve. ��� ������� ������ ����� � ������ sparse
   ����� ����������� ���� ������ �-�� ������� ������� ����������:
    - #mda_order
    - #mda_chol_symb - ������� ����� ���������� ��-�������
    - #mda_chol_num - ������� ����� ���������� ��-�������
    - #mda_chol_solve - ������� ����� ���������� ��-�������
*/

/**
  \example xmda.c
*/


/**
  �������������� ������ ��� ��������, �������������� � MDA.

    - ����:
      - \f$n\f$ - ������� ������������ �������
      - \f$nz\f$ - ����� ��������� ��������� ��� ����������
    - �����: 
      - \f$IA\f$ - ������ ����� \f$n+1\f$
      - \f$JA\f$ - ������ ����� \f$2nz\f$
      - \f$D\f$ - ������ ����� \f$n\f$ 
      - \f$P\f$ - ������ ����� \f$n\f$ 
      - \f$IP\f$ - ������ ����� \f$n\f$
      - \f$M\f$ - ������ ����� \f$n\f$ 
      - \f$L\f$ - ������ ����� \f$n\f$
*/ 
extern void mda_create(
	size_t n,
	size_t nz,
	size_t **IA,
	size_t **JA,
	size_t **D,
	size_t **P,
	size_t **IP,
	size_t **M,
	size_t **L);

/**
  ������������ ������, ���������� �������� #mda_create.
*/ 
extern void mda_free(
	size_t *IA,
	size_t *JA,
	size_t *D,
	size_t *P,
	size_t *IP,
	size_t *M,
	size_t *L);



/**	
  ��������������� ������������� ��� ������ MDA.
  
    - ����:
      - \f$n\f$ - ������� ������� \f$A\f$
      - \f$IA\f$ - ������ �� RR(C)O-������������� �������� ������� \f$A\f$
    - �����:
      - \f$D\f$ - ������ �������� �����
      - \f$P\f$ - ������������������ (�������������) ������������ �����
      -	\f$IP\f$ - �������� ������������ (\f$IP[P[i]] = i\f$ � �������� \f$P[IP[i]] = i\f$)
*/
extern void mda_init(
	size_t n,
	size_t *IA,
	size_t *D,
	size_t *P,
	size_t *IP);

/**	
  �������������� ������������� ������������� � ������ � ������������� ��� ������ MDA.
  
  �������������� ������������� RR(U)U ������������� ������� \f$S\f$
  � ������ ������������� RR(C)O, � ����� ��������������� ������������� ���
  ������ MDA: ���������� �������� ����� \f$D\f$ � ������������� ��������� �������� 
  \f$P\f$ � \f$IP\f$. ��� ������������� ���� ������� #mda_init �������� ��� �� �������.
    - ����:
      - \f$n\f$ - ������� ������� \f$A\f$
      - \f$IS\f$, \f$JS\f$ - RR(U)U-������������� �������� ������������ ������� \f$A\f$
    - �����:
      - \f$IA\f$, \f$JA\f$ - RR(C)O-������������� �������� ������������ ������� \f$A\f$
      - \f$D\f$ - ������ �������� �����
      - \f$P\f$ - ������������������ (�������������) ������������ �����
      -	\f$IP\f$ - �������� ������������ (\f$IP[P[i]] = i\f$ � �������� \f$P[IP[i]] = i\f$)
*/
extern void mda_convert(
	size_t n,
	size_t *IS,
	size_t *JS,
	size_t *IA,
	size_t *JA,
	size_t *D,
	size_t *P,
	size_t *IP);

/*	
  ���� �������� MDA.
size_t mda_iterate(
	size_t h,
	size_t m,
	size_t n,
	size_t *IA,
	size_t *JA,
	size_t *M,
	size_t *L,
	size_t *D,
	size_t *P,
	size_t *IP);
*/

/**	
  ������������ ������� MDA.

    - ����:
      - \f$n\f$ - ������� ������� \f$A\f$
      - \f$IA\f$, \f$JA\f$ - RR(C)O-������������� �������� ������������ ������� \f$A\f$
      - \f$M\f$ - ��������������� ������ ����� \f$n\f$
      - \f$L\f$ - ��������������� ������ ����� \f$n\f$
      - \f$D\f$ - ������ �������� �����, ������������ �������� #mda_init ��� #mda_convert
      - \f$P\f$ - ������������������ ������������ �����, ������������ �������� #mda_init ��� #mda_convert
      -	\f$IP\f$ - �������� ������������, ������������ �������� #mda_init ��� #mda_convert
    - �����:
      - ������� ���������� ����� ��������� ��������� � ������� ����������� ����� ����������
        ���������
      - \f$D\f$ - ������ �������� �����
      - \f$P\f$ - ��������� ������������
      -	\f$IP\f$ - �������� ������������
*/
extern size_t mda_order(
	size_t n,
	size_t *IA,
	size_t *JA,
	size_t *M,
	size_t *L,
	size_t *D,
	size_t *P,
	size_t *IP);

/**
  ������������� ���������� ���������.
  ����� ����� ������: �� xspmda__.	
  ������� ������ ������������ � sparse
    - ����:
      - \f$n\f$ - ������� ������� \f$A\f$
      - \f$unz\f$ - 
      - \f$IS\f$, \f$JS\f$ - RR(U)O-������������� �������� ������������ ������� \f$A\f$
      - \f$M\f$ - ��������������� ������ ����� \f$n\f$
      - \f$L\f$ - ��������������� ������ ����� \f$n\f$
      - \f$D\f$ - ������ �������� �����, ������������ �������� #mda_order
      -	\f$IP\f$ - �������� ������������, ������������ �������� #mda_order
    - �����:
      - \f$IU\f$, \f$JU\f$ - RR(U)O-������������� �������� ���������� ��������� ������� \f$A\f$
*/
extern void mda_chol_symb(
	size_t n,
	size_t unz,
	size_t *IS,
	size_t *JS,
	size_t *IU,
	size_t *JU,
	size_t *M,
	size_t *L,
	size_t *D,
	size_t *IP);

/**	
  ��������� ���������� ���������.
  ������� ������ ������������ � sparse
    - ����:
      - \f$n\f$ - ������� ������� \f$A\f$
      - \f$IS\f$, \f$JS\f$, \f$SN\f$, \f$SD\f$ - RR(U)O-������������� ������������ 
          ������� \f$A\f$
      - \f$IU\f$, \f$JU\f$ - RR(U)O-������������� �������� ���������� ��������� ������� \f$A\f$,
          ������������ �������� #mda_chol_symb
      - \f$M\f$ - ��������������� ������ ����� \f$n\f$
      - \f$L\f$ - ��������������� ������ ����� \f$n\f$
      -	\f$P\f$ - ������������, ������������ �������� #mda_order
    - �����:
      - \f$UN\f$, \f$UD\f$ - �������� ��������� ���������� ��������� ������� \f$A\f$ 
          ��� ���������� � �� ��������� ��������������
*/
extern void mda_chol_num(
	size_t n,
	size_t *IS,
	size_t *JS,
	size_t *IU,
	size_t *JU,
	size_t *M,
	size_t *L,
	size_t *P,
	double *SN,
	double *SD,
	double *UN,
	double *UD);

/*	
  ������� ������� �� ������ ����������� ���������� ���������.

  �������� ������� \f$Ax = b\f$ �� ������ ����������� ���������� ���������.
  ���������� �� #sp_chol_solve �������������� ������������� ����������
  ������� �������� ������������ \f$P\f$
  
    - ����:
      - \f$n\f$ - ������� ������� \f$A\f$
      - \f$IU\f$, \f$JU\f$, \f$UN\f$, \f$UD\f$ - RR(U)O-������������� ���������� ��������� 
                                                 ������� \f$A\f$, ������������ ���������
                                                 #mda_chol_symb � #mda_chol_num
      -	\f$P\f$ - ������������, ������������ �������� #mda_order
      - \f$b\f$ - ������ ����� �������
    - �����:
      - \f$x\f$ - ��������� �������
*/
extern void mda_chol_solve(
	size_t n,
	size_t *IU,
	size_t *JU,
	size_t *P,
	double *UN,
	double *UD,
	double *b,
	double *x);

#endif
