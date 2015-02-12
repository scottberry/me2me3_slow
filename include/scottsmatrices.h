#ifndef SCOTTSMATRICES_H
#define SCOTTSMATRICES_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#define NEW(type) ((type *)calloc((size_t)1,sizeof(type)))
#define NEW_A(num,type) ((type *)calloc((size_t)(num),sizeof(type)))

#define TYPE_C_VEC 11
#define TYPE_C_MAT 12
#define TYPE_I_VEC 21
#define TYPE_I_MAT 22
#define TYPE_D_VEC 31
#define TYPE_D_MAT 32

typedef struct {
  char id;
  unsigned long len;
  signed char *el;
} C_VEC;

typedef struct {
  char id;
  unsigned long rows;
  unsigned long  cols;
  signed char **el;
} C_MAT;

typedef struct {
  char id;
  unsigned long len;
  signed long *el;
} I_VEC;

typedef struct {
  char id;
  unsigned long rows;
  unsigned long  cols;
  signed long **el;
} I_MAT;

typedef struct {
  char id;
  unsigned long len;
  double *el;
} D_VEC;

typedef struct {
  char id;
  unsigned long rows;
  unsigned long  cols;
  double **el;
} D_MAT;

C_VEC *c_vec_get(unsigned long);
int c_vec_free(C_VEC *);
int c_vec_print(FILE *,C_VEC *);

I_VEC *i_vec_get(unsigned long);
int i_vec_free(I_VEC *);
int i_vec_print(FILE *,I_VEC *);

D_VEC *d_vec_get(unsigned long);
int d_vec_free(D_VEC *);
int d_vec_print(FILE *,D_VEC *);

C_MAT *c_mat_get(unsigned long, unsigned long);
int c_mat_free(C_MAT *);
int c_mat_print(FILE *,C_MAT *);

I_MAT *i_mat_get(unsigned long, unsigned long);
int i_mat_free(I_MAT *);
int i_mat_print(FILE *,I_MAT *);

D_MAT *d_mat_get(unsigned long, unsigned long);
int d_mat_free(D_MAT *);
int d_mat_print(FILE *,D_MAT *);

#endif
