#include "scottsmatrices.h"

C_VEC *c_vec_get(unsigned long len) {
  C_VEC *vector;

  vector = NEW(C_VEC);
  vector->id = TYPE_C_VEC;
  vector->len = len;

  if ( (vector->el = (signed char *)calloc(len,sizeof(char))) == 
      (signed char *)NULL ) {
    fprintf(stderr,"c_vec_get failed\n");
    exit(-1);
  }
  return(vector);
}

int c_vec_free(C_VEC *vector) {

  if ( vector==(C_VEC *)NULL || (signed char)(vector->len) < 0 ) 
    return(-1);

  if ( vector->el != (signed char *)NULL )
    free((signed char *)(vector->el));

  free((signed char *)vector);
  return(0);
}

int c_vec_print(FILE *fp, C_VEC *vector) {
  unsigned long i;

  if ( vector==(C_VEC *)NULL ) {
    fprintf(stderr,"i_vec_print: vector is NULL");
    return(-1);
  }
  for (i=0;i<vector->len;i++) {
    fprintf(fp,"%d",vector->el[i]);
    if (i!=(vector->len-1)) fprintf(fp,"\t");
  }
  fprintf(fp,"\n");
  return(1);
}

I_VEC *i_vec_get(unsigned long len) {
  I_VEC *vector;

  vector = NEW(I_VEC);
  vector->id = TYPE_I_VEC;
  vector->len = len;

  if ( (vector->el = (long *)calloc(len,sizeof(long))) == 
      (long *)NULL ) {
    fprintf(stderr,"i_vec_get failed\n");
    exit(-1);
  }
  return(vector);
}

int i_vec_free(I_VEC *vector) {

  if ( vector==(I_VEC *)NULL || (long)(vector->len) < 0 ) 
    return(-1);

  if ( vector->el != (long *)NULL )
    free((signed char *)(vector->el));

  free((char *)vector);
  return(0);
}

int i_vec_print(FILE *fp, I_VEC *vector) {
  unsigned long i;

  if ( vector==(I_VEC *)NULL ) {
    fprintf(stderr,"i_vec_print: vector is NULL");
    return(-1);
  }
  for (i=0;i<vector->len;i++) {
    fprintf(fp,"%ld ",vector->el[i]);
    if (i!=(vector->len-1)) fprintf(fp,"\t");
  }
  fprintf(fp,"\n");
  return(1);
}

D_VEC *d_vec_get(unsigned long len) {
  D_VEC *vector;

  vector = NEW(D_VEC);
  vector->id = TYPE_D_VEC;
  vector->len = len;

  if ( (vector->el = (double *)calloc(len,sizeof(double))) == 
      (double *)NULL ) {
    fprintf(stderr,"d_vec_get failed\n");
    exit(-1);
  }
  return(vector);
}

int d_vec_free(D_VEC *vector) {

  if ( vector==(D_VEC *)NULL || (double)(vector->len) < 0 ) 
    return(-1);

  if ( vector->el != (double *)NULL )
    free((char *)(vector->el));

  free((char *)vector);
  return(0);
}

int d_vec_print(FILE *fp, D_VEC *vector) {
  unsigned long i;

  if ( vector==(D_VEC *)NULL ) {
    fprintf(stderr,"d_vec_print: vector is NULL");
    return(-1);
  }
  for (i=0;i<vector->len;i++) {
    fprintf(fp,"%f ",vector->el[i]);
    if (i!=(vector->len-1)) fprintf(fp,"\t");
  }
  fprintf(fp,"\n");
  return(1);
}

C_MAT *c_mat_get(unsigned long m, unsigned long n) {
  C_MAT *matrix;
  unsigned long i;

  matrix = NEW(C_MAT);
  matrix->id = TYPE_C_MAT;
  matrix->rows = m;
  matrix->cols = n;

  if ( (matrix->el = (signed char **)calloc(m,sizeof(char *))) == 
      (signed char **)NULL ) {
    fprintf(stderr,"c_mat_get failed\n");
    exit(-1);
  }

  for (i=0;i<m;i++) {
    if( (matrix->el[i]=NEW_A(n,signed char)) == (signed char *)NULL ) {
      fprintf(stderr,"c_mat_get failed\n");
      exit(-1);
    }
  }
  return(matrix);
}

int c_mat_free(C_MAT *matrix) {
   unsigned long i;

  if ( matrix==(C_MAT *)NULL || (signed char)(matrix->cols) < 0 || 
       (signed char)(matrix->rows) < 0 ) 
    return(-1);

  for (i=0;i<matrix->rows;i++) 
    free((char *)(matrix->el[i]));

  if (matrix->el != (signed char **)NULL )
    free((char *)(matrix->el));

  free((char *)matrix);
  return(0);
}

int c_mat_print(FILE *fp, C_MAT *matrix) {
  unsigned long i,j;

  if ( matrix==(C_MAT *)NULL ) {
    fprintf(stderr,"c_mat_print: matrix is NULL");
    return(-1);
  }
  for (i=0;i<matrix->rows;i++) {
    for (j=0;j<matrix->cols;j++) {
      fprintf(fp,"%d",matrix->el[i][j]);
      if (j!=(matrix->cols-1)) fprintf(fp,"\t");
    }
    fprintf(fp,"\n");
  }
  return(1);
}

I_MAT *i_mat_get(unsigned long m, unsigned long n) {
  I_MAT *matrix;
  unsigned long i;

  matrix = NEW(I_MAT);
  matrix->id = TYPE_I_MAT;
  matrix->rows = m;
  matrix->cols = n;

  if ( (matrix->el = (long **)calloc(m,sizeof(long *))) == 
      (long **)NULL ) {
    fprintf(stderr,"i_mat_get failed\n");
    exit(-1);
  }

  for (i=0;i<m;i++) {
    if( (matrix->el[i]=NEW_A(n,long)) == (long *)NULL ) {
      fprintf(stderr,"i_mat_get failed\n");
      exit(-1);
    }
  }
  return(matrix);
}

int i_mat_free(I_MAT *matrix) {
  long unsigned i;

  if ( matrix==(I_MAT *)NULL || (long)(matrix->cols) < 0 || 
       (long)(matrix->rows) < 0 ) 
    return(-1);

  for (i=0;i<matrix->rows;i++) 
    free((char *)(matrix->el[i]));

  if (matrix->el != (long **)NULL )
    free((char *)(matrix->el));

  free((char *)matrix);
  return(0);
}

int i_mat_print(FILE *fp, I_MAT *matrix) {
  unsigned long i,j;

  if ( matrix==(I_MAT *)NULL ) {
    fprintf(stderr,"i_mat_print: matrix is NULL");
    return(-1);
  }
  for (i=0;i<matrix->rows;i++) {
    for (j=0;j<matrix->cols;j++) {
      fprintf(fp,"%ld",matrix->el[i][j]);
      if (j!=(matrix->cols-1)) fprintf(fp,"\t");
    }
    fprintf(fp,"\n");
  }
  return(1);
}

D_MAT *d_mat_get(unsigned long m, unsigned long n) {
  D_MAT *matrix;
  unsigned long i;

  matrix = NEW(D_MAT);
  matrix->id = TYPE_D_MAT;
  matrix->rows = m;
  matrix->cols = n;

  if ( (matrix->el = (double **)calloc(m,sizeof(double *))) == 
      (double **)NULL ) {
    fprintf(stderr,"d_mat_get failed\n");
    exit(-1);
  }

  for (i=0;i<m;i++) {
    if( (matrix->el[i]=NEW_A(n,double)) == (double *)NULL ) {
      fprintf(stderr,"d_mat_get failed\n");
      exit(-1);
    }
  }
  return(matrix);
}

int d_mat_free(D_MAT *matrix) {
  long unsigned i;

  if ( matrix==(D_MAT *)NULL || (double)(matrix->cols) < 0 || 
       (double)(matrix->rows) < 0 ) 
    return(-1);

  for (i=0;i<matrix->rows;i++) 
    free((char *)(matrix->el[i]));

  if (matrix->el != (double **)NULL )
    free((char *)(matrix->el));

  free((char *)matrix);
  return(0);
}

int d_mat_print(FILE *fp, D_MAT *matrix) {
  unsigned long i,j;

  if ( matrix==(D_MAT *)NULL ) {
    fprintf(stderr,"d_mat_print: matrix is NULL");
    return(-1);
  }
  for (i=0;i<matrix->rows;i++) {
    for (j=0;j<matrix->cols;j++) {
      fprintf(fp,"%f",matrix->el[i][j]);
      if (j!=(matrix->cols-1)) fprintf(fp,"\t");
    }
    fprintf(fp,"\n");
  }
  return(1);
}

  
