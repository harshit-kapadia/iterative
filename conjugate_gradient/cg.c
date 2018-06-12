#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include <math.h>
#include <time.h>

void mat_vect_CSR(int *IA, int *JA, double *val, double *x, double *y, int M) ;
void mat_vect_CSR_2(int *IA, int *JA, double *val, double **x, double *y, int M, int colm) ;
void mat_vect_CSC_modified(int *IA, int *JA, double *val, double *x, double *y, int M) ;
void mat_vect(double **A, double *x, double *y/*output*/, int m/*row*/, int n) ;

void mat_vect_CSR_CSC(int *IA, int *JA, double *val, double *x, double *y, int M) ;
double dot_product(double *x, double *y/*output*/, int n) ;
double copy_vect(double *x /* data to copy */, double *y/*output*/, int n) ;
void diff_vect(double *x /* data to copy */, double *y/*output*/, double *diff, int n) ;
void sum_vect(double *x /* data to copy */, double *y/*output*/, double *sum, int n) ;
void scalar_vector(double a, double *x, double *a_x, int n) ;

int main(int argc, char *argv[])
{
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, nz;
  int i, *I, *J;
  double *val;

  if (argc < 2)
  {
      fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
      exit(1);
  }
  else
  {
      if ((f = fopen(argv[1], "r")) == NULL)
          exit(1);
  }

  if (mm_read_banner(f, &matcode) != 0)
  {
      printf("Could not process Matrix Market banner.\n");
      exit(1);
  }

  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
      mm_is_sparse(matcode))
  {
      printf("Sorry, this application does not support ");
      printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
      exit(1);
  }

  /* find out size of sparse matrix .... */

  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
      exit(1);

  /* reseve memory for matrices */

  I = (int *)malloc(nz * sizeof(int));
  J = (int *)malloc(nz * sizeof(int));
  val = (double *)malloc(nz * sizeof(double));

  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  for (i = 0; i < nz; i++)
  {
      fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
      I[i]--; /* adjust from 1-based to 0-based */
      J[i]--;
  }

  if (f != stdin)
      fclose(f);

  printf("Matrix read successfully. Starting conversion to CSR format...");

  
  int j, d1, d2 ;
  double d3 ;
  /*Sorting in required COO format*/
  for(i=0; i<nz; ++i){
      for(j=(i+1); j<nz; ++j){
      if(I[i]>I[j]){
          d1 = I[i] ;        d2 = J[i] ;        d3 = val[i] ;
          I[i] = I[j] ;        J[i] = J[j] ;        val[i] = val[j] ;
          I[j] = d1 ;        J[j] = d2 ;        val[j] = d3 ;
      }
      else if(I[i]==I[j]){
          if(J[i]>J[j]){
          d1 = I[i] ;        d2 = J[i] ;        d3 = val[i] ;
          I[i] = I[j] ;        J[i] = J[j] ;        val[i] = val[j] ;
          I[j] = d1 ;        J[j] = d2 ;        val[j] = d3 ;
          }
      }
      }
  }

  int *IA, *JA ;
  IA = (int *) malloc((M+1) * sizeof(int));
  JA = (int *) malloc(nz * sizeof(int));

  IA[0] = 0 ;   JA[0] = J[0] ;
  j = 1 ;
  for(i=1; i<nz; ++i){
      if(I[i]!=I[i-1]){
      IA[j] = i ;
      j = j + 1 ;
      }
      JA[i] = J[i] ;
  }
  if(j!=M){
      printf("\n Error while converting to CSR format! \n") ;
      return 0 ;
  }
  else
      IA[M] = nz ;

  printf("\n Conversion successful! \n") ;

  // FILE *fp_test ;
  // fp_test = fopen("matrix-csr.txt", "w") ;
  // for(i=0; i<M; ++i){
  //   for(j=IA[i]; j<IA[i+1]; ++j)
  //     fprintf(fp_test, "%d  %d  %20.19g\n", i+1, JA[j]+1, val[j]);
  //     // fprintf(fp_test, "%d  \n", i+1 >= JA[j]+1); // checking whether row index >= column, giving lower triangular mat.
  // }
  // fclose(fp_test) ;
  
  /* Conversion to CSR format complete */


  double *x0, *b, *x_star /* , *xm */ ;
  double norm_initial_residual ;

  x0 = (double*)malloc((M)*sizeof(double)) ;
  b = (double*)malloc((M)*sizeof(double)) ;
  // xm = (double*)malloc((M)*sizeof(double)) ;
  x_star = (double*)malloc((M)*sizeof(double)) ;

  for(i=0; i<M; ++i)
    x_star[i] = 1 ;    /* Exact Solution Vector */

  mat_vect_CSR_CSC(IA, JA, val, x_star, b, M) ;   /* Getting RHS, b = A*x_star */
  // for(i=0; i<M; ++i)
  //   printf(" %lf\n", b[i]) ;

  for(i=0; i<M; ++i)
    x0[i] = 0.0 ;

  norm_initial_residual = dot_product(b, b, M) ;  /* r0 = b - A*x0 = b - A*0 = b */
  norm_initial_residual = sqrt(norm_initial_residual) ;
  // printf(" %lf \n", norm_initial_residual) ;

  double *r, *p, *A_e, *temp, *x, *e, *A_p, *alpha_p, *alpha_A_p, *beta_p ;
  double rho, beta, alpha ;
  double e_normA, r_norm2 ;
  // int m = 100000 ; /* arbitrary maximum iteration value, m>=1 */
  r = (double*)malloc((M)*sizeof(double)) ;
  p = (double*)malloc((M)*sizeof(double)) ;
  A_e = (double*)malloc((M)*sizeof(double)) ;
  temp = (double*)malloc((M)*sizeof(double)) ;
  x = (double*)malloc((M)*sizeof(double)) ;
  e = (double*)malloc((M)*sizeof(double)) ;
  A_p = (double*)malloc((M)*sizeof(double)) ;
  alpha_p = (double*)malloc((M)*sizeof(double)) ;
  alpha_A_p = (double*)malloc((M)*sizeof(double)) ;
  beta_p = (double*)malloc((M)*sizeof(double)) ;
  // e_normA = (double*)malloc((m)*sizeof(double)) ;
  // r_norm2 = (double*)malloc((m)*sizeof(double)) ;

  diff_vect(x,x_star,e,M) ; // e = x-x_star
  mat_vect_CSR_CSC(IA, JA, val, e, A_e, M) ;   /* Getting A_e = A*e */
  e_normA = dot_product(A_e,e,M) ; /* returns (A_e,e) */
  e_normA = sqrt(e_normA) ;

  r_norm2 = norm_initial_residual ; // r_norm2 = (r,r)

  copy_vect(x0,x,M) ; // copy x0 in x
  copy_vect(b,r,M) ; /* r0 = b - A*x0 = b - A*0 = b */
  copy_vect(r,p,M) ; /* p = r = r0 */
  rho = norm_initial_residual * norm_initial_residual ;

  // beta = 10.0 ; 
  int iteration = 0 ;
  double tol = 1e-8 , check = 1.0 ; /* arbitrary value > 1e-08 for 1st iteration of the loop */

  FILE *fp_cg ;
  fp_cg = fopen("cg.txt", "w") ;
  fprintf(fp_cg, "Iteration-index\t\te_normA\t\tr_norm2\n") ;
  fprintf(fp_cg, "%d\t\t%.10lf\t\t%.10lf\n", iteration, e_normA, r_norm2) ;
  fclose(fp_cg) ;

  while(check > tol){
    iteration += 1 ;
    printf("\n Starting Iteration %d \n", iteration) ;

    mat_vect_CSR_CSC(IA, JA, val, p, A_p, M) ;
    alpha = dot_product(A_p, p, M) ;
    alpha = rho / alpha ;
    
    scalar_vector(alpha, p, alpha_p, M) ; // alpha*p
    sum_vect(x, alpha_p, x, M) ; // x = x + alpha_p

    diff_vect(x, x_star, e, M) ; // e = x-x_star
    mat_vect_CSR_CSC(IA, JA, val, e, A_e, M) ;   /* Getting A_e = A*e */
    e_normA = dot_product(A_e, e, M) ; /* returns (A_e,e) */
    e_normA = sqrt(e_normA) ;

    scalar_vector(alpha, A_p, alpha_A_p, M) ; // computes alpha*A*p
    diff_vect(r, alpha_A_p, r, M) ;

    beta = rho ;
    rho = dot_product(r, r, M) ;
    beta = rho / beta ;

    r_norm2 = sqrt(rho) ;

    scalar_vector(beta, p, beta_p, M) ; // computes beta*p
    sum_vect(r, beta_p, p, M) ;

    fp_cg = fopen("cg.txt", "a") ;
    fprintf(fp_cg, "%d\t\t%.10lf\t\t%.10lf\n", iteration, e_normA, r_norm2) ;
    fclose(fp_cg) ;

    check = r_norm2 / norm_initial_residual ;
    printf("check : %.10lf", check) ;
  }

  printf("Final iteration number is %d.", iteration) ;

  free(r);
  free(p);
  free(A_e);
  free(temp);
  free(x);
  free(e);
  free(A_p);
  free(alpha_p);
  free(alpha_A_p);
  free(beta_p);
  free(x0);
  free(b);
  free(x_star);
  free(IA);
  free(JA);
  free(I);
  free(J);
  free(val);

  return 0;
}


void mat_vect_CSR(int *IA, int *JA, double *val, double *x, double *y, int M){
  int i, j, i1, i2 ;
  for(i=0; i<M; ++i){
    y[i] = 0 ;
    i1 = IA[i] ;
    i2 = IA[i+1] - 1 ;
    for(j=i1; j<=i2; ++j)
      y[i] = y[i] + ( val[j] * x[JA[j]] ) ;
  }
}

void mat_vect_CSR_2(int *IA, int *JA, double *val, double **x, double *y, int M, int colm){
  int i, j, i1, i2 ;
  for(i=0; i<M; ++i){
    y[i] = 0 ;
    i1 = IA[i] ;
    i2 = IA[i+1] - 1 ;
    for(j=i1; j<=i2; ++j)
      y[i] = y[i] + ( val[j] * x[JA[j]][colm] ) ;
  }
}

void mat_vect_CSC_modified(int *IA, int *JA, double *val, double *x, double *y, int M){
  int i, j, i1, i2 ;
  for(i=0; i<M; ++i)
    y[i] = 0 ;
  for(i=0; i<M; ++i){
    i1 = IA[i] ;
    i2 = IA[i+1] - 1 ;
    for(j=i1; j<i2; ++j)    /* Ignoring the diagonal --> not considering j = i2 */
      y[JA[j]] = y[JA[j]] + ( val[j] * x[i] ) ;
  }
}

void mat_vect_CSR_CSC(int *IA, int *JA, double *val, double *x, double *y, int M){
  double *y1, *y2 ;
  y1 = (double*)malloc((M)*sizeof(double)) ;
  y2 = (double*)malloc((M)*sizeof(double)) ;
  
  mat_vect_CSR(IA, JA, val, x, y1, M) ; // compute A*x (lower triangular A)
  mat_vect_CSC_modified(IA, JA, val, x, y2, M) ; // compute A^(transpose)*x upon omitting diagonal
  sum_vect(y1,y2,y,M) ;
}

void mat_vect(double **A, double *x, double *y/*output*/, int m/*row*/, int n){
  int i, j ;
  for(i=0; i<m; ++i){
    y[i] = 0 ;
    for(j=0; j<n; ++j)
      y[i] = y[i] + ( A[i][j] * x[j] ) ;
  }
}

double dot_product(double *x, double *y/*output*/, int n){
  double ret = 0.0;
  int i;
  for(i=0; i<n; ++i)
    ret = ret + ( x[i] * y[i] ) ;

  return ret;
}

/* copy from x to y */
double copy_vect(double *x /* data to copy */, double *y/*output*/, int n){
  int i;
  for(i=0; i<n; ++i)
    y[i] = x[i] ;
}

/* diff = x - y */
void diff_vect(double *x /* data to copy */, double *y/*output*/, double *diff, int n){
  int i;
  for(i=0; i<n; ++i)
    diff[i] = x[i] - y[i] ;
}

/* sum = x + y */
void sum_vect(double *x /* data to copy */, double *y/*output*/, double *sum, int n){
  int i;
  for(i=0; i<n; ++i)
    sum[i] = x[i] + y[i] ;
}

void scalar_vector(double a, double *x, double *a_x, int n){
  int i;
  for(i=0; i<n; ++i)
    a_x[i] = ( a * x[i] ) ;
}