
int read_sparse_to_CSR(int argc, char *argv[], int *IA, int *JA, double *val);

void mat_vect_CSR(int *IA, int *JA, double *val, double *x, double *y, int M);
void mat_vect_CSR_2(int *IA, int *JA, double *val, double **x, double *y, int M, int colm);
void mat_vect_CSC_modified(int *IA, int *JA, double *val, double *x, double *y, int M);
void mat_vect(double **A, double *x, double *y /*output*/, int m /*row*/, int n);

void mat_vect_CSR_CSC(int *IA, int *JA, double *val, double *x, double *y, int M);
double dot_product(double *x, double *y /*output*/, int n);
double copy_vect(double *x /* data to copy */, double *y /*output*/, int n);
void diff_vect(double *x /* data to copy */, double *y /*output*/, double *diff, int n);
void sum_vect(double *x /* data to copy */, double *y /*output*/, double *sum, int n);
void scalar_vector(double a, double *x, double *a_x, int n);
double norm(double *vector, int n);
void print_vector_double(double *vector, int n);
void print_vector_int(int *vector, int n);

double power_iteration(int *IA, int *JA, double *valA, double *q, int M);
double lanczos(int *IA, int *JA, double *valA, double *q, int M, int dim);



//****************************************************************************************//
double power_iteration(int *IA, int *JA, double *valA, double *q, int M)
{
	int i;
	double tol = 1e-10; /* change tolerance for Lanczos method, as needed (else 1e-8 for STEP-1) */
	double z_norm, eig = 0.0, eig_old = 0.0, res = 1.0;

	double *z = malloc(M * sizeof(double));

	double eig_true = 9598.6080894852857; /* max. eigenvalue of s3rmt3m3.mtx */

	/* Comment file creation when timing */
	// FILE *f1;
	// f1 = fopen("power-iteration.txt", "w");
	// fprintf(f1, "Iteration\t\tResidual\n");

	int iteration = 1;

	// clock_t start = clock(); /* Comment when timing for Lanczos method */

	// needed only for iteration=1; so writing outside while to reduce one if() statement
	for (i = 0; i < M; ++i)
		z[i] = 0.0;
	mat_vect_CSR_CSC(IA, JA, valA, q, z, M); /* Getting z = A*q */

	while (res > tol)	{
	// while (iteration <= 5000) {   /* useful for STEP-2 */

		z_norm = norm(z, M);

		for (i = 0; i < M; ++i)
		{
			q[i] = z[i] / z_norm;
			z[i] = 0.0; /* making it zero as we will be computing z for next iteration via mat. vect. product */
		}

		mat_vect_CSR_CSC(IA, JA, valA, q, z, M); /* Getting z = A*q */

		eig_old = eig; /* only comment when running for some iterations and want to time */
		eig = dot_product(q, z, M); /* computing ( q^T * (A*q) ) */

		res = fabs(eig - eig_old); /* only comment when running for some iterations and want to time */

		// fprintf(f1, "%d\t\t%e\n", iteration, res); /* Comment file print when timing */
		iteration++;
	}

	// clock_t end = clock(); /* Comment when timing for Lanczos method */
	// printf("Power Iteration time : %e\t error: %e\n", (end - start)/(double)CLOCKS_PER_SEC, eig_true - eig);

	free(z);

	return eig;
}

//****************************************************************************************//
double lanczos(int *IA, int *JA, double *valA, double *q, int M, int dim)
{
	int i, j, index;

	double *w = malloc(M * sizeof(double));
	double *v_old = malloc(M * sizeof(double));
	double *v_new = malloc(M * sizeof(double));
	double *product1 = malloc(M * sizeof(double));
	double *product2 = malloc(M * sizeof(double));
	double eig = 0.0;
	double beta_old = 0.0, beta_new = 0.0, alpha = 0.0;

	double *val = malloc((2 * dim - 1) * sizeof(double));
	int *IT = malloc((dim + 1) * sizeof(int));
	int *JT = malloc((2 * dim - 1) * sizeof(int));

	double eig_true = 9598.6080894852857; /* max. eigenvalue of s3rmt3m3.mtx */


	for (i = 0; i < M; ++i)
	{
		v_old[i] = 0.0;
		v_new[i] = q[i];
	}

	clock_t start = clock();

	for (i = 0; i < dim; ++i)
	{
		beta_old = beta_new;

		for (j = 0; j < M; ++j)
		{
			product1[j] = 0.0;
			product2[j] = 0.0;
			w[j] = 0.0;
		}

		mat_vect_CSR_CSC(IA, JA, valA, v_new, product1, M); /* product1 = A*v_new */
		scalar_vector(beta_new, v_old, product2, M); // beta_new * v_old
		diff_vect(product1, product2, w, M);	// w = product1 - product2

		alpha = dot_product(w, v_new, M);	/* returns (w,v_new) */

		for (j = 0; j < M; ++j)
		{
			product1[j] = 0.0;
		}

		scalar_vector(alpha, v_new, product1, M); // alpha * v_new
		diff_vect(w, product1, w, M);	// w = w - product1

		beta_new = norm(w, M);

		for (j = 0; j < M; ++j)
		{
			v_old[j] = v_new[j];
			v_new[j] = w[j] / beta_new;
		}

		if (i == 0)
		{
			val[i] = alpha;
			JT[i] = 0;
			IT[i] = 0;
			index = 1;
		}
		else
		{
			IT[i] = index;
			// IT[i + 1] = 3 * dim - 2;
			JT[index] = i - 1;
			JT[index + 1] = i;
			val[index] = beta_old;
			val[index + 1] = alpha;
			index = index + 2;
		}
	}
	IT[dim] = 2*dim-1 ; /* last entry in IT */

	/* For testing */
	// printf("\n IT : \n");
	// print_vector_int(IT, dim+1);

	// printf("\n JT : \n");
	// print_vector_int(JT, (2 * dim - 1));

	// printf("\n val : \n");
	// print_vector_double(val, (2 * dim - 1));

	/* For testing */
	// FILE *fp_test ;
	// fp_test = fopen("T-matrix-csr.txt", "w") ;
	// for(i=0; i<dim; ++i){
    // 	for(j=IT[i]; j<IT[i+1]; ++j)
    // 		fprintf(fp_test, "%d  %d  %20.19g\n", i+1, JT[j]+1, val[j]);
    //    		// fprintf(fp_test, "%d  \n", i+1 >= JA[j]+1); // checking whether row index >= column, giving lower triangular mat.
	// }
	// fclose(fp_test) ;


	double *q_T = malloc(dim * sizeof(double));
	double sqrt_dim = sqrt(dim);

	for (i = 0; i < dim; i++)
	{
		q_T[i] = 1.0 / sqrt_dim;
	}

	eig = power_iteration(IT, JT, val, q_T, dim) ;

	clock_t end = clock();
	printf("Lanczos time : %e\t error: %e\n", (end - start) / (double)CLOCKS_PER_SEC, eig_true - eig);


	free(w);
	free(v_old);
	free(v_new);
	free(product1);
	free(product2);
	free(q_T);

	return eig;
}

//****************************************************************************************//
int read_sparse_to_CSR(int argc, char *argv[], int *IA, int *JA, double *val)
{
	int ret_code;
	MM_typecode matcode;
	FILE *f;
	int M, N, nz;
	int i, *I, *J;

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

	int j, d1, d2;
	double d3;
	/*Sorting in required COO format*/
	for (i = 0; i < nz; ++i)
	{
		for (j = (i + 1); j < nz; ++j)
		{
			if (I[i] > I[j])
			{
				d1 = I[i];
				d2 = J[i];
				d3 = val[i];
				I[i] = I[j];
				J[i] = J[j];
				val[i] = val[j];
				I[j] = d1;
				J[j] = d2;
				val[j] = d3;
			}
			else if (I[i] == I[j])
			{
				if (J[i] > J[j])
				{
					d1 = I[i];
					d2 = J[i];
					d3 = val[i];
					I[i] = I[j];
					J[i] = J[j];
					val[i] = val[j];
					I[j] = d1;
					J[j] = d2;
					val[j] = d3;
				}
			}
		}
	}

	IA[0] = 0;
	JA[0] = J[0];
	j = 1;
	for (i = 1; i < nz; ++i)
	{
		if (I[i] != I[i - 1])
		{
			IA[j] = i;
			j = j + 1;
		}
		JA[i] = J[i];
	}
	if (j != M)
	{
		printf("\n Error while converting to CSR format! \n");
		return 0;
	}
	else
		IA[M] = nz;

	printf("\n Conversion successful! \n\n");
}

//****************************************************************************************//
void mat_vect_CSR(int *IA, int *JA, double *val, double *x, double *y, int M)
{
	int i, j, i1, i2;
	for (i = 0; i < M; ++i)
	{
		y[i] = 0;
		i1 = IA[i];
		i2 = IA[i + 1] - 1;
		for (j = i1; j <= i2; ++j)
			y[i] += (val[j] * x[JA[j]]);
	}
}

//****************************************************************************************//
void mat_vect_CSR_2(int *IA, int *JA, double *val, double **x, double *y, int M, int colm)
{
	int i, j, i1, i2;
	for (i = 0; i < M; ++i)
	{
		y[i] = 0;
		i1 = IA[i];
		i2 = IA[i + 1] - 1;
		for (j = i1; j <= i2; ++j)
			y[i] += (val[j] * x[JA[j]][colm]);
	}
}

//****************************************************************************************//
void mat_vect_CSC_modified(int *IA, int *JA, double *val, double *x, double *y, int M)
{
	int i, j, i1, i2;
	for (i = 0; i < M; ++i)
		y[i] = 0;
	for (i = 0; i < M; ++i)
	{
		i1 = IA[i];
		i2 = IA[i + 1] - 1;
		for (j = i1; j < i2; ++j) /* Ignoring the diagonal --> not considering j = i2 */
			y[JA[j]] += (val[j] * x[i]);
	}
}

//****************************************************************************************//
void mat_vect_CSR_CSC(int *IA, int *JA, double *val, double *x, double *y, int M)
{
	double *y1, *y2;
	y1 = (double *)malloc((M) * sizeof(double));
	y2 = (double *)malloc((M) * sizeof(double));

	mat_vect_CSR(IA, JA, val, x, y1, M);					// compute A*x (lower triangular A)
	mat_vect_CSC_modified(IA, JA, val, x, y2, M); // compute A^(transpose)*x upon omitting diagonal
	sum_vect(y1, y2, y, M);
}

//****************************************************************************************//
void mat_vect(double **A, double *x, double *y /*output*/, int m /*row*/, int n)
{
	int i, j;
	for (i = 0; i < m; ++i)
	{
		y[i] = 0;
		for (j = 0; j < n; ++j)
			y[i] += (A[i][j] * x[j]);
	}
}

//****************************************************************************************//
double dot_product(double *x, double *y /*output*/, int n)
{
	double ret = 0.0;
	int i;
	for (i = 0; i < n; ++i)
		ret += (x[i] * y[i]);

	return ret;
}

//****************************************************************************************//
/* copy from x to y */
double copy_vect(double *x /* data to copy */, double *y /*output*/, int n)
{
	int i;
	for (i = 0; i < n; ++i)
		y[i] = x[i];
}

//****************************************************************************************//
/* diff = x - y */
void diff_vect(double *x /* data to copy */, double *y /*output*/, double *diff, int n)
{
	int i;
	for (i = 0; i < n; ++i)
		diff[i] = x[i] - y[i];
}

//****************************************************************************************//
/* sum = x + y */
void sum_vect(double *x /* data to copy */, double *y /*output*/, double *sum, int n)
{
	int i;
	for (i = 0; i < n; ++i)
		sum[i] = x[i] + y[i];
}

//****************************************************************************************//
void scalar_vector(double a, double *x, double *a_x, int n)
{
	int i;
	for (i = 0; i < n; ++i)
		a_x[i] = (a * x[i]);
}

//****************************************************************************************//
double norm(double *vector, int n)
{
	int i;
	double NORM = 0.0;

	for (i = 0; i < n; ++i)
	{
		NORM += vector[i] * vector[i];
	}
	return sqrt(NORM);
}

//****************************************************************************************//
void print_vector_double(double *vector, int n)
{
	for(int i=0; i<n; ++i)
		printf(" %lf\n", vector[i]) ;
}

//****************************************************************************************//
void print_vector_int(int *vector, int n)
{
	for(int i=0; i<n; ++i)
		printf(" %d\n", vector[i]) ;
}