#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include <math.h>
#include <time.h>
#include "functions.h"


int main(int argc, char *argv[])
{
	// int M = 675, N = 675, nz = 1965 ; // for NOS6 matrix
	int M = 5357, N = 5357, nz = 106526; // for S3RMT3M3 matrix

	int *IA, *JA;
	double *val;
	IA = (int *)malloc((M + 1) * sizeof(int));
	JA = (int *)malloc(nz * sizeof(int));
	val = (double *)malloc(nz * sizeof(double));

	read_sparse_to_CSR(argc, argv, IA, JA, val);


	int i ;
	double *q = malloc(M * sizeof(double));

	for (i = 0; i < M; i++)
	{
		q[i] = 1.0 / sqrt(M);
	}

	/* Uncomment Power Iteration OR Lanczos method line, as needed */
	// printf("Power Iteration: max. eigenvalue = %.10lf\n\n", power_iteration(IA, JA, val, q, M));
	printf("Lanczos: max. eigenvalue = %.10lf\n\n", lanczos(IA, JA, val, q, M, 100));

	free(IA);
	free(JA);
	free(val);

	return 0;
}