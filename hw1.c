/* Author: Seoyun Lee
*  File Name: hw1.c
*  Date: 2/13/14
*  
*  Description: This program includes the calling sequence
*  void inv_double_gs(double *a, int n, double *u, double *b)
*  that inverts a given matrix of size nxn, where a points to
*  an array of length n that contains the matrix in row order
*  (so first n elements are the first row of a, etc.)
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void inv_double_gs(double *a, int n, double *u, double *b);
void normalize(double *a, int col, int n);
void transpose(double *a, double *u, int n);
double msize(double *a, int col, int n);
double innerproduct(double *u, double *a, int colu, int cola, int n);
void projection(double *u, double *a, int colu, int cola, int n, double *proj);
void makeg(double *g, int col, int n, double sizea);
void transpose(double *u, double *newu, int n);
void multiply(double *g, double *u, double *b, int n);
void multiplyg(double *g, int colg, int n, double constant, double *projg);
double constant(double *u, double *a, int colu, int cola, int n);

int main(){
	return 0;
}

/*
	*a: points to an array on size n*n that contains the elements of the matrix to be inverted
	n: size of the matrix
	*u: points to a after it has been orthogonalized by column (memory allocated by user)
	*b: points to the inverse of the matrix a (memory allocated by user)
	given a pointer a to an n x n matrix it will give the matrix orthogonalized by column
	via gram schmidt in u and the inverse in b
*/
void inv_double_gs(double *a, int n, double *u, double *b) {

	int i;
	int j;
	for(i = 0; i < n * n; i++) {
		u[i] = a[i];
	}
	/* copy a into u */

	double *g;
	g = malloc(n * n * sizeof(*g));

	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			if(i == j)
				g[i*n + j] = 1;
			else
				g[i*n + j] = 0;
		}
	}
	/* g contains the identity matrix now */

	for(i = 0; i < n; i++) {
		makeg(g, i, n, msize(u, i, n));
		normalize(u, i, n); /* normalize column i */
		for(j = i + 1; j < n; j++) {
			double *proj;
			proj = malloc(n * sizeof(*proj));
			projection(u, u, i, j, n, proj); /* project column j of u onto column i of u */

			double *projg;
			projg = malloc(n * sizeof(*projg));
			double c;
			c = constant(u, u, i, j, n);
			multiplyg(g, i, n, c, projg);

			int x;
			for(x = 0; x < n; x++) {
				u[x * n + j] = u[x * n + j] - proj[x];
				g[x * n + j] = g[x * n + j] - projg[x];
			}
			free(proj);
			free(projg);
		}
	}
	/* u now points to orthogonalized a */

	double *tu;
	tu = malloc(n * n * sizeof(*tu));
	transpose(u, tu, n);

	multiply(g, tu, b, n);

	free(tu);
	free(g);
}

/*
	multiplies a contant onto the colg-th column of g then stores this
	in projg
*/
void multiplyg(double *g, int colg, int n, double constant, double *projg) {
	int i;
	for(i = 0; i < n; i++) {
		projg[i] = constant * g[i * n + colg];
	}
}

/*
	the constant we use in projections
*/
double constant(double *u, double *a, int colu, int cola, int n) {
	double top;
	top = innerproduct(u, a, colu, cola, n);
	double bottom;
	bottom = innerproduct(u, u, colu, colu, n);
	double c;
	c = top / bottom;

	return c;
}

/*
	multiply two matrices together
*/
void multiply(double *g, double *u, double *b, int n) {
	int i;
	int j;
	int k;
	double sum;
	sum = 0;
	for(i = 0; i < n; i++) { /* the row of g we're on */
		for(j = 0; j < n; j++) { /* the column of u we're on */
			for(k = 0; k < n; k++) { /* going down the row or column */
				sum = sum + g[i * n + k] * u[k * n  + j];
			}
			b[i * n + j] = sum;
			sum = 0;
		}
	}
}

/*
	divides the col-th column of g by a certain size
*/
void makeg(double *g, int col, int n, double sizea) {
	int i;
	for(i = 0; i < n; i++) {
		g[i * n + col] = g[i * n + col] / sizea;
	}
}

void normalize(double *a, int col, int n) {
	double fsize;
	fsize = msize(a, col, n);
	
	int i;
	for(i = 0; i < n; i++) {
		a[i * n + col] = a[i * n + col] / fsize;
	}
}

void transpose(double *u, double *newu, int n) {
	int i;
	int j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			newu[i * n + j] = u[j * n + i];
		}
	}
}

/*
	length of vector
*/
double msize(double *a, int col, int n) {
	double sqsize;
	sqsize = innerproduct(a, a, col, col, n);
	double size;
	size = sqrt(sqsize);

	return size;
}

double innerproduct(double *u, double *a, int colu, int cola, int n) {
	double ip;
	ip = 0;
	int i;
	for(i = 0; i < n; i++) {
		ip = ip + u[colu + n * i] * a[cola + n * i];
	}

	return ip;
}

/*
	projects the cola-th column of a onto the colu-th column of u
	and stores the projection in proj
*/
void projection(double *u, double *a, int colu, int cola, int n, double *proj) {
	double ip;
	ip = innerproduct(u, a, colu, cola, n);
	double bottom;
	bottom = innerproduct(u, u, colu, colu, n);

	int i;
	for(i = 0; i < n; i++) {
		proj[i] = (ip / bottom) * u[i * n + colu];
	}
}

/*
void multiplymatrices(double *in1, double *in2, double *out, int n) {
	double **rows = malloc(n * sizeof(double *));
	double **cols = malloc(n * sizeof(double *));
	int i;
	int j;

	for(i = 0; i < n; i++) {
		rows[i] = malloc(n * sizeof(double));
		for(j = 0; j < n; j++) {
			rows[i][j] = in1[i * n + j];
		}
	}

	for(j = 0; j < n; j++) {
		cols[j] = malloc(n * sizeof(double));
		for(i = 0; i < n; i++) {
			cols[j][i] = in2[i * n + j];
		}
	}

	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			out[i * n + j] = ijentry(rows[i], cols[j], n);
 		}
	}
 
}
*/



