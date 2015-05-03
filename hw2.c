/*
Author: Seoyun Lee						Class: CPSC 440
Homework Assignment: 2

Uses Givens rotations to take a matrix A and transform it into an
upper Hessenberg matrix U*A*U^T = B. Uses simplified multiplication in
order to speed up the computation and keep it in O(n^3) time
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void upperhes(double *a, int n, double *u, double *b);
double angle(double x, double y);
void makeux(double *ux, int first, int second, int n, double angle);
void copy(double *u, double *copy, int n);
void transpose(double *ux, double *uxt, int n);
double ijentry(double *rowi, double *colj, int n);
void multiplyrows(double *in1, int first, int second, double *in2, int n);
void multiplycols(double *in1, int first, int second, double *in2, int n);

int main() {
	int n;
	n = 2;
	double *a;
	a = malloc(n * n * sizeof(double));
	a[0] = 3;
	a[1] = 2;
	a[2] = 2;
	a[3] = 3;
	double *b;
	b = malloc(n * n * sizeof(double));
	double *u;
	u = malloc(n * n * sizeof(double));
	upperhes(a, n, u, b);
	int i;
	for(i = 0; i < n*n; i++) {
		printf("%f ", b[i]);
	}
}

/*
	*a: points to an n*n size array with entries of a matrix a(1,1) ... a(n,n)
	n: size of the matrix (n x n)
	*u: points to the array of doubles size n*n, and is orthogonal, created
		by composition of Givens rotations (memory allocated by user)
	*b: points to the array of doubles size n*n, is the upper Hessenberg
		matrix (memory allocated by user)
*/
void upperhes(double *a, int n, double *u, double *b) {
	copy(a, b, n); /* *b now points to n*n entry corresponding to *a */
	
	int i;
	int j;

	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			if(i == j)
				u[i * n + j] = 1;
			else
				u[i * n + j] = 0;
		}
	}

	for(j = 0; j < n - 2; j++) {
		for(i = n - 1; i > j + 1; i--) {
			double ang;
			ang = angle(b[n * (i - 1) + j], b[n * i + j]);
			double *ux;
			ux = malloc(n * n * sizeof(double));
			makeux(ux, i - 1, i, n, ang);
			multiplyrows(ux, i - 1, i, u, n);
			double *ut;
			ut = malloc(n * n * sizeof(double));
			transpose(ux, ut, n);
			multiplyrows(ux, i - 1, i, b, n);
			multiplycols(b, i - 1, i, ut, n);
			free(ux);
			free(ut);
		}
	}
}

/*
	x: the value of the element above the one we want to turn 0
	y: the value we want to make 0
	returns: the angle needed for the Givens rotation to make y 0
*/
double angle(double x, double y) {
	double ang;
	ang = atan(- y / x);
	return ang;
}

/*
	*ux: points to an array of size n*n, contains the Givens rotation
			needed to change the element in the second row to 0
	first: the row above the row that gets changed to 0
	second: the row that gets changed to 0, used to place the sin, cos
	n: size of matrix
	angle: the angle of rotation
*/
void makeux(double *ux, int first, int second, int n, double angle) {
	int i; /* rows */
	int j; /* columns */
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			if(i == first && j == first)
				ux[i * n + j] = cos(angle);
			else if(i == first && j == second)
				ux[i * n + j] = -sin(angle);
			else if(i == second && j == first)
				ux[i * n + j] = sin(angle);
			else if(i == second && j == second)
				ux[i * n + j] = cos(angle);
			else if(i == j)
				ux[i * n + j] = 1;
			else
				ux[i * n + j] = 0;
		}
	}
}

/*
	*old: points to an array of size n*n with entries of a matrix
	*new: points to an array of size n*n, memory allocated by calling
			function, contains the same thing as old after calling
	n: size of matrix
*/
void copy(double *old, double *new, int n) {
	int i;
	for(i = 0; i < n * n; i++) {
		new[i] = old[i];
	}
}

/*
	*u: pointer to the matrix u
	*ut: pointer to the transpose of u
	n: size of matrix
*/
void transpose(double *u, double *ut, int n) {
	int i;
	int j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			ut[i * n + j] = u[j * n + i];
		}
	}
}

/*
	*rowi: pointer to the ith row from a matrix, length n
	*colj: pointer to the jth column from a matrix, lengthn
	n: size of the square matrix
	returns: the ij-th entry of the new matrix after multiplication
*/
double ijentry(double *rowi, double *colj, int n) {
	double ans;
	ans = 0;
	int i;
	for(i = 0; i < n; i++) {
		ans = ans + rowi[i] * colj[i];
	}
	return ans;
}

/*
	*in1: pointer to the first matrix
	first: the number of the first row that is to be multiplied
	second: the number of the second row that is to be multiplied
	*in2: pointer to the second matrix that later contains in1*in2
	n: size of matrix (n*n)
*/
void multiplyrows(double *in1, int first, int second, double *in2, int n) {
	int i;
	int j;
	double ** cols = malloc(n * sizeof(double *));
	for(j = 0; j < n; j++) {
		cols[j] = malloc(n * sizeof(double));
		for(i = 0; i < n; i++) {
			cols[j][i] = in2[i * n + j];
		}
	}

	double *row1;
	row1 = malloc(n * sizeof(double));
	double *row2;
	row2 = malloc(n * sizeof(double));
	for(i = 0; i < n; i++) {
		row1[i] = in1[first * n + i];
		row2[i] = in1[second * n + i];
	}
	
	for(j = 0; j < n; j++) {
		in2[first * n + j] = ijentry(row1, cols[j], n);
		in2[second * n + j] = ijentry(row2, cols[j], n);
		free(cols[j]);
	}
	free(row2);
	free(row1);
	free(cols);
}

/*
	*in1: pointer to the first matrix that later contains in1*in2
	first: the number of the first column that is to be multiplied
	second: the number of the second column that is to be multiplied
	*in2: pointer to the second matrix
	n: size of matrix (n*n)
*/
void multiplycols(double *in1, int first, int second, double *in2, int n) {
	double ** rows = malloc(n * sizeof(double *));
	int i;
	int j;
	for(i = 0; i < n; i++) {
		rows[i] = malloc(n * sizeof(double));
		for(j = 0; j < n; j++) {
			rows[i][j] = in1[i * n + j];
		}
	}

	double *col1;
	col1 = malloc(n * sizeof(double));
	double *col2;
	col2 = malloc(n * sizeof(double));

	for(i = 0; i < n; i++) {
		col1[i] = in2[i * n + first];
		col2[i] = in2[i * n + second];
	}

	for(i = 0; i < n; i++) {
		in1[i * n + first] = ijentry(rows[i], col1, n);
		in1[i * n + second] = ijentry(rows[i], col2, n);
		free(rows[i]);
	}

	free(col2);
	free(col1);
	free(rows);
}

