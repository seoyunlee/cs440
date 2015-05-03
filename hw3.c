/*
Author: Seoyun Lee						Class: CPSC 440
Homework Assignment: 3

Finds the eigenvalues of a matrix given a symmetric tridiagonal
form.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void qr_symmetric(double *a, int n, double *b);
void copy(double *a, int n, double *b);
void iteration(double *b, int n);
void addback(double *b, int n, double mu);
void shift(double *b, int n, double mu);
void givens(double *b, int n, int j, double *q);
void mult1(double *b, double *q, int j, int n, double *bcopy);
void transpose(double *q, double *qt, int n, int j);
void mult2(double *b, double *qt, int j, int n, double *bcopy);
void submatrix(double *b, int n, int m, double *bnew);

int main() {
	return 0;
}

void qr_symmetric(double *a, int n, double *b) {
	copy(a, n, b);
	int i;
	for(i = n; i > 1; i--) {
		double *bnew;
		bnew = malloc(i * i * sizeof(double));
		
		submatrix(b, n, i, bnew);

		while(bnew[1] != 0.000) {
			iteration(bnew, i);
		}

		int j;
		int k;
		for(j = n - i; j < n; j++) {
			for(k = n - i; k < n; k++) {
				b[j * n + k] = bnew[(j - (n - i)) * i + (k - (n - i))];
			}
		}
		free(bnew);
	}
}

/*
	takes a pointer to a matrix a of size n*n and
	copies it to b
*/
void copy(double *a, int n, double *b) {
	int i;
	for(i = 0; i < n*n; i++) {
		b[i] = a[i];
	}
}

void iteration(double *b, int n) {
	double *track;
	track = malloc(n * n * sizeof(double));
	copy(b, n, track);

	double mu;
	mu = b[0];

	shift(b, n, mu);
	shift(track, n, mu);

	int j;
	for(j = n - 2; j >= 0; j--) {
		double *q;
		q = malloc(n * n * sizeof(double));
		givens(b, n, j, q);
		double *bcopy;
		bcopy = malloc(2 * n * sizeof(double));
		mult1(b, q, j, n, bcopy);
		free(bcopy);
		free(q);
	}

	for(j = n - 2; j >= 0; j--) {
		double *q;
		q = malloc(n * n * sizeof(double));
		givens(track, n, j, q);

		double *trackcopy;
		trackcopy = malloc(2 * n * sizeof(double));
		mult1(track, q, j, n, trackcopy);

		double *qt;
		qt = malloc(n * n * sizeof(double));
		transpose(q, qt, n, j);

		double *bcopy;
		bcopy = malloc(2 * n * sizeof(double));
		mult2(b, qt, j, n, bcopy);

		free(trackcopy);
		free(bcopy);
		free(qt);
		free(q);
	} 

	free(track);

	addback(b, n, mu);
}

/*
	takes the value of mu that we subtracted before and adds it to the
	diagonal of b
*/
void addback(double *b, int n, double mu) {
	int i;
	int j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			if(i == j) {
				b[i * n + j] = b [i * n + j] + mu;
			}
		}
	}
}

/*
	takes the value at b[0] and subtracts it from the diagonal of b
*/
void shift(double *b, int n, double mu) {
	int i;
	int j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			if(i == j) {
				b[i * n + j] = b[i * n + j] - mu;
			}
		}
	}
}

void givens(double *b, int n, int j, double *q) {
	double angle;
	double x;
	x = b[n * j + j + 1];
	double y;
	y = b[n * (j + 1) + j + 1];
	angle = atan(x / y);
	q[j * n + j] = cos(angle);
	q[j * n + j + 1] = -sin(angle);
	q[(j + 1) * n + j] = sin(angle);
	q[(j + 1) * n + j + 1] = cos(angle);
}

void mult1(double *b, double *q, int j, int n, double *bcopy) {
	int k;
	for(k = 0; k < n; k++) {
		bcopy[k] = q[j * n + j] * b[j * n + k] 
					+ q[j * n + j + 1] * b[(j + 1) * n + k];
		bcopy[n + k] =q[(j + 1) * n + j] * b[j * n + k] 
						+ q[(j + 1) * n + j + 1] * b[(j + 1) * n + k];
	}

	for(k = 0; k < n; k++) {
		b[j * n + k] = bcopy[k];
		b[(j + 1) * n + k] = bcopy[n + k];
	}
}

void transpose(double *q, double *qt, int n, int j) {
	qt[j * n + j] = q[j * n + j];
	qt[j * n + j + 1] = q[(j + 1) * n + j];
	qt[(j + 1) * n + j] = q[j * n + j + 1];
	qt[(j + 1) * n + j + 1] = q[(j + 1) * n + j + 1];
}

void mult2(double *b, double *qt, int j, int n, double *bcopy) {
	int i;
	for(i = 0; i < n; i++) {
		bcopy[i] = b[i * n + j] * qt[j * n + j]
					+ b[i * n + j + 1] * qt[(j + 1) * n + j];
		bcopy[n + i] = b[i * n + j] * qt[j * n + j + 1]
						+ b[i * n + j + 1] * qt[(j + 1) * n + j + 1];
	}

	for(i = 0; i < n; i++) {
		b[i * n + j] = bcopy[i];
		b[i * n + j + 1] = bcopy[n + i];
	}
}

void submatrix(double *b, int n, int m, double *bnew) {
	int i;
	int j;
	for(i = 0; i < m; i++) {
		for(j = 0; j < m; j++) {
			bnew[i * m + j] = b[(n - m + i) * n + (n - m + j)];
		}
	}
}
