#include <cmath>
#include "inter_app.h"

#define EPS 1e-16
#define N 4

double BuildingLagrange(const double *x, const double *f, double x0) {
	double result = 0.0;
	double sum = 0.0;
	double prod = 0.0;
	    for (int i = 0; i < 4; i++) {
        	sum = 0.0;
	        for (int j = 0; j < 4; j++) {
	            if (j == i)  {continue; }
	            prod = 1.0;
	            for (int k = 0; k < 4; k++) {
	                if (k == i || k == j)  {continue;}
	                prod *= (x0 - x[k]) / (x[i] - x[k]);
	            }
	            sum += prod / (x[i] - x[j]);
	        }
	        result += f[i] * sum;
	    }
	    return result;
}


// extra = n * f[x_i, x_i+1] + n * w_i + n * d_i
int BuildingAkima(int n, const double *x, const double *f, double *a, double *extra) {
	double w_l = 0.0;
	double w_r = 0.0;
	double* divided_diff = &extra[0];
	double *w = &extra[n];
	double *d = &extra[2*n];
        if (n < 4)  {return -1;  }
	for(int i = 0; i < n - 1; i++) {
		divided_diff[i] = (f[i+1] - f[i])/(x[i + 1] - x[i]);
	}
	for(int i = 1; i < n - 1; i++) {
		w[i] = fabs(divided_diff[i] - divided_diff[i - 1]);
	}
	for(int i = 2; i < n - 2; i++) {
		w_l = w[i-1];
		w_r = w[i+1];
		if(w_l * w_l + w_r * w_r > EPS) {
			d[i] = (w_r * divided_diff[i -1]  + w_l * divided_diff[i -1])/
										(w_r + w_l);
		}
		else {
			d[i] =  ((x[i+1] - x[i]) * divided_diff[i-1] +
					(x[i] - x[i-1]) * divided_diff[i]) / (x[i+1] - x[i-1]);
		}
	}


	double x_left[4]  = {x[0], x[1], x[2], x[3]};
	double f_left[4]  = {f[0], f[1], f[2], f[3]};
	d[0] = BuildingLagrange(x_left, f_left, x[0]);
	d[1] = BuildingLagrange(x_left, f_left, x[1]);

	double x_right[4] = {x[n-4], x[n-3], x[n-2], x[n-1]};
	double f_right[4] = {f[n-4], f[n-3], f[n-2], f[n-1]};
	d[n-2] = BuildingLagrange(x_right, f_right, x[n-2]);
	d[n-1] = BuildingLagrange(x_right, f_right, x[n-1]);

	for(int i = 0; i < n - 1; i++) {
		a[4*i] = f[i];
		a[4*i + 1] = d[i];
		a[4*i + 2] = (3*divided_diff[i] - 2*d[i] - d[i+1])/(x[i+1] - x[i]);
		a[4*i + 3] = (d[i] + d[i+1] - 2*divided_diff[i])/(pow(x[i + 1] - x[i], 2));
	}
	return 0;
}

double EvaluationAkima(double x, double a, double b, int n, const double *X, const double *A) {
	double dx = 0.0;
	int i = 0;
	int left = 0;
	int right = 0;
	int middle = 0;
	(void)a;
	(void)b;
        if (x <= X[0]) {
	        dx = x - X[0];
	        return (A[0] + A[1]*dx + A[2]*dx*dx + A[3]*dx*dx*dx);
        }
        if (x >= X[n -1]) {
	        dx = x - X[n-2];
	        i = n - 2;
	        return (A[4*i] + A[4*i+1]*dx + A[4*i+2]*dx*dx + A[4*i+3]*dx*dx*dx);
        }

         left = 0;
	 right = n - 2;
	 while (left < right) {
	         middle = (left + right) / 2;
	         if (X[middle+1] < x) {
	              left = middle + 1;
		 }
	         else {
	                right = middle;
		}
	}
	i = left;
	dx = x -X[i];
	return (A[4*i] + A[4*i+1]*dx + A[4*i+2]*dx*dx + A[4*i+3]*dx*dx*dx);
}


int BuildingSplines(int n, const double *x, const double *f, double *a, double *extra) {
	if(n < 4) { return -1; }
	double*divided_diff = &extra[0];
	double x0 = 0.0;
	double f_x0 = 0.0;
	double x_n1 = 0.0;
	double f_xn1 = 0.0;
	double div_f1 = 0.0;
	double div_f2 = 0.0;
	double*d = &extra[n];
	double*upper = &extra[2*n]; // n - 3
	double*diag = &extra[3*n - 3];// n - 2
	double*lower = &extra[4*n - 5]; // n- 3
	double*b = &extra[5*n - 8]; // n - 2
	double temp = 0.0;
	int t = 0;
	for(int i = 0; i < n - 1; i++) {
		divided_diff[i] = (f[i+1] - f[i]) / (x[i+1] - x[i]);
	}
	x0 = x[0] - (x[1] - x[0]);
	f_x0 = f[0] - (f[1] - f[0]);
	div_f1 = (f[0] - f_x0)/(x[0] - x0);
	d[0] = ((x[1] -x[0])* div_f1  +  (x[0] - x0) * divided_diff[0]) / (x[1] - x0);

	x_n1 = x[n - 1] + (x[n - 1] - x[n - 2]);
	f_xn1 = f[n - 1] + (f[n - 1] - f[n - 2]);
	div_f2 = (f_xn1 - f[n - 1])/(x_n1 - x[n -1]);
	d[n - 1] = ((x_n1 -x[n - 1])* divided_diff[n - 2]+
				(x[n -1] - x[n - 2]) * div_f2)/(x_n1 - x[n-2]);
	// have to be changed
	for(int i = 1; i < n - 1; i++) {
		t = i - 1;
		upper[t] = x[i] - x[i -1];
		b[t] = 3*divided_diff[i - 1]*(x[i + 1] - x[i]) +
				3*divided_diff[i]*(x[i] - x[i - 1]);
		if (t < n - 2) {
			lower[t] = x[i + 1] - x[i];
			diag[t] = 2*(x[i + 1] - x[i - 1]);
	        }
	}
	b[0] -= (x[2] - x[1]) * d[0];
	b[n - 3] -= (x[n-1] - x[n -2]) * d[n - 1];
	// gauss
	for(int i = 1; i < n - 2; i++) {
		temp = lower[i]/diag[i -1];
		lower[i] = 0;
		diag[i] -= temp*upper[i-1];
		b[i] -= temp*b[i-1];
	}
	d[n -2] = b[n - 3]/diag[n - 3];
	for(int i = n - 4; i >= 0; i--) {
		d[i + 1] = (b[i] - upper[i] * d[i + 2])/diag[i];
	}
	for(int i = 0; i < n - 1; i++) {
		a[4*i] = f[i];
		a[4*i + 1] = d[i];
		a[4*i + 2] = (3*divided_diff[i] - 2*d[i] - d[i+1])/(x[i+1] - x[i]);
		a[4*i + 3] = (d[i] + d[i+1] - 2*divided_diff[i])/(pow(x[i + 1] - x[i], 2));
	}
	return 0;
}

double EvaluationSplines(double x, double a, double b, int n, const double *X, const double *A) {
	double dx = 0.0;
	int i = 0;
	int left = 0;
	int right = 0;
	int middle = 0;
	(void)a;
	(void)b;
        if (x <= X[0]) {
	        dx = x - X[0];
	        return (A[0] + A[1]*dx + A[2]*dx*dx + A[3]*dx*dx*dx);
        }
        if (x >= X[n -1]) {
	        dx = x - X[n-2];
	        i = n - 2;
	        return (A[4*i] + A[4*i+1]*dx + A[4*i+2]*dx*dx + A[4*i+3]*dx*dx*dx);
        }


         left = 0;
	 right = n - 2;
	 while (left < right) {
         middle = (left + right) / 2;
         if (X[middle+1] < x) {
              left = middle + 1;
	 }
         else {
                right = middle;
		}
	}
	i = left;
	dx = x -X[i];
	return (A[4*i] + A[4*i+1]*dx + A[4*i+2]*dx*dx + A[4*i+3]*dx*dx*dx);
}
