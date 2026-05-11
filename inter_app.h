#ifndef INTER_APP_H
#define INTER_APP_H

int BuildingAkima(int n, const double *x, const double *f, double *a, double *extra);

double EvaluationAkima(double x, double a, double b, int n, const double *X, const double *A);

double BuildingLagrange(const double *x, const double *f, double x0);

int BuildingSplines(int n, const double *x, const double *f, double *a, double *extra);

double EvaluationSplines(double x, double a, double b, int n, const double *X, const double *A);


#endif
