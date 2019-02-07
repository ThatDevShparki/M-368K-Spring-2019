// iterativeLA.h
// Matt Murray

#ifndef IterativeLA_Included
#define IterativeLA_Included

#include "matrix.h"

enum state {
    SUCCESS,
    WONT_STOP,
    BAD_SDD,
    BAD_DIAGONAL,
    BAD_DATA,
    BAD_W,
    BAD_METHOD,
    NOT_SYMMETRIC
};


void monitor_headers(const int xn);
void monitor_log(const int xn, const int iter, const Vector x, const double l2error);

bool is_sdd(const Matrix& A, int& n);

state jacobi(const Matrix& A, const Vector& b, Vector& x, int& maxIter, double tol);
state gauss_seidel(const Matrix& A, const Vector& b, Vector& x, int& maxIter, double tol);
state sor(const Matrix& A, const Vector& b, Vector& x, const double& w, int& maxIter, double tol);

state conj_grad(const Matrix& A, const Vector&b, Vector& x, int& maxIter, double tol);
state jacobi_precg(const Matrix& A, const Vector&b, Vector& x, int& maxIter, double tol);
state ssor_precg(const Matrix& A, const Vector&b, Vector& x, const double w, int& maxIter, double tol);

#endif
