// iterativeLA.cpp
// Matt Murray

#include <iostream>
#include <iomanip>
using namespace std;

#include <math.h>

#include "iterativeLA.h"
#include "matrix.h"

 #define MONITOR

void monitor_headers(int xn) {
#ifdef MONITOR
    cout << "\n" << endl;
    cout << setw(7) << left << "Iter."
         << setw(10*xn) << left << "X"
         << setw(20) << left << "l2 Error" << endl;
    
    string xnbar = "=========";
    for (int i = 0; i < xn-1; i++){
        xnbar += "==========";
    }
    cout << setw(7) << left << "======"
         << setw(10*xn) << left << xnbar
         << setw(20) << left << "===================" << endl;
#endif
}

void monitor_log(const int xn, const int iter, const Vector x, const double l2error) {
#ifdef MONITOR
    cout << setw(7) << left << iter
         << setw(10*xn) << left << x
         << setw(20) << left << l2error << endl;
//    cout << "Iter.: " << iter << ", " << flush;
//    cout << "X: " << x << ", " << flush;
//    cout << "l2 Error: " << l2error << endl;
#endif
}


// Returns true if matrix is Strictly Diagonally Dominant and false otherwise
bool is_sdd(const Matrix& A, int& n) {
    bool sdd = true;
    for (int i=0; i<n; i++) {
        int diag = fabs(A(i,i));
        int row_sum = 0;
        for (int j=0; j<n; j++){
            if (j==i) continue;
            row_sum += fabs(A(i,j));
        }
        if (row_sum >= diag) {
            sdd = false;
        }
    }
    return sdd;
}

// Jacobi Iteration Algorithm
state jacobi(const Matrix& A, const Vector& b, Vector& x, int& maxIter, double tol) {
    
    // Check data
    int n = A.n(0);
    if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA;
    if(tol <= 0) return BAD_DATA;
    if(maxIter <= 0) maxIter = 1;
    
    for(int i=0; i<n; i++) {
        if(A(i,i) == 0) return BAD_DIAGONAL;
    }
    
    // Check for Strictly Diagonally Dominant
    if (!is_sdd(A, n)) return BAD_SDD;
    
    // Apply Jacobi Iteration
    Vector y(x);
    
    for(int iter=0; iter<maxIter; iter++) {
        
        // Get new x
        for(int i=0; i<n; i++) {
            double sum = 0;
            for(int j=0; j<n; j++) {
                if(j==i) continue;
                sum += A(i,j)*y(j);
            }
            x(i) = ( -sum + b(i) ) / A(i,i);
        }
        
        // Check error tolerance
        y -= x;
        double l2error = l2norm(y) / (l2norm(x)+1e-16);
        
        monitor_log(n, iter+1, x, l2error);
        
        if( l2error <= tol) {
            maxIter = iter+1;
            return SUCCESS;
        }
        y = x;
    }
    
    return WONT_STOP;
}

// Gauss-Seidel Iteration Algorithm
state gauss_seidel(const Matrix& A, const Vector& b, Vector& x, int& maxIter, double tol) {
    
    // Check data
    int n = A.n(0);
    if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA;
    if(tol <= 0) return BAD_DATA;
    if(maxIter <= 0) maxIter = 1;
    
    for(int i=0; i<n; i++) {
        if(A(i,i) == 0) return BAD_DIAGONAL;
    }
    
    // Check for Strictly Diagonally Dominant
    if (!is_sdd(A, n)) return BAD_SDD;
    
    // Apply Gauss-Seidel Iteration
    Vector y(x);
    
    for(int iter=0; iter<maxIter; iter++) {
        
        // Get new x
        for(int i=0; i<n; i++) {
            double sum_1 = 0;
            for(int j=0; j<i; j++) {
                sum_1 += A(i,j)*x(j);
            }
            double sum_2 = 0;
            for(int j=i+1; j<n; j++) {
                sum_2 += A(i,j)*y(j);
            }
            x(i) = (b(i) - sum_1 - sum_2) / A(i,i);
        }
        
        // Check error tolerance
        y -= x;
        double l2error = l2norm(y) / (l2norm(x)+1e-16);
        
        monitor_log(n, iter+1, x, l2error);
        
        if( l2error <= tol) {
            maxIter = iter+1;
            return SUCCESS;
        }
        y = x;
    }
    
    return WONT_STOP;
}

// Successive Over-Relaxation Iteration Algorithm
state sor(const Matrix& A, const Vector& b, Vector& x, const double& w, int& maxIter, double tol) {
    
    // Check data
    int n = A.n(0);
    if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA;
    if(tol <= 0) return BAD_DATA;
    if(maxIter <= 0) maxIter = 1;
    
    for(int i=0; i<n; i++) {
        if(A(i,i) == 0) return BAD_DIAGONAL;
    }
    
    // Check for Strictly Diagonally Dominant
    if (!is_sdd(A, n)) return BAD_SDD;
    
    // Check w between 0 and 2
    if (w <= 0 || w >= 2) return BAD_W;
    
    // Apply Gauss-Seidel Iteration
    Vector y(x);
    
    for(int iter=0; iter<maxIter; iter++) {
        
        // Get new x
        for(int i=0; i<n; i++) {
            double sum_1 = 0;
            for(int j=0; j<i; j++) {
                sum_1 += A(i,j)*x(j);
            }
            double sum_2 = 0;
            for(int j=i+1; j<n; j++) {
                sum_2 += A(i,j)*y(j);
            }
            x(i) = (1-w)*y(i) + w*(b(i) - sum_1 - sum_2)/A(i,i);
        }
        
        // Check error tolerance
        y -= x;
        double l2error = l2norm(y) / (l2norm(x)+1e-16);

        monitor_log(n, iter+1, x, l2error);
        
        if( l2error <= tol) {
            maxIter = iter+1;
            return SUCCESS;
        }
        y = x;
    }
    
    return WONT_STOP;
}


/* =================================================== */

// (Unconditioned) Conjugate Gradient Method
state conj_grad(const Matrix& A, const Vector& b, Vector& x, int& maxIter, double tol) {
    
    // CHECK DATA
    
    int n = A.n(0);
    
    if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA;
    if(tol <= 0) return BAD_DATA;
    if(maxIter <= 0) maxIter = 1;
    
    for(int i=0; i<n-1; i++) {
        for(int j=i+1; j<n; j++) {
            if(A(i,j) != A(j,i) ) return NOT_SYMMETRIC;
        }
    }
    
    // INITIALIZE CONJUGATE GRADIENTS
    
    // Set initial residual r = b - Ax
    Vector r(n);
    
    matVecMult(A,x,r);
    r-=b; r*=(-1);
    
    double alpha = scDot(r,r);
    
    // Set initial search direction d = r
    Vector d(r); // Creates d and sets d = r
    
    double tolSq = tol*tol;
    
    // CONJUGATE GRADIENT LOOP
    
    for(int iter=0; iter<maxIter; iter++) {
        
        if(scDot(d,d) <= tolSq) {
            maxIter = iter+1;
            return SUCCESS;
        }
        
        // Set u = Ad
        Vector u(n);
        
        matVecMult(A,d,u);
        
        // Update x = x + td and r = r - tu
        double t = alpha / scDot(d,u);
        
        for(int i=0; i<n; i++) {
            x(i) += t*d(i);
            r(i) -= t*u(i);
        }
        
        // Get new search direction d = r + s*d;
        double beta = scDot(r,r);
        
        monitor_log(n, iter+1, x, sqrt(beta));
        
        if(beta <= tolSq) {
            maxIter = iter+1;
            return SUCCESS;
        }
        
        double s = beta / alpha;
        
        for(int i=0; i<n; i++) {
            d(i) = r(i) + s*d(i);
        }
        
        alpha = beta;
    }
    
    return WONT_STOP;
}


// Jacbi Pre-Conditioned Conjugate Gradient Method
state jacobi_precg(const Matrix& A, const Vector& b, Vector& x, int& maxIter, double tol) {
    
    // CHECK DATA
    
    int n = A.n(0);
    
    if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA;
    if(tol <= 0) return BAD_DATA;
    if(maxIter <= 0) maxIter = 1;
    
    for(int i=0; i<n-1; i++) {
        for(int j=i+1; j<n; j++) {
            if(A(i,j) != A(j,i) ) return NOT_SYMMETRIC;
        }
    }
    
    // INITIALIZE CONJUGATE GRADIENTS
    
    // Set initial residual r = b - Ax
    Vector r(n);
    matVecMult(A,x,r);
    r-=b; r*=(-1);

    Vector z(n);
    for (int i = 0; i < n; i++){
        z(i) = r(i)/A(i, i);
    }
    
    double alpha = scDot(r,z);
    
    // Set initial search direction d = r
    Vector d(z); // Creates d and sets d = r
    
    double tolSq = tol*tol;
    
    // CONJUGATE GRADIENT LOOP
    
    for(int iter=0; iter<maxIter; iter++) {
        
        if(scDot(d,d) <= tolSq) {
            maxIter = iter+1;
            return SUCCESS;
        }
        
        // Set u = Ad
        Vector u(n);
        matVecMult(A,d,u);
        
        // Update x = x + td and r = r - tu
        double t = alpha / scDot(d,u);
        
        for(int i=0; i<n; i++) {
            x(i) += t*d(i);
            r(i) -= t*u(i);
        }
        for(int i = 0; i < n; i++){
            z(i) = r(i)/A(i,i);
        }
        
        // Get new search direction d = r + s*d;
        double beta = scDot(r,z);
        
        monitor_log(n, iter+1, x, sqrt(beta));
        
        if(beta <= tolSq) {
            maxIter = iter+1;
            return SUCCESS;
        }
        
        double s = beta / alpha;
        
        for(int i=0; i<n; i++) {
            d(i) = z(i) + s*d(i);
        }
        
        alpha = beta;
    }
    
    return WONT_STOP;
}


// Symmetric SOR Pre-Conditioned Conjugate Gradient Method
state ssor_precg(const Matrix& A, const Vector& b, Vector& x, const double w, int& maxIter, double tol) {
    // CHECK DATA
    
    int n = A.n(0);
    
    if(A.n(1) != n || b.n() != n || x.n() != n) return BAD_DATA;
    if(tol <= 0) return BAD_DATA;
    if(maxIter <= 0) maxIter = 1;
    
    for(int i=0; i<n-1; i++) {
        for(int j=i+1; j<n; j++) {
            if(A(i,j) != A(j,i) ) return NOT_SYMMETRIC;
        }
    }
    
    if (w <= 0 || w >= 2) return BAD_W;
    
    // INITIALIZE CONJUGATE GRADIENTS
    
    // Set initial residual r = b - Ax
    Vector r(n);
    matVecMult(A,x,r);
    r-=b; r*=(-1);
    
    Vector c(n);
    for (int i = 0; i < n; i++){
        c(i) = r(i);
        for (int j = 0; j < i; j++){
            c(i) -= w*c(j)*(A(i, j)/A(j, j));
        }
    }
    
    Vector z(n);
    for (int i = n-1; i >= 0; i--){
        z(i) = c(i)/A(i,i);
        for (int j = n-1; j > i; j--){
            z(i) -= w*z(j)*(A(i, j)/A(j, j));
        }
    }
    
    double alpha = scDot(r,z);
    
    // Set initial search direction d = r
    Vector d(z); // Creates d and sets d = r
    
    double tolSq = tol*tol;
    
    // CONJUGATE GRADIENT LOOP
    
    for(int iter=0; iter<maxIter; iter++) {
        
        if(scDot(d,d) <= tolSq) {
            maxIter = iter+1;
            return SUCCESS;
        }
        
        // Set u = Ad
        Vector u(n);
        matVecMult(A,d,u);
        
        // Update x = x + td and r = r - tu
        double t = alpha / scDot(d,u);
        
        for(int i=0; i<n; i++) {
            x(i) += t*d(i);
            r(i) -= t*u(i);
        }
        for (int i = 0; i < n; i++){
            c(i) = r(i);
            for (int j = 0; j < i; j++){
                c(i) -= w*c(j)*(A(i, j)/A(j, j));
            }
        }
        for (int i = n-1; i >= 0; i--) {
            z(i) = c(i)/A(i,i);
            for (int j = n-1; j > i; j--){
                z(i) -= w*z(j)*(A(i, j)/A(j, j));
            }
        }
        
        // Get new search direction d = r + s*d;
        double beta = scDot(r,z);
        
        monitor_log(n, iter+1, x, sqrt(beta));
        
        if(beta <= tolSq) {
            maxIter = iter+1;
            return SUCCESS;
        }
        
        double s = beta / alpha;
        
        for(int i=0; i<n; i++) {
            d(i) = z(i) + s*d(i);
        }
        
        alpha = beta;
    }
    
    return WONT_STOP;
}
