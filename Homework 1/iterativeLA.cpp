// iterativeLA.cpp
// Matt Murray

#include <iostream>
using namespace std;

#include <math.h>

#include "iterativeLA.h"
#include "matrix.h"

#define MONITOR


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
#ifdef MONITOR
        cout << "x: " << x <<" , Iter: " << iter+1 << " , l2-error: " << l2error << endl;
#endif
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
#ifdef MONITOR
        cout << "x: " << x <<" , Iter: " << iter+1 << " , l2-error: " << l2error << endl;
#endif
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
#ifdef MONITOR
        cout << "w: " << w << " , x: " << x <<" , Iter: " << iter+1 << " , l2-error: " << l2error << endl;
#endif
        if( l2error <= tol) {
            maxIter = iter+1;
            return SUCCESS;
        }
        y = x;
    }
    
    return WONT_STOP;
}
