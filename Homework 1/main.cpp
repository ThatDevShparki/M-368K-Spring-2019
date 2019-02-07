// main.cpp
// Matt Murray

#include <iostream>
using namespace std;

#include "iterativeLA.h"

int main() {
    
    int n;
    cout << "Enter size: " << flush;
    cin >> n;
    
    Matrix A(n,n);
    Vector x(n);
    Vector b(n);
    
    x = 0;
    
    cout << "Enter A by rows: " << flush;
    cin >> A;
    
    cout << "Enter b: " << flush;
    cin >> b;
    
    int maxIter;
    double tolerance;
    cout << "Enter maxIter and tolerance: " << flush;
    cin >> maxIter >> tolerance;
    
    int method;
    cout << "Which method would you like to use to solve?" << endl;
    cout << "   (1) Jacobi Iteration" << endl;
    cout << "   (2) Gauss-Seidel Iteration" << endl;
    cout << "   (3) Successive Over-Relaxation Iteration" << endl;
    cout << "Enter method number:" << endl;
    cin >> method;
    
    double w;
    state s = BAD_METHOD;
    switch(method) {
        // Jacobi Iteration
        case 1:
            s = jacobi(A,b,x,maxIter,tolerance);
            break;
        case 2:
            s = gauss_seidel(A,b,x,maxIter,tolerance);
            break;
        case 3:
            cout << "Enter omega:" << endl;
            cin >> w;
            s = sor(A, b, x, w, maxIter, tolerance);
            break;
    }
    
    switch(s) {
        case WONT_STOP:
            cout << "ERROR: Exceeded maximum number of iterations." << endl;
            return 1;
        case BAD_SDD:
            cout << "ERROR: Expected Strictly Diagonally Dominant Matrix." << endl;
            return 1;
        case BAD_DIAGONAL:
            cout << "ERROR: A diagonal entry of A was 0." << endl;
            return 1;
        case BAD_METHOD:
            cout << "ERROR: Unknown method. Expected values 1, 2, or 3." << endl;
            return 1;
        case BAD_W:
            cout << "ERROR: Omega (w) must be between 0 and 2 exclusively." << endl;
            return 1;
        default:s
            cout << "ERROR: Unspecified." << endl;
            return 1;
        case SUCCESS:
            cout << "The solution is:" << endl;
            cout << x << endl;
            
            Vector y(n);
            matVecMult(A,x,y);
            y -= b;
            cout << "The number of iterations is: " << maxIter << endl;
            cout << "The max-norm of residual is: " << maxNorm(y) << endl;
            cout << "The residual is: " << endl;
            cout << y << endl;
            return 0;
    }
}

