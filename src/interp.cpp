#include "interp.h"

/* idx = getNearGLNodeIdx(x, m)
* Get the index of the Gauss-Legendre node of order m
* nearest and less than the point x on the interval [a, b]
* x : evaluation point
* m : order of near nodes (less than order of x)
* a : lower bound of interval
* b : upper bound of interval
* idx : Index of nearest/less Gauss-Legendre node
*/
int Interp::getNearGLNodeIdx(
    const double xi, const int m, const double a = -1.0, const double b = 1.0) {

    const double leng = b - a;
    const double mid = (a + b)/2.0;
    assert(leng > 0);

    const double x = 2.0*(xi - mid) / leng; // Change to interval [-1,1]

    const int idx = m - floor(((4.0*m+2.0) * acos(x) / PI + 1.0)/4.0) - 1;

    assert(idx >= -1 && idx < m);

    return idx;
}

/* evalLagrangeBasis(x, xs, k)
* Evaluate the Lagrange basis polynomial taking on 1 at xs[k]
* and 0 at xs[j] for j \neq k, at the point x
* x  : evaluation point
* xs : interpolation nodes
* k  : index of basis function \in {0,1,...,order}
*/
double Interp::evalLagrangeBasis(
    const double x, const realVec& xs, const int k) {

    // assert(k < xs.size());

    double product = 1.0;

    for (int j = 0; j < xs.size(); ++j) {
        if (j == k) continue;

        product *= (x - xs[j]) / (xs[k] - xs[j]);
    }

    return product;
}

/* evalTrigBasis(x, xs, k)
* Evaluate the trigonometric basis function taking on 1 at xs[k]
* and 0 at xs[j] for j \neq k, at the point x
* x  : evaluation point
* xs : interpolation nodes
* k  : index of basis function \in {0,1,...,N-1}
*/
double Interp::evalTrigBasis(
    const double x, const realVec& xs, const int k) {

    const int N = xs.size();
    assert(k <= N-1);

    const double diff = x - xs[k];
    const double denom =
        N%2 ? N*sin(diff/2.0) : N*tan(diff/2.0);

    return denom ?
        sin(N*diff/2.0) / denom :
        1.0;
}