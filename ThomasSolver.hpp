#ifndef THOMAS_SOLVER_HPP
#define THOMAS_SOLVER_HPP

#include <vector>

// ---------------------------------------------------------------------------
//  Tridiagonal solve in O(n). Pure linear algebra: if the word "option" ever
//  appears in this file, something has gone wrong.
//
//      | b0 c0             | |x0|   |d0|
//      | a1 b1 c1          | |x1| = |d1|
//      |    a2 b2 c2       | |..|   |..|
//      |          a_{n-1} b_{n-1} | |x_{n-1}|   |d_{n-1}|
//
//  a[0] and c[n-1] are ignored. Both b and d are OVERWRITTEN; on return d
//  holds the solution.
//
//  No pivoting. Stable when the matrix is diagonally dominant, which the
//  BS discretisation satisfies as long as the cell Peclet number stays below
//  2, i.e. |b(S)| * h < 2 * a(S). Violate that and you need upwinding.
// ---------------------------------------------------------------------------
namespace thomas {

    void SolveInPlace(const std::vector<double>& a,
        std::vector<double>& b,
        const std::vector<double>& c,
        std::vector<double>& d);

}  // namespace thomas

#endif // THOMAS_SOLVER_HPP