#ifndef GRID_HPP
#define GRID_HPP

#include <cstddef>
#include <vector>

#include "Domain.hpp"

// ---------------------------------------------------------------------------
//  Spatial stencil coefficients at an interior node i.
//
//  On a NON-uniform grid the standard uniform formulas are wrong. With
//      hm = S_i - S_{i-1}      (backward spacing)
//      hp = S_{i+1} - S_i      (forward  spacing)
//  the three-point approximations are
//
//      V_S  (S_i) ~ dm * V_{i-1} + d0 * V_i + dp * V_{i+1}
//      V_SS (S_i) ~ em * V_{i-1} + e0 * V_i + ep * V_{i+1}
//
//  which reduce to the familiar 1/(2h), 1/h^2 forms when hm == hp.
// ---------------------------------------------------------------------------
struct Stencil
{
    // first derivative
    double dm{ 0.0 };
    double d0{ 0.0 };
    double dp{ 0.0 };
    // second derivative
    double em{ 0.0 };
    double e0{ 0.0 };
    double ep{ 0.0 };
};

// ---------------------------------------------------------------------------
//  Abstract spatial/temporal discretisation.
//
//  Knows only about the solution domain and how many nodes to place in it.
//  It has no idea what PDE is being solved or what the payoff is.
// ---------------------------------------------------------------------------
class Grid
{
public:
    virtual ~Grid() = default;

    // Spatial nodes, strictly increasing: S_[0] == Smin, S_[M] == Smax
    const std::vector<double>& S() const { return S_; }

    // Time nodes, strictly increasing: t_[0] == 0, t_[N] == T
    const std::vector<double>& t() const { return t_; }

    // Number of intervals (number of nodes = number of intervals + 1)
    std::size_t M() const { return S_.size() - 1; }
    std::size_t N() const { return t_.size() - 1; }

    // Spacings around interior node i
    double hMinus(std::size_t i) const { return S_[i] - S_[i - 1]; }
    double hPlus(std::size_t i) const { return S_[i + 1] - S_[i]; }

    // Time step n -> n+1
    double dt(std::size_t n) const { return t_[n + 1] - t_[n]; }

    // Finite-difference weights at interior node i (1 <= i <= M-1)
    Stencil StencilAt(std::size_t i) const;

    // Index of the node closest to S. Used by Price() before interpolating.
    std::size_t NearestIndex(double S) const;

    // Linear interpolation of a solution vector V (size M+1) at S
    double Interpolate(const std::vector<double>& V, double S) const;

    // True if some node sits on S to within tol. Barrier grids must satisfy
    // this at S == B, otherwise convergence drops from O(h^2) to O(h).
    bool HasNodeAt(double S, double tol = 1e-10) const;

protected:
    // Derived classes fill these in their constructor, then call Validate().
    std::vector<double> S_;
    std::vector<double> t_;

    void BuildUniformTimeNodes(double T, std::size_t N);
    void Validate() const;   // throws if nodes are not strictly increasing
};



// ---------------------------------------------------------------------------
//  Uniform grid. The obvious baseline; use it as the reference in the
//  convergence study and to sanity-check anything fancier.
// ---------------------------------------------------------------------------
class UniformGrid final : public Grid
{
public:
    UniformGrid(const Domain& domain, std::size_t M, std::size_t N);
};

// ---------------------------------------------------------------------------
//  Tavella-Randall sinh grid: clusters nodes around `focus`.
//
//      S(xi) = focus + density * sinh( c1 + xi * (c2 - c1) ),   xi in [0,1]
//      c1 = asinh( (Smin - focus) / density )
//      c2 = asinh( (Smax - focus) / density )
//
//  `density` controls the clustering:
//      small (e.g. 0.05 * focus)  ->  aggressive clustering
//      large (e.g. 10  * focus)   ->  effectively uniform
//
//  The xi-grid is built piecewise-uniform so that `focus` lands EXACTLY on a
//  node. For a vanilla option pass focus = K (resolves the payoff kink); for a
//  barrier pass focus = B (the barrier MUST be a node).
//
//  Ref: Tavella & Randall, "Pricing Financial Instruments: The Finite
//       Difference Method" (2000), Ch. 5.
// ---------------------------------------------------------------------------
class SinhGrid final : public Grid
{
private:
    double focus_;
    double density_;

public:
    SinhGrid(const Domain& domain,
        std::size_t M,
        std::size_t N,
        double focus,
        double density);

    double Focus()   const { return focus_; }
    double Density() const { return density_; }
};

#endif // GRID_HPP