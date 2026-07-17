#include "Grid.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// ===========================================================================
//  Grid
// ===========================================================================

Stencil Grid::StencilAt(std::size_t i) const
{
    if (i == 0 || i >= M())
        throw std::out_of_range("StencilAt: interior nodes only (1 <= i <= M-1)");

    const double hm = hMinus(i);
    const double hp = hPlus(i);
    const double hs = hm + hp;

    Stencil s;

    // V_S
    s.dm = -hp / (hm * hs);
    s.d0 = (hp - hm) / (hm * hp);
    s.dp = hm / (hp * hs);

    // V_SS
    s.em = 2.0 / (hm * hs);
    s.e0 = -2.0 / (hm * hp);
    s.ep = 2.0 / (hp * hs);

    return s;
}

std::size_t Grid::NearestIndex(double S) const
{
    // S_ is sorted, so binary search rather than a linear scan
    auto it = std::lower_bound(S_.begin(), S_.end(), S);

    if (it == S_.begin())  return 0;
    if (it == S_.end())    return S_.size() - 1;

    const std::size_t hi = static_cast<std::size_t>(it - S_.begin());
    const std::size_t lo = hi - 1;

    return (S - S_[lo] <= S_[hi] - S) ? lo : hi;
}

double Grid::Interpolate(const std::vector<double>& V, double S) const
{
    if (V.size() != S_.size())
        throw std::invalid_argument("Interpolate: V size does not match grid");

    if (S <= S_.front()) return V.front();
    if (S >= S_.back())  return V.back();

    auto it = std::lower_bound(S_.begin(), S_.end(), S);
    const std::size_t hi = static_cast<std::size_t>(it - S_.begin());
    const std::size_t lo = hi - 1;

    const double w = (S - S_[lo]) / (S_[hi] - S_[lo]);
    return (1.0 - w) * V[lo] + w * V[hi];
}

bool Grid::HasNodeAt(double S, double tol) const
{
    const std::size_t i = NearestIndex(S);
    return std::abs(S_[i] - S) <= tol * std::max(1.0, std::abs(S));
}

void Grid::BuildUniformTimeNodes(double T, std::size_t N)
{
    if (T <= 0.0)  throw std::invalid_argument("Grid: T must be positive");
    if (N < 1)     throw std::invalid_argument("Grid: N must be >= 1");

    t_.resize(N + 1);
    const double dt = T / static_cast<double>(N);

    for (std::size_t n = 0; n <= N; ++n)
        t_[n] = static_cast<double>(n) * dt;

    t_[N] = T;   // kill accumulated round-off at the endpoint
}

void Grid::Validate() const
{
    if (S_.size() < 3)
        throw std::invalid_argument("Grid: need at least 3 spatial nodes");
    if (t_.size() < 2)
        throw std::invalid_argument("Grid: need at least 2 time nodes");

    for (std::size_t i = 1; i < S_.size(); ++i)
        if (!(S_[i] > S_[i - 1]))
            throw std::runtime_error("Grid: spatial nodes not strictly increasing");

    for (std::size_t n = 1; n < t_.size(); ++n)
        if (!(t_[n] > t_[n - 1]))
            throw std::runtime_error("Grid: time nodes not strictly increasing");
}

// ===========================================================================
//  UniformGrid
// ===========================================================================

UniformGrid::UniformGrid(const Domain& domain, std::size_t M, std::size_t N)
{
    domain.Validate();
    if (M < 2) throw std::invalid_argument("UniformGrid: M must be >= 2");

    S_.resize(M + 1);
    const double dS = (domain.Smax - domain.Smin) / static_cast<double>(M);

    for (std::size_t i = 0; i <= M; ++i)
        S_[i] = domain.Smin + static_cast<double>(i) * dS;

    S_[M] = domain.Smax;   // exact endpoint

    BuildUniformTimeNodes(domain.T, N);
    Validate();
}

// ===========================================================================
//  SinhGrid
// ===========================================================================

SinhGrid::SinhGrid(const Domain& domain,
    std::size_t M,
    std::size_t N,
    double focus,
    double density)
    : focus_(focus), density_(density)
{
    domain.Validate();
    if (M < 2)
        throw std::invalid_argument("SinhGrid: M must be >= 2");
    if (density <= 0.0)
        throw std::invalid_argument("SinhGrid: density must be positive");
    if (focus < domain.Smin || focus > domain.Smax)
        throw std::invalid_argument("SinhGrid: focus must lie inside the domain");

    const double c1 = std::asinh((domain.Smin - focus) / density);   // <= 0
    const double c2 = std::asinh((domain.Smax - focus) / density);   // >= 0

    if (!(c2 > c1))
        throw std::runtime_error("SinhGrid: degenerate transformation");

    // xi coordinate of the focus point
    const double xiStar = -c1 / (c2 - c1);

    // Split the M intervals either side of xiStar so that focus is a node.
    // Round to the nearest index, but keep at least one interval on each side
    // unless the focus sits on the domain boundary itself.
    std::size_t j = static_cast<std::size_t>(
        std::llround(xiStar * static_cast<double>(M)));

    if (xiStar > 0.0 && j == 0)  j = 1;
    if (xiStar < 1.0 && j == M)  j = M - 1;

    std::vector<double> xi(M + 1);

    if (j == 0)
    {
        for (std::size_t i = 0; i <= M; ++i)
            xi[i] = static_cast<double>(i) / static_cast<double>(M);
    }
    else if (j == M)
    {
        for (std::size_t i = 0; i <= M; ++i)
            xi[i] = static_cast<double>(i) / static_cast<double>(M);
    }
    else
    {
        // piecewise-uniform in xi: [0, xiStar] with j intervals,
        //                          [xiStar, 1] with M-j intervals
        for (std::size_t i = 0; i <= j; ++i)
            xi[i] = xiStar * static_cast<double>(i) / static_cast<double>(j);

        for (std::size_t i = j; i <= M; ++i)
            xi[i] = xiStar + (1.0 - xiStar)
            * static_cast<double>(i - j) / static_cast<double>(M - j);
    }

    xi[0] = 0.0;
    xi[M] = 1.0;

    S_.resize(M + 1);
    for (std::size_t i = 0; i <= M; ++i)
        S_[i] = focus + density * std::sinh(c1 + xi[i] * (c2 - c1));

    // Pin the three points we care about to exact values
    S_[0] = domain.Smin;
    S_[M] = domain.Smax;
    if (j > 0 && j < M)
        S_[j] = focus;

    BuildUniformTimeNodes(domain.T, N);
    Validate();
}