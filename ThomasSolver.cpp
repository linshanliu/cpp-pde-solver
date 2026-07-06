#include "ThomasSolver.hpp"

#include <vector>
#include <stdexcept>
#include <cmath>

std::vector<double> ThomasSolver::Solve(
    const std::vector<double>& lower,
    const std::vector<double>& diag,
    const std::vector<double>& upper,
    const std::vector<double>& rhs)
{
    const std::size_t n = diag.size();

    if (n == 0)
        throw std::invalid_argument("ThomasSolver: empty system.");

    if (rhs.size() != n)
        throw std::invalid_argument("ThomasSolver: rhs size mismatch.");

    if (lower.size() != n - 1)
        throw std::invalid_argument("ThomasSolver: lower diagonal size mismatch.");

    if (upper.size() != n - 1)
        throw std::invalid_argument("ThomasSolver: upper diagonal size mismatch.");

    std::vector<double> cPrime(n - 1);
    std::vector<double> dPrime(n);
    std::vector<double> x(n);

    if (std::fabs(diag[0]) < 1e-14)
        throw std::runtime_error("ThomasSolver: zero pivot encountered.");

    if (n > 1)
        cPrime[0] = upper[0] / diag[0];

    dPrime[0] = rhs[0] / diag[0];

    for (std::size_t i = 1; i < n; ++i)
    {
        double denominator = diag[i] - lower[i - 1] * cPrime[i - 1];

        if (std::fabs(denominator) < 1e-14)
            throw std::runtime_error("ThomasSolver: zero pivot encountered.");

        if (i < n - 1)
            cPrime[i] = upper[i] / denominator;

        dPrime[i] =
            (rhs[i] - lower[i - 1] * dPrime[i - 1]) / denominator;
    }

    x[n - 1] = dPrime[n - 1];

    for (std::size_t i = n - 1; i-- > 0;)
    {
        x[i] = dPrime[i] - cPrime[i] * x[i + 1];
    }

    return x;
}