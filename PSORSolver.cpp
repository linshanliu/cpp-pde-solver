#include "PSORSolver.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

std::vector<double> PSORSolver::Solve(
    const std::vector<double>& lower,
    const std::vector<double>& diag,
    const std::vector<double>& upper,
    const std::vector<double>& rhs,
    const std::vector<double>& payoff,
    double omega,
    std::size_t maxIter,
    double tol)
{
    const std::size_t N = diag.size();

    if (rhs.size() != N || payoff.size() != N)
        throw std::invalid_argument("PSORSolver: rhs/payoff size mismatch.");

    if (lower.size() != N - 1 || upper.size() != N - 1)
        throw std::invalid_argument("PSORSolver: lower/upper size mismatch.");

    if (omega <= 0.0 || omega >= 2.0)
        throw std::invalid_argument("PSORSolver: omega must be in (0, 2).");

    std::vector<double> x = payoff;

    for (std::size_t iter = 0; iter < maxIter; ++iter)
    {
        double error = 0.0;

        for (std::size_t i = 0; i < N; ++i)
        {
            double left = 0.0;
            double right = 0.0;

            if (i > 0)
                left = lower[i - 1] * x[i - 1];

            if (i + 1 < N)
                right = upper[i] * x[i + 1];

            double y = (rhs[i] - left - right) / diag[i];

            double old = x[i];

            // SOR update
            x[i] = old + omega * (y - old);

            // Projection step for American option
            x[i] = std::max(x[i], payoff[i]);

            error = std::max(error, std::abs(x[i] - old));
        }

        if (error < tol)
            break;
    }

    return x;
}