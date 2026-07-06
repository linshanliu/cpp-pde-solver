#ifndef PSORSOLVER_HPP
#define PSORSOLVER_HPP

#include <vector>
#include <cstddef>

class PSORSolver
{
public:
    static std::vector<double> Solve(
        const std::vector<double>& lower,
        const std::vector<double>& diag,
        const std::vector<double>& upper,
        const std::vector<double>& rhs,
        const std::vector<double>& payoff,
        double omega = 1.2,
        std::size_t maxIter = 10000,
        double tol = 1e-8);
};

#endif