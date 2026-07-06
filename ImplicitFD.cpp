#include "ImplicitFD.hpp"

#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "ThomasSolver.hpp"
#include "PSORSolver.hpp"

ImplicitFD::ImplicitFD(
    const BlackScholesPDE& pde,
    const Grid1D& sGrid,
    std::size_t timeSteps)
    : pde_(pde), grid_(sGrid), timeSteps_(timeSteps)
{
    if (timeSteps_ == 0)
        throw std::invalid_argument("ImplicitFD: timeSteps must be positive.");

    dt_ = pde_.T() / static_cast<double>(timeSteps_);
}

double ImplicitFD::Price(double S0) const
{
    const std::size_t M = grid_.NumSteps();
    const double dS = grid_.dx();
    const double Smax = grid_.Max();

    std::vector<double> oldValues(M + 1);
    std::vector<double> newValues(M + 1);

    // Terminal condition: V(S,T) = payoff(S)
    for (std::size_t i = 0; i <= M; ++i)
    {
        oldValues[i] = pde_.Payoff(grid_.X(i));
    }

    // We solve only for interior nodes:
    // i = 1, ..., M - 1.
    const std::size_t nInterior = M - 1;

    for (std::size_t n = 0; n < timeSteps_; ++n)
    {
        double t = pde_.T() - (n + 1) * dt_;

        double Boundary0 = pde_.BoundaryAtZero(t);
        double BoundaryMax = pde_.BoundaryAtMax(Smax, t);

        std::vector<double> lower(nInterior - 1);
        std::vector<double> diag(nInterior);
        std::vector<double> upper(nInterior - 1);
        std::vector<double> rhs(nInterior);

        for (std::size_t j = 0; j < nInterior; ++j)
        {
            std::size_t i = j + 1;
            double S = grid_.X(i);

            double A = pde_.Diffusion(S);
            double B = pde_.Convection(S);
            double C = pde_.Reaction();

            double alpha = A / (dS * dS) - B / (2.0 * dS);
            double beta = -2.0 * A / (dS * dS) + C;
            double gamma = A / (dS * dS) + B / (2.0 * dS);

            // Implicit backward-time scheme:
            // (I - dt L) V_new = V_old
            double a = -dt_ * alpha;
            double b = 1.0 - dt_ * beta;
            double c = -dt_ * gamma;

            diag[j] = b;

            if (j > 0)
                lower[j - 1] = a;

            if (j < nInterior - 1)
                upper[j] = c;

            rhs[j] = oldValues[i];
        }

        // Boundary contributions
        rhs[0] -= (-dt_ *
            (pde_.Diffusion(grid_.X(1)) / (dS * dS)
                - pde_.Convection(grid_.X(1)) / (2.0 * dS)))
            * Boundary0;

        rhs[nInterior - 1] -= (-dt_ *
            (pde_.Diffusion(grid_.X(M - 1)) / (dS * dS)
                + pde_.Convection(grid_.X(M - 1)) / (2.0 * dS)))
            * BoundaryMax;

        /*std::vector<double> interior =
            ThomasSolver::Solve(lower, diag, upper, rhs);*/


        std::vector<double> payoff(nInterior);

        for (std::size_t j = 0; j < nInterior; ++j)
        {
            payoff[j] = pde_.Payoff(grid_.X(j + 1));
        }

        std::vector<double> interior;

        if (pde_.HasEarlyExercise())
        {
            interior = PSORSolver::Solve(
                lower,
                diag,
                upper,
                rhs,
                payoff
            );
        }
        else
        {
            interior = ThomasSolver::Solve(
                lower,
                diag,
                upper,
                rhs
            );
        }


        newValues[0] = Boundary0;
        newValues[M] = BoundaryMax;

        for (std::size_t j = 0; j < nInterior; ++j)
        {
            newValues[j + 1] = interior[j];
        }

        oldValues = newValues;
    }

    return Interpolate(S0, oldValues);
}

double ImplicitFD::Interpolate(
    double S0,
    const std::vector<double>& values) const
{
    if (S0 <= grid_.Min())
        return values.front();

    if (S0 >= grid_.Max())
        return values.back();

    double dS = grid_.dx();

    std::size_t i =
        static_cast<std::size_t>((S0 - grid_.Min()) / dS);

    double SLeft = grid_.X(i);
    double SRight = grid_.X(i + 1);

    double VLeft = values[i];
    double VRight = values[i + 1];

    double w = (S0 - SLeft) / (SRight - SLeft);

    return (1.0 - w) * VLeft + w * VRight;
}