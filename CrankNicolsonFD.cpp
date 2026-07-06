#include "CrankNicolsonFD.hpp"

#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "ThomasSolver.hpp"
#include "PSORSolver.hpp"

CrankNicolsonFD::CrankNicolsonFD(
    const BlackScholesPDE& pde,
    const Grid1D& sGrid,
    std::size_t timeSteps)
    : pde_(pde), sGrid_(sGrid), timeSteps_(timeSteps)
{
    if (timeSteps_ == 0)
        throw std::invalid_argument("CrankNicolsonFD: timeSteps must be positive.");

    dt_ = pde_.T() / static_cast<double>(timeSteps_);
}

double CrankNicolsonFD::Price(double S0) const
{
    const std::size_t M = sGrid_.NumSteps();
    const double dS = sGrid_.dx();
    const double Smax = sGrid_.Max();

    std::vector<double> oldValues(M + 1);
    std::vector<double> newValues(M + 1);

    // Terminal condition: V(S,T) = payoff(S)
    for (std::size_t i = 0; i <= M; ++i)
    {
        oldValues[i] = pde_.Payoff(sGrid_.X(i));
    }

    const std::size_t nInterior = M - 1;

    for (std::size_t n = 0; n < timeSteps_; ++n)
    {
        double tNew = pde_.T() - (n + 1) * dt_;
        double tOld = pde_.T() - n * dt_;

        double leftNew = pde_.BoundaryAtZero(tNew);
        double rightNew = pde_.BoundaryAtMax(Smax, tNew);

        double leftOld = pde_.BoundaryAtZero(tOld);
        double rightOld = pde_.BoundaryAtMax(Smax, tOld);

        std::vector<double> lower(nInterior - 1);
        std::vector<double> diag(nInterior);
        std::vector<double> upper(nInterior - 1);
        std::vector<double> rhs(nInterior);

        for (std::size_t j = 0; j < nInterior; ++j)
        {
            std::size_t i = j + 1;
            double S = sGrid_.X(i);

            double A = pde_.Diffusion(S);
            double B = pde_.Convection(S);
            double C = pde_.Reaction();

            double alpha = A / (dS * dS) - B / (2.0 * dS);
            double beta = -2.0 * A / (dS * dS) + C;
            double gamma = A / (dS * dS) + B / (2.0 * dS);

            // Left-hand side: (I - 0.5 dt L) V_new
            double aL = -0.5 * dt_ * alpha;
            double bL = 1.0 - 0.5 * dt_ * beta;
            double cL = -0.5 * dt_ * gamma;

            // Right-hand side: (I + 0.5 dt L) V_old
            double aR = 0.5 * dt_ * alpha;
            double bR = 1.0 + 0.5 * dt_ * beta;
            double cR = 0.5 * dt_ * gamma;

            diag[j] = bL;

            if (j > 0)
                lower[j - 1] = aL;

            if (j < nInterior - 1)
                upper[j] = cL;

            rhs[j] = bR * oldValues[i];

            if (i > 1)
                rhs[j] += aR * oldValues[i - 1];
            else
                rhs[j] += aR * leftOld;

            if (i < M - 1)
                rhs[j] += cR * oldValues[i + 1];
            else
                rhs[j] += cR * rightOld;

            // Boundary contributions from the left-hand side
            if (i == 1)
                rhs[j] -= aL * leftNew;

            if (i == M - 1)
                rhs[j] -= cL * rightNew;
        }

        /*std::vector<double> interior =
            ThomasSolver::Solve(lower, diag, upper, rhs);*/

        std::vector<double> payoff(nInterior);

        for (std::size_t j = 0; j < nInterior; ++j)
        {
            payoff[j] = pde_.Payoff(sGrid_.X(j + 1));
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

        newValues[0] = leftNew;
        newValues[M] = rightNew;

        for (std::size_t j = 0; j < nInterior; ++j)
        {
            newValues[j + 1] = interior[j];
        }

        oldValues = newValues;
    }

    return Interpolate(S0, oldValues);
}

double CrankNicolsonFD::Interpolate(
    double S0,
    const std::vector<double>& values) const
{
    if (S0 <= sGrid_.Min())
        return values.front();

    if (S0 >= sGrid_.Max())
        return values.back();

    double dS = sGrid_.dx();

    std::size_t i =
        static_cast<std::size_t>((S0 - sGrid_.Min()) / dS);

    double SLeft = sGrid_.X(i);
    double SRight = sGrid_.X(i + 1);

    double VLeft = values[i];
    double VRight = values[i + 1];

    double w = (S0 - SLeft) / (SRight - SLeft);

    return (1.0 - w) * VLeft + w * VRight;
}