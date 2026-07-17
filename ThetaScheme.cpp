#include "ThetaScheme.hpp"
#include "ThomasSolver.hpp"

#include <stdexcept>

ThetaScheme::ThetaScheme(const ConvectionDiffusionPDE& pde,
    const Payoff& payoff,
    const BoundaryConditions& bc,
    const Grid& grid,
    double                        theta)
    : FDMBase(pde, payoff, bc, grid), theta_(theta)
{
    if (theta < 0.0 || theta > 1.0)
        throw std::invalid_argument("ThetaScheme: theta must lie in [0, 1]");

    const std::size_t n = grid.S().size();
    loNew_.assign(n, 0.0); diNew_.assign(n, 0.0); upNew_.assign(n, 0.0);
    loOld_.assign(n, 0.0); diOld_.assign(n, 0.0); upOld_.assign(n, 0.0);
    a_.assign(n, 0.0); b_.assign(n, 0.0); c_.assign(n, 0.0); d_.assign(n, 0.0);
}

void ThetaScheme::StepBack(std::size_t n)
{
    const std::size_t M = grid_.M();
    const double      tOld = grid_.t()[n + 1];   // known layer
    const double      tNew = grid_.t()[n];       // layer we are solving for
    const double      dt = tOld - tNew;

    // Coefficients are frozen at each layer separately. For BS this is
    // redundant (a, b, c do not depend on t) but it is what makes local vol
    // drop straight in later.
    AssembleOperator(tNew, loNew_, diNew_, upNew_);
    AssembleOperator(tOld, loOld_, diOld_, upOld_);

    for (std::size_t i = 1; i < M; ++i)
    {
        a_[i] = -theta_ * dt * loNew_[i];
        b_[i] = 1.0 - theta_ * dt * diNew_[i];
        c_[i] = -theta_ * dt * upNew_[i];

        d_[i] = V_[i] + (1.0 - theta_) * dt
            * (loOld_[i] * V_[i - 1] + diOld_[i] * V_[i] + upOld_[i] * V_[i + 1]);
    }

    ApplyBoundaryRow(0, tNew, tOld, dt);
    ApplyBoundaryRow(M, tNew, tOld, dt);

    a_[0] = 0.0;   // unused by Thomas, zeroed so a debugger shows nothing odd
    c_[M] = 0.0;

    thomas::SolveInPlace(a_, b_, c_, d_);
    V_ = d_;
}

void ThetaScheme::ApplyBoundaryRow(std::size_t i, double tNew, double tOld, double dt)
{
    const auto& S = grid_.S();
    const std::size_t M = grid_.M();
    const bool        lower = (i == 0);

    const BoundaryCondition& cond = lower ? bc_.lower : bc_.upper;

    switch (cond.Type())
    {
    case BCType::Dirichlet:
    {
        if (lower) { b_[0] = 1.0; c_[0] = 0.0; }
        else { a_[M] = 0.0; b_[M] = 1.0; }
        d_[i] = cond.Value(S[i], tNew);
        break;
    }

    case BCType::Neumann:
    {
        // first-order one-sided; the boundary is far from anything we care
        // about, so the extra order is not worth the fill-in
        const double h = lower ? (S[1] - S[0]) : (S[M] - S[M - 1]);
        if (lower) { b_[0] = -1.0 / h; c_[0] = 1.0 / h; }
        else { a_[M] = -1.0 / h; b_[M] = 1.0 / h; }
        d_[i] = cond.Value(S[i], tNew);
        break;
    }

    case BCType::Linearity:
    {
        // V_SS := 0, substituted back into the PDE. Stays tridiagonal.
        // At S = 0 the BS operator has a = b = 0, so this collapses to
        // V_t - r V = 0 all by itself -- the degenerate Fichera boundary.
        double nNew, sNew, nOld, sOld;
        LinearityRow(i, tNew, nNew, sNew);
        LinearityRow(i, tOld, nOld, sOld);

        const std::size_t j = lower ? 1 : M - 1;

        if (lower) { b_[0] = 1.0 - theta_ * dt * sNew; c_[0] = -theta_ * dt * nNew; }
        else { b_[M] = 1.0 - theta_ * dt * sNew; a_[M] = -theta_ * dt * nNew; }

        d_[i] = V_[i] + (1.0 - theta_) * dt * (sOld * V_[i] + nOld * V_[j]);
        break;
    }
    }
}