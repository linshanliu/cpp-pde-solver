#include "FDMBase.hpp"

#include <stdexcept>

FDMBase::FDMBase(const ConvectionDiffusionPDE& pde,
    const Payoff& payoff,
    const BoundaryConditions& bc,
    const Grid& grid)
    : pde_(pde), payoff_(payoff), bc_(bc), grid_(grid), V_(grid.S().size(), 0.0)
{
}

void FDMBase::ApplyTerminal()
{
    const auto& S = grid_.S();
    for (std::size_t i = 0; i < S.size(); ++i)
        V_[i] = payoff_(S[i]);
}

const std::vector<double>& FDMBase::Solve()
{
    ApplyTerminal();

    // n = N-1 down to 0; each call moves from t_[n+1] to t_[n]
    for (std::size_t n = grid_.N(); n-- > 0; )
        StepBack(n);

    return V_;
}

void FDMBase::AssembleOperator(double t,
    std::vector<double>& lo,
    std::vector<double>& di,
    std::vector<double>& up) const
{
    const auto& S = grid_.S();
    const std::size_t M = grid_.M();

    for (std::size_t i = 1; i < M; ++i)
    {
        const Stencil st = grid_.StencilAt(i);

        const double a = pde_.Diffusion(S[i], t);
        const double b = pde_.Convection(S[i], t);
        const double c = pde_.Reaction(S[i], t);

        lo[i] = a * st.em + b * st.dm;
        di[i] = a * st.e0 + b * st.d0 + c;
        up[i] = a * st.ep + b * st.dp;
    }
}

void FDMBase::LinearityRow(std::size_t i, double t,
    double& cNeighbour, double& cSelf) const
{
    const auto& S = grid_.S();
    const std::size_t M = grid_.M();

    const double b = pde_.Convection(S[i], t);
    const double c = pde_.Reaction(S[i], t);

    if (i == 0)
    {
        // forward difference: V_S ~ (V_1 - V_0) / h
        const double h = S[1] - S[0];
        cNeighbour = b / h;
        cSelf = c - b / h;
    }
    else if (i == M)
    {
        // backward difference: V_S ~ (V_M - V_{M-1}) / h
        const double h = S[M] - S[M - 1];
        cNeighbour = -b / h;
        cSelf = c + b / h;
    }
    else
    {
        throw std::out_of_range("LinearityRow: boundary nodes only");
    }
}

std::size_t FDMBase::NearestInteriorIndex(double S0) const
{
    std::size_t i = grid_.NearestIndex(S0);
    if (i == 0)          i = 1;
    if (i >= grid_.M())  i = grid_.M() - 1;
    return i;
}

double FDMBase::Price(double S0) const
{
    return grid_.Interpolate(V_, S0);
}

double FDMBase::Delta(double S0) const
{
    const std::size_t i = NearestInteriorIndex(S0);
    const Stencil     st = grid_.StencilAt(i);
    return st.dm * V_[i - 1] + st.d0 * V_[i] + st.dp * V_[i + 1];
}

double FDMBase::Gamma(double S0) const
{
    const std::size_t i = NearestInteriorIndex(S0);
    const Stencil     st = grid_.StencilAt(i);
    return st.em * V_[i - 1] + st.e0 * V_[i] + st.ep * V_[i + 1];
}