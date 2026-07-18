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

// ---------------------------------------------------------------------------
//
//  DERIVATION (uniform spacing h, for intuition).
//  Expand the two neighbours about S_i, with all derivatives taken at S_i:
//
//      V[i+1] = V + h V' + (h^2/2) V'' + (h^3/6) V''' + ...      (step +h)
//      V[i-1] = V - h V' + (h^2/2) V'' - (h^3/6) V''' + ...      (step -h)
//
//  Under +h vs -h, ODD-order terms (V', V''') flip sign; EVEN-order terms
//  (V'') do not. So:
//
//    * SUBTRACT -> the V'' terms cancel, leaving V':
//          V[i+1] - V[i-1] = 2h V' + O(h^3)
//          => V' ~ ( -1/2h ) V[i-1] + 0 V[i] + ( 1/2h ) V[i+1]      (2nd order)
//      First-derivative weights are ANTISYMMETRIC, centre weight is 0:
//      a slope reads the left/right difference; the centre point drops out.
//
//    * ADD -> the V' terms cancel, leaving V'':
//          V[i+1] + V[i-1] = 2V[i] + h^2 V'' + O(h^4)
//          => V'' ~ ( 1/h^2 ) V[i-1] + ( -2/h^2 ) V[i] + ( 1/h^2 ) V[i+1]
//      Second-derivative weights are SYMMETRIC, centre is the -2 main term:
//      curvature reads how far the centre sits from the average of its sides.
//
//  NON-UNIFORM GRID
//  With unequal spacings
//      hm = S_i - S_{i-1}   (backward),   hp = S_{i+1} - S_i   (forward),
//  the +hp / -hm expansions no longer cancel cleanly (the two V'' coeffs
//  hp^2/2 and hm^2/2 differ), so the weights are solved from the pair of
//  Taylor equations directly:
//
//      V_S  (S_i) ~ dm V_{i-1} + d0 V_i + dp V_{i+1}
//      V_SS (S_i) ~ em V_{i-1} + e0 V_i + ep V_{i+1}
//
//  Note d0 = (hp - hm)/(hm*hp) is now NON-ZERO: on a non-uniform grid the
//  centre point does contribute to the first derivative. All six reduce to
//  the 1/2h and 1/h^2 forms above when hm == hp.
// ---------------------------------------------------------------------------

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