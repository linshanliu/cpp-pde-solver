#ifndef FDM_BASE_HPP
#define FDM_BASE_HPP

#include <cstddef>
#include <vector>

#include "BoundaryCondition.hpp"
#include "ConvectionDiffusionPDE.hpp"
#include "Grid.hpp"
#include "Payoff.hpp"

// ---------------------------------------------------------------------------
//  Glue. Owns nothing but the solution vector; every ingredient is injected.
//
//  Marches V backwards from t = T (payoff) to t = 0, one time step at a time.
//  The concrete stepping rule lives in StepBack().
// ---------------------------------------------------------------------------
class FDMBase
{
protected:
    // One step: t_[n+1] (known) -> t_[n] (unknown).
    virtual void StepBack(std::size_t n) = 0;

    void ApplyTerminal();

    // Fills interior rows 1..M-1 of the spatial operator L at time t:
    //     (L V)_i = lo[i] V_{i-1} + di[i] V_i + up[i] V_{i+1}
    void AssembleOperator(double t,
        std::vector<double>& lo,
        std::vector<double>& di,
        std::vector<double>& up) const;

    // Reduced operator row at boundary node i (0 or M) with V_SS dropped:
    //     (L V)_i = cNeighbour * V_j + cSelf * V_i,   j = 1 or M-1
    void LinearityRow(std::size_t i, double t,
        double& cNeighbour, double& cSelf) const;

    std::size_t NearestInteriorIndex(double S0) const;

    const ConvectionDiffusionPDE& pde_;
    const Payoff& payoff_;
    const BoundaryConditions& bc_;
    const Grid& grid_;

    std::vector<double> V_;

public:
    FDMBase(const ConvectionDiffusionPDE& pde,
        const Payoff& payoff,
        const BoundaryConditions& bc,
        const Grid& grid);

    virtual ~FDMBase() = default;

    // Runs the full backward sweep. Returns V(S, 0) on the grid.
    const std::vector<double>& Solve();

    const std::vector<double>& Solution() const { return V_; }

    double Price(double S0) const;

    // Read straight off the grid -- no bump-and-revalue. Both are evaluated at
    // the node nearest S0, so put S0 on a node if you care about the last
    // digit (SinhGrid with focus = S0 does this for free).
    double Delta(double S0) const;
    double Gamma(double S0) const;
};

#endif // FDM_BASE_HPP