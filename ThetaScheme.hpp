#ifndef THETA_SCHEME_HPP
#define THETA_SCHEME_HPP

#include "FDMBase.hpp"

// ---------------------------------------------------------------------------
//  One class, one knob. With tau = T - t and L the spatial operator:
//
//      (I - theta*dt*L^{new}) V^{new} = (I + (1-theta)*dt*L^{old}) V^{old}
//
//      theta = 0.0  ->  explicit        O(dt),   conditionally stable
//      theta = 0.5  ->  Crank-Nicolson  O(dt^2), A-stable but NOT L-stable
//                                       -> oscillates on the payoff kink
//      theta = 1.0  ->  fully implicit  O(dt),   L-stable, damps everything
//
//  Rannacher = a handful of theta=1 half-steps up front, then theta=0.5.
//  That is a two-line change here and a rewrite in three separate classes,
//  which is the whole reason this is one class.
// ---------------------------------------------------------------------------
class ThetaScheme : public FDMBase
{
private:
    void ApplyBoundaryRow(std::size_t i, double tNew, double tOld, double dt);
    double theta_;

    // scratch, sized M+1, refilled every step
    std::vector<double> loNew_, diNew_, upNew_;
    std::vector<double> loOld_, diOld_, upOld_;
    std::vector<double> a_, b_, c_, d_;
protected:
    void StepBack(std::size_t n) override;
public:
    ThetaScheme(const ConvectionDiffusionPDE& pde,
        const Payoff& payoff,
        const BoundaryConditions& bc,
        const Grid& grid,
        double                        theta);

    double Theta() const { return theta_; }
};

#endif // THETA_SCHEME_HPP