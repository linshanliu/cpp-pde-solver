#ifndef BLACKSCHOLESPDE_HPP
#define BLACKSCHOLESPDE_HPP


#include "ConvectionDiffusionPDE.hpp"

// ---------------------------------------------------------------------------
//  Black-Scholes operator under GBM:  dS = (r - q) S dt + sigma S dW
//
//      V_t + a(S) V_SS + b(S) V_S + c V = 0
//        a(S) = 0.5 sigma^2 S^2     -> Diffusion
//        b(S) = (r - q) S           -> Convection
//        c    = -r                  -> Reaction   (minus sign folded in)
//
//  Note (S, t) in every signature even though nothing here uses t: local vol
//  will, and changing the interface later means touching every call site.
//
//  final: a barrier is a contract feature, not a change of dynamics. It
//  belongs in Domain + BoundaryCondition, not in a subclass of this.
// ---------------------------------------------------------------------------
class BlackScholesPDE final : public ConvectionDiffusionPDE
{
public:
    BlackScholesPDE(double r, double sigma, double q = 0.0);

    double Diffusion(double S, double t) const override;
    double Convection(double S, double t) const override;
    double Reaction(double S, double t) const override;

    double r()     const { return r_; }
    double Sigma() const { return sigma_; }
    double q()     const { return q_; }

private:
    double r_;
    double sigma_;
    double q_;
};

#endif //!BLACKSCHOLESPDE_HPP