#ifndef CONVECTION_DIFFUSION_PDE_HPP
#define CONVECTION_DIFFUSION_PDE_HPP


// ---------------------------------------------------------------------------
//  Generic 1-D parabolic operator, all terms carrying their own sign:
//
//      V_t + a(S,t) V_SS + b(S,t) V_S + c(S,t) V = 0
//
//  Black-Scholes is one instance; local vol and CEV are others. Nothing here
//  knows about strikes, payoffs, exercise style or the solution domain.
// ---------------------------------------------------------------------------
class ConvectionDiffusionPDE
{
public:
    virtual ~ConvectionDiffusionPDE() = default;

    virtual double Diffusion(double S, double t) const = 0;   // a
    virtual double Convection(double S, double t) const = 0;   // b
    virtual double Reaction(double S, double t) const = 0;   // c
};

#endif // CONVECTION_DIFFUSION_PDE_HPP