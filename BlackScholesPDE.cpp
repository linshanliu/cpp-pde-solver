#include "BlackScholesPDE.hpp"

#include <stdexcept>

BlackScholesPDE::BlackScholesPDE(double r, double sigma, double q)
    : r_(r), sigma_(sigma), q_(q)
{
    if (!(sigma > 0.0))
        throw std::invalid_argument("BlackScholesPDE: sigma must be positive");
}

double BlackScholesPDE::Diffusion(double S, double) const
{
    return 0.5 * sigma_ * sigma_ * S * S;
}

double BlackScholesPDE::Convection(double S, double) const
{
    return (r_ - q_) * S;
}

double BlackScholesPDE::Reaction(double, double) const
{
    return -r_;
}