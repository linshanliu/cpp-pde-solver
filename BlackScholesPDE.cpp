#include "BlackScholesPDE.hpp"

BlackScholesPDE::BlackScholesPDE(double K, double T, double r, double sigma, bool isCall, ExerciseType exerciseType, double q)
    : K_(K), T_(T), r_(r), sigma_(sigma), q_(q), exerciseType_(exerciseType), isCall_(isCall)
{
    if (K_ <= 0.0)
        throw std::invalid_argument("Strike must be positive.");

    if (T_ <= 0.0)
        throw std::invalid_argument("Maturity must be positive.");

    if (sigma_ <= 0.0)
        throw std::invalid_argument("Volatility must be positive.");
}


// coefficient of V_SS
double BlackScholesPDE::Diffusion(double S) const
{
    return 0.5 * sigma_ * sigma_ * S * S;
}


// coefficient of V_S
double BlackScholesPDE::Convection(double S) const
{
    return (r_ - q_) * S;
}


// coefficient of V
double BlackScholesPDE::Reaction() const
{
    return -r_;
}


double BlackScholesPDE::Payoff(double S) const
{
    if (isCall_)
        return std::max(S - K_, 0.0);
    else
        return std::max(K_ - S, 0.0);
}


// left boundary: S = 0
double BlackScholesPDE::BoundaryAtZero(double t) const
{
    if (isCall_)
    {
        return 0.0;
    }
    else
    {
        return K_ * std::exp(-r_ * (T_ - t));
    }
}


// right boundary: S = Smax
double BlackScholesPDE::BoundaryAtMax(double Smax, double t) const
{
    if (isCall_)
    {
        return Smax * std::exp(-q_ * (T_ - t))
            - K_ * std::exp(-r_ * (T_ - t));
    }
    else
    {
        return 0.0;
    }
}


bool BlackScholesPDE::HasEarlyExercise() const
{
    return exerciseType_ == ExerciseType::American;
}



double BlackScholesPDE::K() const { return K_; }
double BlackScholesPDE::T() const { return T_; }
double BlackScholesPDE::r() const { return r_; }
double BlackScholesPDE::sigma() const { return sigma_; }
double BlackScholesPDE::q() const { return q_; }
bool BlackScholesPDE::IsCall() const { return isCall_; }
ExerciseType BlackScholesPDE::Exercise() const { return exerciseType_; }