#include "BarrierBlackScholesPDE.hpp"

#include <stdexcept>

BarrierBlackScholesPDE::BarrierBlackScholesPDE(
    double K,
    double T,
    double r,
    double sigma,
    bool isCall,
    double barrier,
    BarrierType barrierType,
    double q)
    : BlackScholesPDE(
        K,
        T,
        r,
        sigma,
        isCall,
        ExerciseType::European,
        q),
    barrier_(barrier),
    barrierType_(barrierType)
{
    if (barrier_ <= 0.0)
        throw std::invalid_argument(
            "BarrierBlackScholesPDE: barrier must be positive.");
}

double BarrierBlackScholesPDE::Barrier() const
{
    return barrier_;
}

BarrierType BarrierBlackScholesPDE::Type() const
{
    return barrierType_;
}

double BarrierBlackScholesPDE::Payoff(double S) const
{
    if (barrierType_ == BarrierType::DownAndOut && S <= barrier_)
        return 0.0;

    if (barrierType_ == BarrierType::UpAndOut && S >= barrier_)
        return 0.0;

    return BlackScholesPDE::Payoff(S);
}

double BarrierBlackScholesPDE::BoundaryAtZero(double t) const
{
    if (barrierType_ == BarrierType::DownAndOut)
    {
        // For Down-and-Out, the left boundary of the grid
        // should be the barrier B, and V(B,t)=0.
        return 0.0;
    }

    return BlackScholesPDE::BoundaryAtZero(t);
}

double BarrierBlackScholesPDE::BoundaryAtMax(double Smax, double t) const
{
    if (barrierType_ == BarrierType::UpAndOut)
    {
        // For Up-and-Out, the right boundary of the grid
        // should be the barrier B, and V(B,t)=0.
        return 0.0;
    }

    return BlackScholesPDE::BoundaryAtMax(Smax, t);
}