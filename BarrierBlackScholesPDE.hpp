#ifndef BARRIERBLACKSCHOLES_PDE_HPP
#define BARRIERBLACKSCHOLES_PDE_HPP

#include "BlackScholesPDE.hpp"

enum class BarrierType
{
    DownAndOut,
    UpAndOut
};

class BarrierBlackScholesPDE : public BlackScholesPDE
{
private:
    double barrier_;
    BarrierType barrierType_;

public:
    BarrierBlackScholesPDE(
        double K,
        double T,
        double r,
        double sigma,
        bool isCall,
        double barrier,
        BarrierType barrierType,
        double q = 0.0);

    double Barrier() const;
    BarrierType Type() const;

    double Payoff(double S) const override;

    double BoundaryAtZero(double t) const override;

    double BoundaryAtMax(double Smax, double t) const override;
};

#endif