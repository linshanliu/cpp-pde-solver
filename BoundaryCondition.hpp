#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include <cmath>

//  Boumdary condition type
enum class BCType
{
    Dirichlet,   // V = value
    Neumann,     // dV/dS = value
    Linearity    // d2V/dS2 = 0
};

// Base class for boundary condition
class BoundaryCondition
{
public:
    virtual ~BoundaryCondition() = default;
    virtual BCType Type() const = 0;
    virtual double Value(double S, double t) const = 0;
};




// Zero boundary condition (for call option at S=0 and knock-out barrier)
class ZeroBC final : public BoundaryCondition
{
public:
    BCType Type() const override { return BCType::Dirichlet; }
    double Value(double, double) const override { return 0.0; }
};




// Discounted constant boundary condition (for put option at S=0 the value K*exp(-r*(T-t)))
class DiscountedConstantBC final : public BoundaryCondition
{
private:
    double amount_;
    double r_;
    double T_;

public:
    DiscountedConstantBC(double amount, double r, double T): amount_(amount), r_(r), T_(T) {}
    BCType Type() const override { return BCType::Dirichlet; }
    double Value(double, double t) const override
    {
        return amount_ * std::exp(-r_ * (T_ - t));
    }
};



// Linearity boundary condition (for large S, the second derivative of option value with respect to S is zero)
class LinearityBC final : public BoundaryCondition
{
public:
    BCType Type() const override { return BCType::Linearity; }
    double Value(double, double) const override { return 0.0; }
};


struct BoundaryConditions
{
    const BoundaryCondition& lower;
    const BoundaryCondition& upper;
};


#endif // BOUNDARYCONDITION_HPP