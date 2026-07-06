#ifndef BLACKSCHOLES_PDE_HPP
#define BLACKSCHOLES_PDE_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>

enum class ExerciseType
{
    European,
    American
};



class BlackScholesPDE
{
protected:
    double K_;        // strike
    double T_;        // maturity
    double r_;        // risk-free rate
    double sigma_;    // volatility
    double q_;        // dividend yield
    bool isCall_;

    ExerciseType exerciseType_;

public:
    BlackScholesPDE(double K, double T, double r, double sigma, bool isCall, ExerciseType exerciseType = ExerciseType::European, double q = 0.0);

    virtual ~BlackScholesPDE() = default;

    // coefficient of V_SS
    double Diffusion(double S) const;

    // coefficient of V_S
    double Convection(double S) const;


    // coefficient of V
    double Reaction() const;


    virtual double Payoff(double S) const;


    // left boundary: S = 0
    virtual double BoundaryAtZero(double t) const;


    // right boundary: S = Smax
    virtual double BoundaryAtMax(double Smax, double t) const;


    virtual bool HasEarlyExercise() const;


    double K() const;
    double T() const;
    double r() const;
    double sigma() const;
    double q() const;
    bool IsCall() const;
    ExerciseType Exercise() const;
};

#endif