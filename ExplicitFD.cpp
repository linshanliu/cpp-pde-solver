#include "ExplicitFD.hpp"

ExplicitFD::ExplicitFD(const BlackScholesPDE& pde, const Grid1D& sGrid, std::size_t timeSteps)
    : pde_(pde), sGrid_(sGrid), timeSteps_(timeSteps)
{
    if (timeSteps_ == 0)
        throw std::invalid_argument("ExplicitFD: timeSteps must be positive.");

    dt_ = pde_.T() / static_cast<double>(timeSteps_);
}

double ExplicitFD::Price(double S0) const
{
    const std::size_t M = sGrid_.NumSteps();
    const double dS = sGrid_.dx();
    const double Smax = sGrid_.Max();


    // create two vectors to represent the values at 2 different time value
    std::vector<double> oldValues(M + 1);
    std::vector<double> newValues(M + 1);


    // Terminal condition (boundary condition t = T): V(S,T) = payoff(S)
    for (std::size_t i = 0; i <= M; ++i)
    {
        double S = sGrid_.X(i);
        oldValues[i] = pde_.Payoff(S);
    }


    // Step backward in time: t = T -> 0
    for (std::size_t n = 0; n < timeSteps_; ++n)
    {
        double t = pde_.T() - (n + 1) * dt_;

        // Boundary conditions (S = 0 and S= S_max)
        newValues[0] = pde_.BoundaryAtZero(t);
        newValues[M] = pde_.BoundaryAtMax(Smax, t);


        // Interior points
        for (std::size_t i = 1; i < M; ++i)
        {
            // Current asset price corresponding to grid point i
            double S = sGrid_.X(i);


            // PDE coefficients:
            // A(S) : diffusion coefficient (second derivative term)
            // B(S) : convection/drift coefficient (first derivative term)
            // C    : reaction/discount coefficient (zero-order term)
            double A = pde_.Diffusion(S);
            double B = pde_.Convection(S);
            double C = pde_.Reaction();


            // Central difference approximation of
            // ∂V/∂S (Delta)
            double delta = (oldValues[i + 1] - oldValues[i - 1]) / (2.0 * dS);


            // Central difference approximation of
            // ∂²V/∂S² (Gamma)
            double gamma = (oldValues[i + 1] - 2.0 * oldValues[i] + oldValues[i - 1]) / (dS * dS);


            // Explicit Euler time stepping:
            //
            // V^{n+1}_i = V^n_i + Δt [ A V_SS + B V_S + C V ]
            //
            // where V_SS and V_S are approximated by finite differences above.
            newValues[i] = oldValues[i] + dt_ * (A * gamma + B * delta + C * oldValues[i]);
        }


        // Apply early exercise constraint for American options:
        //
        // V(S,t) >= payoff(S)
        //
        // After each time step, project the continuation value
        // back onto the admissible region.
        if (pde_.HasEarlyExercise())
        {
            for (std::size_t i = 0; i <= M; ++i)
            {
                double S = sGrid_.X(i);
                newValues[i] = std::max(newValues[i], pde_.Payoff(S));
            }
        }


        oldValues = newValues;
    }

    return Interpolate(S0, oldValues);
}


double ExplicitFD::Interpolate(double S0, const std::vector<double>&values) const
{
    if (S0 <= sGrid_.Min())
        return values.front();

    if (S0 >= sGrid_.Max())
        return values.back();

    // find S0 falls in which interval
    double dS = sGrid_.dx();
    std::size_t i = static_cast<std::size_t>((S0 - sGrid_.Min()) / dS);


    // left right S interval
    double S_left = sGrid_.X(i);
    double S_right = sGrid_.X(i + 1);


    // left right values for the solution at this interval
    double V_left = values[i];
    double V_right = values[i + 1];


    // weight
    double weight = (S0 - S_left) / (S_right - S_left);

    return (1.0 - weight) * V_left + weight * V_right;
}