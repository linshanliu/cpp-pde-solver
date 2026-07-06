#include "ConvergenceTest.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>

#include "Grid1D.hpp"
#include "BlackScholesPDE.hpp"
#include "ExplicitFD.hpp"


double ConvergenceTest::NormalCDF(double x)
{
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}


double ConvergenceTest::BlackScholesClosedForm(double S, double K, double T, double r, double sigma, bool isCall, double q)
{
    double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));

    double d2 = d1 - sigma * std::sqrt(T);

    if (isCall)
    {
        return S * std::exp(-q * T) * NormalCDF(d1) - K * std::exp(-r * T) * NormalCDF(d2);
    }

    return K * std::exp(-r * T) * NormalCDF(-d2) - S * std::exp(-q * T) * NormalCDF(-d1);
}


void ConvergenceTest::Run()
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double sigma = 0.2;

    //  Black Scholes pde object with default parameters
    BlackScholesPDE pde(K, T, r, sigma, true);


    // exact solution
    double exact = BlackScholesClosedForm(S0, K, T, r, sigma, true);


    // 4 different number of steps for state space
    std::vector<int> spatialSteps = {50, 100, 200, 400};

    std::cout << "\n";
    std::cout << "==============================\n";
    std::cout << "Convergence Test\n";
    std::cout << "==============================\n";


    // set width of cout
    std::cout
        << std::setw(10) << "M(spatialSteps)"
        << std::setw(15) << "Price"
        << std::setw(15) << "Error"
        << std::endl;

    for (auto M : spatialSteps)
    {
        // create grid object with xMin = 0, xMax = 300 and number of steps (state space) = M
        Grid1D grid(0.0, 300.0, M);


        // create Explicit finit difference object with the corresponding pde, grid and time step
        // here with choose the time step depending on state step M (why?)
        // For the convergence test, the number of time steps is chosen
        // proportional to M^2, where M is the number of spatial grid
        // intervals.
        //
        // Since the explicit finite difference scheme requires
        // dt = O(dS^2) for stability and dS = O(1/M), choosing
        //
        //     N = 5 * M^2
        //
        // ensures that dt decreases at the same rate as dS^2.
        //
        // This keeps the temporal discretisation error sufficiently
        // small so that the observed error is dominated by the spatial
        // discretisation, allowing the expected second-order spatial
        // convergence to be clearly observed.
        ExplicitFD solver(pde, grid, 5 * M * M);


        double price = solver.Price(S0);

        double error = std::fabs(price - exact);

        std::cout
            << std::setw(10) << M
            << std::setw(15) << price
            << std::setw(15) << error
            << std::endl;
    }
}