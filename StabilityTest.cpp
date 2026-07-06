#include "StabilityTest.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "Grid1D.hpp"
#include "BlackScholesPDE.hpp"
#include "ExplicitFD.hpp"

double StabilityTest::NormalCDF(double x)
{
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double StabilityTest::BlackScholesClosedForm(double S, double K, double T, double r, double sigma, bool isCall, double q)
{
    double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));

    double d2 = d1 - sigma * std::sqrt(T);

    if (isCall)
    {
        return S * std::exp(-q * T) * NormalCDF(d1) - K * std::exp(-r * T) * NormalCDF(d2);
    }

    return K * std::exp(-r * T) * NormalCDF(-d2) - S * std::exp(-q * T) * NormalCDF(-d1);
}

void StabilityTest::Run()
{
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double sigma = 0.2;
    bool isCall = true;

    double Smax = 3.0 * K;
    std::size_t M = 200;

    BlackScholesPDE pde(K, T, r, sigma, isCall);
    Grid1D grid(0.0, Smax, M);

    double exact = BlackScholesClosedForm(S0, K, T, r, sigma, isCall);

    std::vector<std::size_t> timeSteps = {10, 50, 100, 500, 1000, 5000, 20000};

    std::cout << "\n";
    std::cout << "==============================\n";
    std::cout << "Stability Test: Explicit FD\n";
    std::cout << "==============================\n";

    std::cout << "Closed-form price = " << exact << "\n";
    std::cout << "Space steps M     = " << M << "\n";
    std::cout << "dS                = " << grid.dx() << "\n\n";

    std::cout
        << std::setw(10) << "N"
        << std::setw(15) << "dt"
        << std::setw(15) << "Price"
        << std::setw(15) << "Error"
        << std::setw(15) << "Status"
        << std::endl;

    for (std::size_t N : timeSteps)
    {
        double dt = T / static_cast<double>(N);

        try
        {
            ExplicitFD solver(pde, grid, N);

            double price = solver.Price(S0);
            double error = std::fabs(price - exact);

            bool stable =
                std::isfinite(price)
                &&
                std::fabs(price) < 1e6;

            std::string status = stable ? "OK" : "Unstable";

            std::cout
                << std::setw(10) << N
                << std::setw(15) << dt
                << std::setw(15) << price
                << std::setw(15) << error
                << std::setw(15) << status
                << std::endl;
        }
        catch (const std::exception&)
        {
            std::cout
                << std::setw(10) << N
                << std::setw(15) << dt
                << std::setw(15) << "NaN"
                << std::setw(15) << "NaN"
                << std::setw(15) << "Exception"
                << std::endl;
        }
    }
}