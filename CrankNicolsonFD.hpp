#ifndef CRANK_NICOLSON_FD_HPP
#define CRANK_NICOLSON_FD_HPP

#include <vector>
#include <cstddef>

#include "Grid1D.hpp"
#include "BlackScholesPDE.hpp"

class CrankNicolsonFD
{
private:
    const BlackScholesPDE& pde_;
    Grid1D sGrid_;
    std::size_t timeSteps_;
    double dt_;

public:
    CrankNicolsonFD(
        const BlackScholesPDE& pde,
        const Grid1D& sGrid,
        std::size_t timeSteps);

    double Price(double S0) const;

private:
    double Interpolate(
        double S0,
        const std::vector<double>& values) const;
};

#endif