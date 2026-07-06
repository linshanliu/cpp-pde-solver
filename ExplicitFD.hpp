#ifndef EXPLICIT_FD_HPP
#define EXPLICIT_FD_HPP

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "Grid1D.hpp"
#include "BlackScholesPDE.hpp"

class ExplicitFD
{
private:
    const BlackScholesPDE& pde_;
    Grid1D sGrid_;
    std::size_t timeSteps_;
    double dt_;

public:
    ExplicitFD(const BlackScholesPDE& pde, const Grid1D& sGrid, std::size_t timeSteps);
    

    double Price(double S0) const;

private:
    double Interpolate(double S0, const std::vector<double>& values) const;
};

#endif