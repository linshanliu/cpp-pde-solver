#ifndef IMPLICIT_FD_HPP
#define IMPLICIT_FD_HPP

#include <vector>

#include "Grid1D.hpp"
#include "BlackScholesPDE.hpp"


class ImplicitFD
{
private:

    const BlackScholesPDE& pde_;

    Grid1D grid_;

    std::size_t timeSteps_;

    double dt_;

public:

    ImplicitFD(
        const BlackScholesPDE& pde,
        const Grid1D& grid,
        std::size_t timeSteps);

    double Price(double S0) const;

private:
    double Interpolate(double S0, const std::vector<double>& values) const;

};

#endif