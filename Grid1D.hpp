#ifndef GRID1D_HPP
#define GRID1D_HPP

#include <vector>
#include <stdexcept>

class Grid1D
{
private:
    double xMin_;
    double xMax_;
    std::size_t numSteps_;
    double dx_;                          // Step size
    std::vector<double> points_;

public:
    Grid1D(double xMin, double xMax, std::size_t numSteps): xMin_(xMin), xMax_(xMax), numSteps_(numSteps)
    {
        if (xMax <= xMin)
            throw std::invalid_argument("Grid1D: xMax must be greater than xMin.");

        if (numSteps_ == 0)
            throw std::invalid_argument("Grid1D: numSteps must be positive.");

        dx_ = (xMax_ - xMin_) / static_cast<double>(numSteps_);

        points_.resize(numSteps_ + 1);

        for (std::size_t i = 0; i <= numSteps_; ++i)
        {
            points_[i] = xMin_ + i * dx_;
        }
    }


    double X(std::size_t i) const
    {
        if (i >= points_.size())
            throw std::out_of_range("Grid1D: index out of range.");

        return points_[i];
    }


    double dx() const
    {
        return dx_;
    }


    std::size_t NumSteps() const
    {
        return numSteps_;
    }


    std::size_t Size() const
    {
        return points_.size();
    }


    double Min() const
    {
        return xMin_;
    }


    double Max() const
    {
        return xMax_;
    }


    const std::vector<double>& Points() const
    {
        return points_;
    }
};

#endif