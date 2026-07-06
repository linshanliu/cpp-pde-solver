#ifndef THOMAS_SOLVER_HPP
#define THOMAS_SOLVER_HPP

#include <vector>

class ThomasSolver
{
public:

    static std::vector<double> Solve(
        const std::vector<double>& lower,
        const std::vector<double>& diag,
        const std::vector<double>& upper,
        const std::vector<double>& rhs);
};
 
#endif