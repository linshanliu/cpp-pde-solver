#ifndef CONVERGENCE_TEST_HPP
#define CONVERGENCE_TEST_HPP

#include <vector>

class ConvergenceTest
{
public:

    static void Run();

private:

    static double BlackScholesClosedForm(
        double S,
        double K,
        double T,
        double r,
        double sigma,
        bool isCall,
        double q = 0.0);

    static double NormalCDF(double x);
};

#endif