#ifndef STABILITY_TEST_HPP
#define STABILITY_TEST_HPP

class StabilityTest
{
public:
    static void Run();

private:
    static double NormalCDF(double x);

    static double BlackScholesClosedForm(double S, double K, double T, double r, double sigma, bool isCall, double q = 0.0);
};

#endif