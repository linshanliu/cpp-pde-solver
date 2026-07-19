#include <cstdio>

#include "AnalyticBlackScholes.hpp"
#include "BlackScholesPDE.hpp"
#include "BoundaryCondition.hpp"
#include "Domain.hpp"
#include "Grid.hpp"
#include "Payoff.hpp"
#include "ThetaScheme.hpp"
#include "ConvergenceTest.hpp"

int main()
{
    const double K = 100.0, T = 1.0, r = 0.05, sigma = 0.2, q = 0.0, S0 = 100.0;

    BlackScholesPDE    pde(r, sigma, q);
    VanillaPayoff      payoff(K, OptionType::Call);
    Domain             dom{ 0.0, 4.0 * K, T };
    LinearityBC        lo, hi;
    BoundaryConditions bc{ lo, hi };
    UniformGrid        grid(dom, 400, 2000);

    const analytic::BSInputs ref{ S0, K, T, r, sigma, q, OptionType::Call };

    std::printf("           price        delta        gamma\n");
    std::printf("analytic   %10.6f  %10.6f  %10.6f\n",
        analytic::Price(ref), analytic::Delta(ref), analytic::Gamma(ref));

    for (double theta : {1.0, 0.5})
    {
        ThetaScheme fd(pde, payoff, bc, grid, theta);
        fd.Solve();
        std::printf("theta=%.1f   %10.6f  %10.6f  %10.6f\n",
            theta, fd.Price(S0), fd.Delta(S0), fd.Gamma(S0));
    }

    ConvergenceTest::Run();
    return 0;
}