#include "ConvergenceTest.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "AnalyticBlackScholes.hpp"
#include "BlackScholesPDE.hpp"
#include "BoundaryCondition.hpp"
#include "Domain.hpp"
#include "Grid.hpp"
#include "Payoff.hpp"
#include "ThetaScheme.hpp"

// ---------------------------------------------------------------------------
//  A note on reading these tables.
//
//  The total error has two independent parts: a SPATIAL part O(h^2) and a
//  TEMPORAL part that is O(dt^2) for Crank-Nicolson and O(dt) for fully
//  implicit. Whichever part is larger is the one a refinement study "sees".
//
//  So a single refinement direction can MISLEAD. If we refine space and time
//  together while the spatial error dominates, even the first-order implicit
//  scheme looks second order -- because we are watching the spatial O(h^2)
//  term, not its temporal O(dt) term. To observe the temporal order we must
//  hold the grid fixed (spatial error frozen and small) and refine time alone.
//
//  This test therefore does both:
//    (A) spatial refinement  -- double M, scale N with M   -> spatial order
//    (B) temporal refinement -- fix a large M, double N     -> temporal order
//
//  Knowing that the observed order depends on which error dominates is itself
//  the point: it is exactly the kind of thing a model validation has to get
//  right before certifying a convergence claim.
// ---------------------------------------------------------------------------

namespace
{
    double PriceCall(double S0, double K, double T, double r, double sigma, double q,
        std::size_t M, std::size_t N, double theta)
    {
        BlackScholesPDE    pde(r, sigma, q);
        VanillaPayoff      payoff(K, OptionType::Call);
        Domain             dom{ 0.0, 4.0 * K, T };
        LinearityBC        lo, hi;
        BoundaryConditions bc{ lo, hi };
        UniformGrid        grid(dom, M, N);

        ThetaScheme fd(pde, payoff, bc, grid, theta);
        fd.Solve();
        return fd.Price(S0);
    }

    void PrintHeader(const char* col1)
    {
        std::cout << std::left
            << std::setw(8) << col1
            << std::setw(18) << "Price"
            << std::setw(16) << "Error"
            << std::setw(10) << "Ratio"
            << "\n" << std::string(52, '-') << "\n";
    }

    void PrintRow(std::size_t param, double price, double error, double prevError)
    {
        std::cout << std::left
            << std::setw(8) << param
            << std::setw(18) << std::fixed << std::setprecision(8) << price
            << std::setw(16) << std::scientific << std::setprecision(3) << error;
        if (prevError > 0.0)
            std::cout << std::fixed << std::setprecision(2) << std::setw(10) << (prevError / error);
        else
            std::cout << std::setw(10) << "-";
        std::cout << "\n";
    }

    // (A) Refine space; scale N with M so the temporal error stays subordinate.
    void SpatialStudy(const char* name, double theta,
        double S0, double K, double T, double r, double sigma, double q,
        double exact)
    {
        std::cout << "\n[A] Spatial refinement -- " << name
            << " (theta = " << theta << "), N = 20*M\n";
        PrintHeader("M");
        double prev = 0.0;
        for (std::size_t M : { 50, 100, 200, 400 })
        {
            const std::size_t N = 20 * M;
            const double price = PriceCall(S0, K, T, r, sigma, q, M, N, theta);
            const double error = std::fabs(price - exact);
            PrintRow(M, price, error, prev);
            prev = error;
        }
    }

    // (B) Fix a large grid; refine time alone to expose the temporal order.
    //     CN -> ratio ~ 4 (second order); implicit -> ratio ~ 2 (first order),
    //     until the frozen spatial error is reached and the ratio collapses.
    void TemporalStudy(const char* name, double theta,
        double S0, double K, double T, double r, double sigma, double q,
        double exact)
    {
        std::cout << "\n[B] Temporal refinement -- " << name
            << " (theta = " << theta << "), M = 800 fixed\n";
        PrintHeader("N");
        double prev = 0.0;
        for (std::size_t N : { 25, 50, 100, 200 })
        {
            const double price = PriceCall(S0, K, T, r, sigma, q, 800, N, theta);
            const double error = std::fabs(price - exact);
            PrintRow(N, price, error, prev);
            prev = error;
        }
    }
}

void ConvergenceTest::Run()
{
    const double S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, sigma = 0.2, q = 0.0;

    const analytic::BSInputs ref{ S0, K, T, r, sigma, q, OptionType::Call };
    const double exact = analytic::Price(ref);

    std::cout << "\n============================================================\n";
    std::cout << "Convergence Test -- European call vs closed-form Black-Scholes\n";
    std::cout << "============================================================\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Closed-form price: " << exact << "\n";

    // Spatial order: both schemes are O(h^2) in space, so both show ratio ~ 4.
    SpatialStudy("Crank-Nicolson", 0.5, S0, K, T, r, sigma, q, exact);
    SpatialStudy("Fully implicit", 1.0, S0, K, T, r, sigma, q, exact);

    // Temporal order: CN is O(dt^2) -> ratio ~ 4; implicit is O(dt) -> ratio ~ 2.
    TemporalStudy("Crank-Nicolson", 0.5, S0, K, T, r, sigma, q, exact);
    TemporalStudy("Fully implicit", 1.0, S0, K, T, r, sigma, q, exact);

    std::cout << "\n";
}