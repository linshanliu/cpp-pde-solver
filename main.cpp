#include <iostream>
#include <iomanip>
#include <chrono>

#include "Grid1D.hpp"
#include "BlackScholesPDE.hpp"
#include "ExplicitFD.hpp"
#include "ImplicitFD.hpp"
#include "CrankNicolsonFD.hpp"
#include "BarrierBlackScholesPDE.hpp"

#include "ConvergenceTest.hpp"
#include "StabilityTest.hpp"

//int main()
//{
//    try
//    {
//        //--------------------------------------------------
//        // Pricing Example with default parameters
//        //--------------------------------------------------
//
//        double S0 = 100.0;
//        double K = 100.0;
//        double T = 1.0;
//        double r = 0.05;
//        double sigma = 0.20;
//
//
//        // create pde object with default parameters
//        BlackScholesPDE pde(K, T, r, sigma, true);
//
//
//        // create grid object with xMin = 0, xMax = 3K and numSteps = 200 (space state)
//        Grid1D grid(0.0, 3.0 * K, 200);
//
//        // create ExplicitFD solver with corresponding pde, grid and time step =50000
//        ExplicitFD solver(pde, grid, 50000);
//
//        double price = solver.Price(S0);
//
//        std::cout << "=================================\n";
//        std::cout << "European Call Pricing Example\n";
//        std::cout << "=================================\n";
//
//        std::cout << "Price = "
//            << price
//            << std::endl;
//
//        //--------------------------------------------------
//        // Tests
//        //--------------------------------------------------
//
//        ConvergenceTest::Run();
//
//        // Convergence observation:
//        //
//        // As the number of spatial grid intervals M doubles,
//        // the pricing error decreases by approximately a factor of 4:
//        //
//        //     0.0631 -> 0.0151 -> 0.00385 -> 0.00095
//        //
//        // This behaviour is consistent with second-order convergence
//        // in space, since the asset price derivatives are approximated
//        // using central finite differences whose truncation error is
//        // O(dS^2).
//        //
//        // The time discretization error is kept sufficiently small by
//        // choosing the number of time steps proportional to M^2,
//        // ensuring that the spatial discretization error dominates.
//
//        StabilityTest::Run();
//    }
//    catch (const std::exception& e)
//    {
//        std::cerr
//            << "Exception: "
//            << e.what()
//            << std::endl;
//    }
//
//    return 0;
//}


//int main()
//{
//    double S0 = 100.0;
//    double K = 100.0;
//    double T = 1.0;
//    double r = 0.05;
//    double sigma = 0.20;
//
//    std::size_t M = 200;
//    std::size_t N = 5000;
//
//    double Smax = 3.0 * K;
//
//    BlackScholesPDE pde(
//        K,
//        T,
//        r,
//        sigma,
//        true
//    );
//
//    Grid1D grid(
//        0.0,
//        Smax,
//        M
//    );
//
//    ExplicitFD explicitSolver(
//        pde,
//        grid,
//        N
//    );
//
//    ImplicitFD implicitSolver(
//        pde,
//        grid,
//        N
//    );
//
//    double explicitPrice =
//        explicitSolver.Price(S0);
//
//    double implicitPrice =
//        implicitSolver.Price(S0);
//
//    std::cout << std::fixed
//        << std::setprecision(6);
//
//    std::cout << "========================\n";
//    std::cout << "European Call Option\n";
//    std::cout << "========================\n";
//
//    std::cout << "S0    = " << S0 << '\n';
//    std::cout << "K     = " << K << '\n';
//    std::cout << "T     = " << T << '\n';
//    std::cout << "r     = " << r << '\n';
//    std::cout << "sigma = " << sigma << '\n';
//
//    std::cout << "\n";
//
//    std::cout << "Explicit FD Price : "
//        << explicitPrice
//        << '\n';
//
//    std::cout << "Implicit FD Price : "
//        << implicitPrice
//        << '\n';
//
//    return 0;
//}


//int main()
//{
//    double S0 = 100.0;
//    double K = 100.0;
//    double T = 1.0;
//    double r = 0.05;
//    double sigma = 0.20;
//
//    bool isCall = true;
//
//    std::size_t M = 200;
//    std::size_t N = 5000;
//
//    double Smax = 3.0 * K;
//
//    BlackScholesPDE pde(K, T, r, sigma, isCall);
//
//    Grid1D grid(0.0, Smax, M);
//
//    ExplicitFD explicitSolver(pde, grid, N);
//    ImplicitFD implicitSolver(pde, grid, N);
//    CrankNicolsonFD cnSolver(pde, grid, N);
//
//    double explicitPrice = explicitSolver.Price(S0);
//    double implicitPrice = implicitSolver.Price(S0);
//    double cnPrice = cnSolver.Price(S0);
//
//    std::cout << std::fixed << std::setprecision(6);
//
//    std::cout << "========================\n";
//    std::cout << "European Call Option\n";
//    std::cout << "========================\n";
//
//    std::cout << "Explicit FD       : " << explicitPrice << '\n';
//    std::cout << "Implicit FD       : " << implicitPrice << '\n';
//    std::cout << "Crank-Nicolson FD : " << cnPrice << '\n';
//
//    return 0;
//}


// While the explicit finite-difference method is computationally cheap per time step,
// its conditional stability often requires a prohibitively large number of time steps. 
// The implicit scheme offers unconditional stability, whereas the Crank¨CNicolson scheme 
// provides the best balance between stability, accuracy, and computational efficiency, 
// making it the preferred method for practical option pricing applications.
//int main()
//{
//    double S0 = 100.0;
//    double K = 100.0;
//    double T = 1.0;
//    double r = 0.05;
//    double sigma = 0.20;
//
//    bool isCall = true;
//
//    std::size_t M = 200;
//    std::size_t N = 5000;
//
//    Grid1D grid(
//        0.0,
//        3.0 * K,
//        M
//    );
//
//    BlackScholesPDE pde(
//        K,
//        T,
//        r,
//        sigma,
//        isCall
//    );
//
//    ExplicitFD explicitSolver(
//        pde,
//        grid,
//        N
//    );
//
//    ImplicitFD implicitSolver(
//        pde,
//        grid,
//        N
//    );
//
//    CrankNicolsonFD cnSolver(
//        pde,
//        grid,
//        N
//    );
//
//    //---------------------------------------
//    // Explicit
//    //---------------------------------------
//
//    auto start = std::chrono::high_resolution_clock::now();
//
//    double explicitPrice =
//        explicitSolver.Price(S0);
//
//    auto end = std::chrono::high_resolution_clock::now();
//
//    auto explicitTime =
//        std::chrono::duration_cast<
//        std::chrono::microseconds>(
//            end - start
//        );
//
//    //---------------------------------------
//    // Implicit
//    //---------------------------------------
//
//    start = std::chrono::high_resolution_clock::now();
//
//    double implicitPrice =
//        implicitSolver.Price(S0);
//
//    end = std::chrono::high_resolution_clock::now();
//
//    auto implicitTime =
//        std::chrono::duration_cast<
//        std::chrono::microseconds>(
//            end - start
//        );
//
//    //---------------------------------------
//    // Crank Nicolson
//    //---------------------------------------
//
//    start = std::chrono::high_resolution_clock::now();
//
//    double cnPrice =
//        cnSolver.Price(S0);
//
//    end = std::chrono::high_resolution_clock::now();
//
//    auto cnTime =
//        std::chrono::duration_cast<
//        std::chrono::microseconds>(
//            end - start
//        );
//
//    //---------------------------------------
//    // Output
//    //---------------------------------------
//
//    std::cout << std::fixed
//        << std::setprecision(6);
//
//    std::cout << "\n";
//    std::cout << "=============================================\n";
//    std::cout << "Method Comparison\n";
//    std::cout << "=============================================\n";
//
//    std::cout
//        << std::setw(20) << "Method"
//        << std::setw(15) << "Price"
//        << std::setw(15) << "Time(us)"
//        << std::endl;
//
//    std::cout
//        << std::setw(20) << "Explicit FD"
//        << std::setw(15) << explicitPrice
//        << std::setw(15) << explicitTime.count()
//        << std::endl;
//
//    std::cout
//        << std::setw(20) << "Implicit FD"
//        << std::setw(15) << implicitPrice
//        << std::setw(15) << implicitTime.count()
//        << std::endl;
//
//    std::cout
//        << std::setw(20) << "Crank Nicolson"
//        << std::setw(15) << cnPrice
//        << std::setw(15) << cnTime.count()
//        << std::endl;
//
//    return 0;
//}


//int main()
//{
//    //---------------------------------------
//    // Market Parameters
//    //---------------------------------------
//
//    double S0 = 100.0;
//    double K = 100.0;
//    double T = 1.0;
//    double r = 0.05;
//    double sigma = 0.20;
//
//    //---------------------------------------
//    // Grid Parameters
//    //---------------------------------------
//
//    std::size_t M = 200;
//    std::size_t N = 50000;
//
//    Grid1D grid(0.0, 3.0 * K, M);
//
//    //---------------------------------------
//    // European Put
//    //---------------------------------------
//
//    BlackScholesPDE europeanPut(K, T, r, sigma, false, ExerciseType::European);
//
//    ExplicitFD europeanSolver(europeanPut, grid, N);
//
//    double europeanPrice =
//        europeanSolver.Price(S0);
//
//    // Implicit FD
//    ImplicitFD europeanImpSolver(europeanPut, grid, N);
//
//    double europeanImpPrice =
//        europeanImpSolver.Price(S0);
//
//    //CrankNicolson FD
//    CrankNicolsonFD europeanCNSolver(europeanPut, grid, N);
//
//    double europeanCNPrice =
//        europeanCNSolver.Price(S0);
//
//
//    //---------------------------------------
//    // American Put
//    //---------------------------------------
//
//    BlackScholesPDE americanPut(K, T, r, sigma, false, ExerciseType::American);
//
//    ExplicitFD americanSolver(americanPut, grid, N);
//
//    double americanPrice =
//        americanSolver.Price(S0);
//
//
//    //Implicit
//    ImplicitFD americanImpSolver(americanPut, grid, N);
//
//    double americanImpPrice =
//        americanImpSolver.Price(S0);
//
//    // Crank Nicolson FD
//    CrankNicolsonFD americanCNSolver(americanPut, grid, N);
//
//    double americanCNPrice =
//        americanCNSolver.Price(S0);
//    //---------------------------------------
//    // Output
//    //---------------------------------------
//
//    std::cout << std::fixed
//        << std::setprecision(6);
//
//    std::cout << "\n";
//    std::cout << "=================================\n";
//    std::cout << "American Option Test\n";
//    std::cout << "=================================\n";
//
//    std::cout
//        << "European Put Price = "
//        << europeanPrice
//        << std::endl;
//
//    std::cout
//        << "American Put Price = "
//        << americanPrice
//        << std::endl;
//
//    std::cout << "\n";
//
//    if (americanPrice >= europeanPrice)
//    {
//        std::cout
//            << "PASS: American premium detected."
//            << std::endl;
//    }
//    else
//    {
//        std::cout
//            << "FAIL: American option cheaper than European."
//            << std::endl;
//    }
//
//    //Implicit
//    std::cout << std::fixed
//        << std::setprecision(6);
//
//    std::cout << "\n";
//    std::cout << "=================================\n";
//    std::cout << "American Option Test\n";
//    std::cout << "=================================\n";
//
//    std::cout
//        << "European Put Price = "
//        << europeanImpPrice
//        << std::endl;
//
//    std::cout
//        << "American Put Price = "
//        << americanImpPrice
//        << std::endl;
//
//    std::cout << "\n";
//
//    if (americanImpPrice >= europeanImpPrice)
//    {
//        std::cout
//            << "PASS: American premium detected."
//            << std::endl;
//    }
//    else
//    {
//        std::cout
//            << "FAIL: American option cheaper than European."
//            << std::endl;
//    }
//
//    //Crank Nicolson
//    std::cout << std::fixed
//        << std::setprecision(6);
//
//    std::cout << "\n";
//    std::cout << "=================================\n";
//    std::cout << "American Option Test\n";
//    std::cout << "=================================\n";
//
//    std::cout
//        << "European Put Price = "
//        << europeanCNPrice
//        << std::endl;
//
//    std::cout
//        << "American Put Price = "
//        << americanCNPrice
//        << std::endl;
//
//    std::cout << "\n";
//
//    if (americanCNPrice >= europeanCNPrice)
//    {
//        std::cout
//            << "PASS: American premium detected."
//            << std::endl;
//    }
//    else
//    {
//        std::cout
//            << "FAIL: American option cheaper than European."
//            << std::endl;
//    }
//
//    return 0;
//}

int main()
{
    //--------------------------------------------------
    // Market parameters
    //--------------------------------------------------

    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double sigma = 0.20;

    //--------------------------------------------------
    // Vanilla Call
    //--------------------------------------------------

    BlackScholesPDE vanillaPDE(
        K,
        T,
        r,
        sigma,
        true
    );

    Grid1D vanillaGrid(
        0.0,
        500.0,
        1000
    );

    CrankNicolsonFD vanillaSolver(
        vanillaPDE,
        vanillaGrid,
        1000
    );

    double vanillaPrice =
        vanillaSolver.Price(S0);

    //--------------------------------------------------
    // Down-and-Out Call
    //--------------------------------------------------

    double barrier = 96.0;

    BarrierBlackScholesPDE barrierPDE(
        K,
        T,
        r,
        sigma,
        true,
        barrier,
        BarrierType::DownAndOut
    );

    Grid1D barrierGrid(
        barrier,
        500.0,
        1000
    );

    CrankNicolsonFD barrierSolver(
        barrierPDE,
        barrierGrid,
        1000
    );

    double barrierPrice =
        barrierSolver.Price(S0);

    //--------------------------------------------------
    // Output
    //--------------------------------------------------

    std::cout << std::fixed
        << std::setprecision(6);

    std::cout << std::endl;
    std::cout << "============================"
        << std::endl;

    std::cout << "Barrier Option Test"
        << std::endl;

    std::cout << "============================"
        << std::endl;

    std::cout << "Vanilla Call Price      = "
        << vanillaPrice
        << std::endl;

    std::cout << "Down-and-Out Call Price = "
        << barrierPrice
        << std::endl;

    std::cout << std::endl;

    if (barrierPrice < vanillaPrice)
    {
        std::cout
            << "PASS: Barrier option is cheaper."
            << std::endl;
    }
    else
    {
        std::cout
            << "FAIL: Unexpected result."
            << std::endl;
    }

    return 0;
}