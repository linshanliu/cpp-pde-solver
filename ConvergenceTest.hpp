#ifndef CONVERGENCE_TEST_HPP
#define CONVERGENCE_TEST_HPP

// ---------------------------------------------------------------------------
//  Convergence study for the finite-difference pricer.
//
//  Prices a European call on a sequence of refined grids and compares each
//  result against the closed-form Black-Scholes value (AnalyticBlackScholes,
//  kept strictly outside the solver). Reports the error and, crucially, the
//  RATIO of successive errors -- that ratio is the actual evidence of the
//  convergence order:
//
//      Crank-Nicolson (theta = 0.5): second order  -> error ratio ~ 4
//      Fully implicit (theta = 1.0): first  order  -> error ratio ~ 2
//
//  Grid refinement doubles the number of spatial intervals M each time. The
//  number of time steps N is scaled with M so that the spatial error, not the
//  temporal one, dominates -- see the note in Run().
// ---------------------------------------------------------------------------
class ConvergenceTest
{
public:
    static void Run();
};

#endif // CONVERGENCE_TEST_HPP