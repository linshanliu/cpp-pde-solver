#ifndef ANALYTIC_BLACKSCHOLES_HPP
#define ANALYTIC_BLACKSCHOLES_HPP

#include <algorithm>
#include <cmath>

#include "Payoff.hpp"   // OptionType

// ---------------------------------------------------------------------------
//  Closed-form European vanilla under GBM. Exists ONLY as a benchmark: every
//  FD result is judged against this. Never call it from the solver.
//
//      d1 = [ ln(S/K) + (r - q + sigma^2/2) T ] / (sigma sqrt(T))
//      d2 = d1 - sigma sqrt(T)
//      C  = S e^{-qT} N(d1) - K e^{-rT} N(d2)
//      P  = K e^{-rT} N(-d2) - S e^{-qT} N(-d1)
// ---------------------------------------------------------------------------
namespace analytic {

    inline double NormCdf(double x)
    {
        static constexpr double kInvSqrt2 = 0.7071067811865475244;
        return 0.5 * std::erfc(-x * kInvSqrt2);
    }

    inline double NormPdf(double x)
    {
        static constexpr double kInvSqrt2Pi = 0.3989422804014326779;
        return kInvSqrt2Pi * std::exp(-0.5 * x * x);
    }

    struct BSInputs
    {
        double S, K, T, r, sigma, q;
        OptionType type;
    };

    inline double D1(const BSInputs& p)
    {
        const double v = p.sigma * std::sqrt(p.T);
        return (std::log(p.S / p.K) + (p.r - p.q + 0.5 * p.sigma * p.sigma) * p.T) / v;
    }

    inline double Price(const BSInputs& p)
    {
        if (p.T <= 0.0 || p.sigma <= 0.0)
            return p.type == OptionType::Call ? std::max(p.S - p.K, 0.0)
            : std::max(p.K - p.S, 0.0);

        const double d1 = D1(p);
        const double d2 = d1 - p.sigma * std::sqrt(p.T);
        const double dfq = std::exp(-p.q * p.T);
        const double dfr = std::exp(-p.r * p.T);

        if (p.type == OptionType::Call)
            return p.S * dfq * NormCdf(d1) - p.K * dfr * NormCdf(d2);
        else
            return p.K * dfr * NormCdf(-d2) - p.S * dfq * NormCdf(-d1);
    }

    inline double Delta(const BSInputs& p)
    {
        const double dfq = std::exp(-p.q * p.T);
        const double d1 = D1(p);
        return p.type == OptionType::Call ? dfq * NormCdf(d1)
            : -dfq * NormCdf(-d1);
    }

    inline double Gamma(const BSInputs& p)
    {
        const double dfq = std::exp(-p.q * p.T);
        return dfq * NormPdf(D1(p)) / (p.S * p.sigma * std::sqrt(p.T));
    }

}  // namespace analytic

#endif // ANALYTIC_BLACKSCHOLES_HPP