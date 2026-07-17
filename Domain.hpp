#ifndef DOMAIN_HPP
#define DOMAIN_HPP



#include <stdexcept>

// ---------------------------------------------------------------------------
//  Solution region  (S, t) in [Smin, Smax] x [0, T].
//
//  vanilla        : Smin = 0,  Smax ~ 3-5 x K
//  down-and-out   : Smin = B
//  up-and-out     : Smax = B
//
//  Smax is a truncation of the real half-line [0, inf), so it is a modelling
//  choice, not a fact about the contract.
// ---------------------------------------------------------------------------
struct Domain
{
    double Smin{ 0.0 };
    double Smax{ 0.0 };
    double T{ 0.0 };

    void Validate() const
    {
        if (!(Smax > Smin)) throw std::invalid_argument("Domain: need Smax > Smin");
        if (Smin < 0.0)     throw std::invalid_argument("Domain: Smin must be >= 0");
        if (!(T > 0.0))     throw std::invalid_argument("Domain: T must be positive");
    }
};

#endif // !DOMAIN_HPP