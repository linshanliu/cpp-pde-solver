// Payoff.hpp
#ifndef PAYOFF_HPP
#define PAYOFF_HPP

#include <algorithm>

enum class OptionType { Call, Put };


// Base class for option payoff
class Payoff
{
public:
    virtual ~Payoff() = default;
    virtual double operator()(double S) const = 0;
};




// Derived class for vanilla option payoff
class VanillaPayoff final : public Payoff
{

private:
    double K_;
    OptionType type_;


public:
	// Constructor
    VanillaPayoff(double K, OptionType type) : K_(K), type_(type) {}

	// Override the operator() to calculate the payoff
    double operator()(double S) const override
    {
        return type_ == OptionType::Call ? std::max(S - K_, 0.0) : std::max(K_ - S, 0.0);
    }

	// Getters
    double Strike() const { return K_; }
    OptionType Type() const { return type_; }

};




// Derived class for digital option payoff
class DigitalPayoff final : public Payoff
{
private:
    double K_;
    OptionType type_;

public:
	// Constructor
    DigitalPayoff(double K, OptionType type) : K_(K), type_(type) {}

	// Override the operator() to calculate the payoff
    double operator()(double S) const override
    {
        bool in = (type_ == OptionType::Call) ? (S > K_) : (S < K_);
        return in ? 1.0 : 0.0;
    }
};

#endif