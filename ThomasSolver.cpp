#include "ThomasSolver.hpp"

#include <cmath>
#include <stdexcept>

namespace thomas {

    void SolveInPlace(const std::vector<double>& a,
        std::vector<double>& b,
        const std::vector<double>& c,
        std::vector<double>& d)
    {
        const std::size_t n = b.size();
        if (n == 0)
            return;
        if (a.size() != n || c.size() != n || d.size() != n)
            throw std::invalid_argument("thomas: size mismatch");

        constexpr double kTiny = 1e-300;

        for (std::size_t i = 1; i < n; ++i)
        {
            if (std::abs(b[i - 1]) < kTiny)
                throw std::runtime_error("thomas: zero pivot in forward sweep");

            const double m = a[i] / b[i - 1];
            b[i] -= m * c[i - 1];
            d[i] -= m * d[i - 1];
        }

        if (std::abs(b[n - 1]) < kTiny)
            throw std::runtime_error("thomas: zero pivot at last row");

        d[n - 1] /= b[n - 1];

        for (std::size_t i = n - 1; i-- > 0; )
            d[i] = (d[i] - c[i] * d[i + 1]) / b[i];
    }

}  // namespace thomas