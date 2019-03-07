//
// Created by joao on 03/03/2019.
//

#ifndef LINALG_BASIC_MATRIX_H
#define LINALG_BASIC_MATRIX_H

#include "matrix_base.h"
#include <numeric>

namespace numeric {

    template<typename __data_t>
    class matrix : public matrix_base<__data_t> {
    public:
        using data_t = __data_t;
        using base = typename matrix_base<data_t>::base;
        using matrix_base<data_t>::matrix_base;
    };

    template<typename __data_t>
    class square_matrix : public matrix<__data_t> {
    public:
        using data_t = __data_t;
        using base = typename matrix<data_t>::base;
    private:
        using matrix<data_t>::matrix;

        static bool is_zero(const data_t &value,
                            std::enable_if_t<not std::is_integral_v<data_t>> eps) {
            if constexpr (std::is_integral_v<data_t>) {
                if (value == data_t(0))
                    return true;
            } else {
                assert(eps > data_t(0));
                if (std::abs(value) < eps)
                    return true;
            }
            return false;
        }

    public:

        square_matrix(square_matrix<data_t>&& rhs) {
            assert(this->_size.n == this->_size.m);
        }

        square_matrix(const square_matrix<data_t>& rhs) {
            assert(this->_size.n == this->_size.m);
        }

        explicit constexpr square_matrix(size_t n)
                : matrix<data_t>({n, n}) {}

        vector_view<data_t> get_diagonal() const {
            auto d = vector_view<data_t>(this->_size.n);
            for (size_t i = 0u; i < this->_size.n; ++i) {
                d[i] = (*this)(i, i);
            }
            return d;
        }

        data_t get_trace() const {
            const auto d = get_diagonal();
            return std::accumulate(d.begin(), d.end(), data_t(0));
        }

        static square_matrix identity(size_t n) {
            auto i = square_matrix(n);
            for (size_t j = 0u; j < n; ++j) {
                i(j, j) = data_t(1);
            }
            return i;
        }

        bool is_upper_triangular(std::enable_if_t<not std::is_integral_v<data_t>> eps = data_t(1e-2)) const {
            for (size_t j = 0u; j < this->_size.m; ++j) {
                for (size_t i = j + 1u; i < this->_size.n; ++i) {
                    if (not is_zero((*this)(i, j)))
                        return false;
                }
            }
            return true;
        }

        bool is_lower_triangular(std::enable_if_t<not std::is_integral_v<data_t>> eps = data_t(1e-2)) const {
            for (size_t j = 0u; j < this->_size.m; ++j) {
                for (size_t i = 0; i < j; ++i) {
                    if (not is_zero((*this)(i, j)))
                        return false;
                }
            }
            return true;
        }

        bool is_diagonal(std::enable_if_t<not std::is_integral_v<data_t>> eps = data_t(1e-2)) const {
            for (size_t j = 0u; j < this->_size.m; ++j) {
                for (size_t i = 0; i < this->_size.n; ++i) {
                    if (i == j)
                        continue;
                    if (not is_zero((*this)(i, j)))
                        return false;
                }
            }
            return true;
        }
    };

}
#endif //LINALG_BASIC_MATRIX_H
