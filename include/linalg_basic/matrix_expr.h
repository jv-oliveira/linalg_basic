//
// Created by joao on 03/03/2019.
//

#ifndef LINALG_BASIC_MATRIX_EXPR_H
#define LINALG_BASIC_MATRIX_EXPR_H

#include "forward_decl.h"
#include "utils.h"

#include <type_traits>
#include <cstddef>
#include <functional>
#include <ostream>
#include <typeinfo>
#include <cxxabi.h>
#include <string.h>

namespace linalg {

    struct dim_t {
        size_t n;
        size_t m;

        bool operator==(const dim_t &rhs) const {
            return this->m == rhs.m &&
                   this->n == rhs.n;
        }

        bool operator!=(const dim_t &rhs) const {
            return not(rhs == *this);
        }

        friend std::ostream &operator<<(std::ostream &os, const dim_t &dim) {
            os << "{" << dim.n << ", " << dim.m << "}";
            return os;
        }
    };

    template<typename data_t>
    class matrix_expr {
        template<typename data_type>
        class iterator_base;

        static std::true_type is_derived_from_matrix_expr_impl(const matrix_expr *impl);

        static std::false_type is_derived_from_matrix_expr_impl(...);

    public:
        template<class Derived>
        using is_derived_from_matrix_expr_t = decltype(is_derived_from_matrix_expr_impl(std::declval<Derived *>()));

        template<class Derived>
        static inline constexpr bool is_derived_from_matrix_expr_v = is_derived_from_matrix_expr_t<Derived>::value;

        using iterator = iterator_base<data_t>;
        using const_iterator = iterator_base<const data_t>;

        virtual data_t operator()(size_t i, size_t j) const = 0;

        virtual dim_t get_size() const = 0;

        virtual size_t get_raw_size() const = 0;

        const_iterator cbegin() const { return const_iterator(*this); }

        const_iterator cend() const { return const_iterator(*this, get_size()); }

        iterator begin() { return iterator(*this); }

        iterator end() { return iterator(*this, get_size()); }

        const_iterator begin() const { return cbegin(); }

        const_iterator end() const { return cend(); }

        template<typename Derived, typename = std::enable_if_t<is_derived_from_matrix_expr_v<Derived> > >
        friend std::ostream &operator<<(std::ostream &os, const Derived &expr) {
            auto mangled_str = typeid(Derived).name();
            std::unique_ptr<char, decltype(&::free)> str(abi::__cxa_demangle(mangled_str, nullptr, nullptr, nullptr),
                                                         &::free);
            const auto &conv_expr = static_cast<const matrix_expr &>(expr);
            os << (str ? str.get() : mangled_str) << " of size " << conv_expr.get_size() << "\n";
            os << "{\n";
            auto[n, m] = conv_expr.get_size();
            for (size_t i = 0u; i < n; ++i) {
                os << "\t{";
                for (size_t j = 0u; j < m; ++j) {
                    os << conv_expr(i, j);
                    if (j + 1 < m)
                        os << ", ";
                }
                os << "}";
                if (i + 1 < n)
                    os << ",";
                os << "\n";
            }
            os << "}";
            return os;
        }

    private:
        template<typename data_type>
        class iterator_base : public std::iterator<std::random_access_iterator_tag, data_type> {
            using base_class = std::iterator<std::random_access_iterator_tag, data_type>;
            using iterator_category = typename base_class::iterator_category;
            using value_type        = typename base_class::value_type;
            using difference_type   = typename base_class::difference_type;
            using pointer           = typename base_class::pointer;
            using reference         = typename base_class::reference;

            matrix_expr &_matrix;
            dim_t _index{0, 0};
            pointer _ptr;

            void increment_index(difference_type i = 1) {
                _index.n = (_index.n + _index.m + i) / _matrix.get_size().m;
                _index.m = (_index.m + i) % _matrix.get_size().m;
                _ptr = &_matrix(_index.n, _index.m);
            }

            void decrement_index(difference_type i = 1) {
                _index.n = (_index.n + _index.m - i) / _matrix.get_size().m;
                _index.m = (_index.m - i) % _matrix.get_size().m;
                _ptr = &_matrix(_index.n, _index.m);
            }

        private:
            iterator_base(matrix_expr &matrix, const dim_t &index = {0, 0}, const data_type *ptr = nullptr)
                    :
                    _matrix(matrix),
                    _index(index),
                    _ptr(ptr) {
                if (_index != _matrix.get_size()) {
                    if (_ptr == nullptr)
                        _ptr = &_matrix(index.n, index.m);
                    assert(_ptr == &matrix(index.n, index.m));
                } else {
                    ptr = nullptr;
                }
            }

            iterator_base(const iterator_base &it) = default;

            iterator_base(iterator_base &&it) = default;

            ~iterator_base() = default;

            iterator_base &operator=(const iterator_base &it) = default;

            iterator_base &operator=(pointer ptr) {
                for (size_t i = 0u; i < _matrix.get_size().n; ++i) {
                    for (size_t j = 0u; j < _matrix.get_size().m; ++j) {
                        if (ptr == &_matrix(i, j)) {
                            _ptr = ptr;
                            _index = {i, j};
                        }
                    }
                }
                return (*this);
            }

            operator bool() const { return _ptr ? true : false; }

            iterator_base<data_type> &operator+=(const difference_type &movement) {
                increment_index(movement);
                return (*this);
            }

            iterator_base<data_type> &operator-=(const difference_type &movement) {
                increment_index(movement);
                return (*this);
            }

            iterator_base<data_type> &operator++() {
                increment_index();
                return (*this);
            }

            iterator_base<data_type> &operator--() {
                decrement_index();
                return (*this);
            }

            iterator_base<data_type> operator++(difference_type) {
                auto temp(*this);
                temp.increment_index();
                return temp;
            }

            iterator_base<data_type> operator--(difference_type) {
                auto temp(*this);
                temp.decrement_index();
                return temp;
            }

            iterator_base<data_type> operator+(const difference_type &movement) {
                auto temp(*this);
                temp.increment_index(movement);
                return temp;
            }

            iterator_base<data_type> operator-(const difference_type &movement) {
                auto temp(*this);
                temp.decrement_index(movement);
                return temp;
            }

            difference_type operator-(const iterator_base<data_type> &it) {
                return it._index.n * it._index.m -
                       this->_index.n * this->_index.m;
            }

            reference operator*() { return *_ptr; }

            const reference operator*() const { return *_ptr; }

            pointer operator->() { return _ptr; }

            pointer get_ptr() const { return _ptr; }

            const pointer get_const_ptr() const { return _ptr; }
        };
    };

    template<typename T1, typename T2, typename std::enable_if_t<std::conjunction_v<is_data_type_t<T1>, is_data_type_t<T2>>, int> = 0>
    bool operator==(const matrix_expr<T1> &a, const matrix_expr<T2> &b) {
        if (a.get_size() != b.get_size())
            return false;
        auto[n, m] = a.get_size();
        for (size_t i = 0u; i < n; ++i) {
            for (size_t j = 0u; j < m; ++j) {
                if (a(i, j) != b(i, j))
                    return false;
            }
        }
        return true;
    }

    template<typename E1, typename E2, typename func>
    class matrix_arith_expr
            : public matrix_expr<typename std::common_type_t<typename E1::data_t, typename E2::data_t>> {
        const E1 &_u;
        const E2 &_v;
        const func _functor;

        using data_t = typename std::common_type_t<typename E1::data_t, typename E2::data_t>;
    public:
        matrix_arith_expr(const E1 &u, const E2 &v) : _u(u), _v(v), _functor() {
            static_assert(matrix_expr<typename E1::data_t>::template is_derived_from_matrix_expr_v<E1>);
            static_assert(matrix_expr<typename E2::data_t>::template is_derived_from_matrix_expr_v<E2>);
            assert(u.get_raw_size() == v.get_raw_size());
        }

        virtual data_t operator()(size_t i, size_t j) const { return _functor(_u(i, j), _v(i, j)); }

        dim_t get_size() const { return _v.get_size(); }

        constexpr size_t get_raw_size() const { return _v.get_raw_size(); }
    };

    template<typename E1, typename E2>
    using matrix_sum_expr = matrix_arith_expr<E1, E2, std::plus<>>;

    template<typename E1, typename E2>
    using matrix_sub_expr = matrix_arith_expr<E1, E2, std::minus<>>;

    template<typename E1, typename E2, typename std::enable_if_t<std::conjunction_v<
            typename matrix_expr<typename E1::data_t>::template is_derived_from_matrix_expr_t<E1>,
            typename matrix_expr<typename E2::data_t>::template is_derived_from_matrix_expr_t<E2>>, int> = 0>
    auto operator+(const E1 &u, const E2 &v) {
        return matrix_sum_expr<E1, E2>(u, v);
    }

    template<typename E1, typename E2, typename std::enable_if_t<std::conjunction_v<
            typename matrix_expr<typename E1::data_t>::template is_derived_from_matrix_expr_t<E1>,
            typename matrix_expr<typename E2::data_t>::template is_derived_from_matrix_expr_t<E2>>, int> = 0>
    auto operator-(const E1 &u, const E2 &v) {
        return matrix_sub_expr<E1, E2>(u, v);
    }

    template<typename E1, typename E2>
    class matrix_mult_expr : public matrix_expr<typename std::common_type_t<typename E1::data_t, typename E2::data_t>> {
        const E1 &_u;
        const E2 &_v;

        using base = matrix_mult_expr<E1, E2>;
        using data_t = typename std::common_type_t<typename E1::data_t, typename E2::data_t>;

        const dim_t _size;
        const size_t _m;
    public:
        matrix_mult_expr(const E1 &u, const E2 &v) : _u(u), _v(v), _size{u.get_size().n, v.get_size().m},
                                                     _m(u.get_size().m) {
            static_assert(matrix_expr<typename E1::data_t>::template is_derived_from_matrix_expr_v<E1>);
            static_assert(matrix_expr<typename E2::data_t>::template is_derived_from_matrix_expr_v<E2>);
            const auto[u_n, u_m] = u.get_size();
            const auto[v_n, v_m] = v.get_size();
            assert(u_m == v_n);
        }

        constexpr data_t operator()(size_t i, size_t j) const {
            auto res = data_t(0);
            for (size_t k = 0u; k < _m; ++k) {
                res += _u(i, k) * _v(k, j);
            }
            return res;
        }

        dim_t get_size() const { return _size; }

        constexpr size_t get_raw_size() const { return _size.n * _size.m; }
    };

    template<typename E1, typename E2, typename std::enable_if_t<std::conjunction_v<
            typename matrix_expr<typename E1::data_t>::template is_derived_from_matrix_expr_t<E1>,
            typename matrix_expr<typename E2::data_t>::template is_derived_from_matrix_expr_t<E2>>, int> = 0>
    auto operator*(const E1 &u, const E2 &v) {
        return matrix_mult_expr<E1, E2>(u, v);
    }

    template<typename scalar_t, typename matrix_t, typename func>
    class matrix_scalar_expr : public matrix_expr<typename std::common_type_t<scalar_t, typename matrix_t::data_t>> {
        const scalar_t &_scalar;
        const matrix_t &_v;
        const func _functor;

        using base = matrix_scalar_expr<scalar_t, matrix_t, func>;
        using data_t = typename std::common_type_t<scalar_t, typename matrix_t::data_t>;
    public:
        matrix_scalar_expr(const scalar_t &scalar, const matrix_t &v) : _scalar(scalar), _v(v), _functor() {
            static_assert(matrix_expr<typename matrix_t::data_t>::template is_derived_from_matrix_expr_v<matrix_t>);
        }

        constexpr data_t operator()(size_t i, size_t j) const { return _functor(_v(i, j), _scalar); }

        dim_t get_size() const { return _v.get_size(); }

        constexpr size_t get_raw_size() const { return _v.get_raw_size(); }
    };

    template<typename scalar_t, typename matrix_t, typename std::enable_if_t<std::conjunction_v<
            is_data_type_t<scalar_t>,
            typename matrix_expr<typename matrix_t::data_t>::template is_derived_from_matrix_expr_t<matrix_t>>, int> = 0>
    auto operator*(const scalar_t &k, const matrix_t& m) {
        return matrix_scalar_expr<scalar_t, matrix_t, std::multiplies<>>(k, m);
    }

    template<typename scalar_t, typename matrix_t, typename std::enable_if_t<std::conjunction_v<
            is_data_type_t<scalar_t>,
            typename matrix_expr<typename matrix_t::data_t>::template is_derived_from_matrix_expr_t<matrix_t>>, int> = 0>
    auto operator*(const matrix_t& m, const scalar_t &k) {
        return matrix_scalar_expr<scalar_t, matrix_t, std::multiplies<>>(k, m);
    }

    template<typename scalar_t, typename matrix_t, typename std::enable_if_t<std::conjunction_v<
            is_data_type_t<scalar_t>,
            typename matrix_expr<typename matrix_t::data_t>::template is_derived_from_matrix_expr_t<matrix_t>>, int> = 0>
    auto operator/(const matrix_t& m, const scalar_t &k) {
        return matrix_scalar_expr<scalar_t, matrix_t, std::divides<>>(k, m);;
    }

    template<typename matrix_t>
    class matrix_transpose_expr : public matrix_expr<typename matrix_t::data_t> {
        using data_t = typename matrix_t::data_t;

        const matrix_t &_matrix;
    public:
        matrix_transpose_expr(const matrix_t &matrix) : _matrix(matrix) {
            static_assert(
                    std::is_convertible_v<matrix_t, matrix_expr<typename matrix_t::data_t >>);
        }

        constexpr const data_t &operator()(size_t i, size_t j) const { return _matrix(j, i); }

        constexpr data_t &operator()(size_t i, size_t j) { return _matrix(j, i); }

        const dim_t &get_size() const {
            auto[n, m] = _matrix.get_size();
            return {m, n};
        }

        constexpr size_t get_raw_size() const { return _matrix.get_raw_size(); }
    };

}
#endif //LINALG_BASIC_MATRIX_EXPR_H
