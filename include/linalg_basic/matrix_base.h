//
// Created by joao on 19/02/2019.
//

#ifndef LINALG_BASIC_MATRIX_BASE_H
#define LINALG_BASIC_MATRIX_BASE_H

#include <utility>
#include <cstddef>
#include <cmath>
#include <complex>
#include <numeric>
#include <initializer_list>
#include <memory>
#include <algorithm>
#include <functional>
#include <string>
#include <ostream>
#include <cassert>


#include "forward_decl.h"
#include "matrix_expr.h"
#include "index_exception.h"


namespace numeric {

    template<typename __data_t>
    class matrix_base : public matrix_expr<matrix_base<__data_t>, __data_t> {
        template<typename, typename>
        friend
        class matrix_expr;

    public:
        using data_t = __data_t;
        using base = matrix_base<data_t>;
        using dim_t = typename matrix_expr<matrix_base<data_t>, data_t>::dim_t;
        using column_view = vector_view<data_t>;
        using line_view = transposed_vector_view<data_t>;
    private:

        template<typename T>
        bool is_same_storage_data(std::shared_ptr<T> storage) {
            for (int i = 0; i < this->get_raw_size(); ++i) {
                if (this->_storage[i] != storage[i])
                    return false;
            }
            return true;
        }

    protected:

        virtual data_t get_element(size_t i, size_t j) const {
            return _storage[i * _size.m + j];
        }

        virtual data_t &get_element(size_t i, size_t j) {
            return _storage[i * _size.m + j];
        }

        dim_t _size;
        std::shared_ptr<data_t[]> _storage;

        matrix_base() = default;

    public:
        matrix_base(const dim_t &size) : _size(size), _storage(new data_t[_size.n * _size.m]) {}

        template<typename T, typename std::enable_if_t<is_data_type_v<T>, int> = 0>
        matrix_base(std::initializer_list<std::initializer_list<T>> groovy_matrix)
                : _size {
                groovy_matrix.size(),
                groovy_matrix.begin()->size()
        },
                  _storage(new data_t[_size.n * _size.m]) {
            auto storage_ptr = _storage.get();
            for (const auto &line : groovy_matrix)
                storage_ptr = std::copy(line.begin(), line.end(), storage_ptr);
        }

        constexpr matrix_base(size_t n, size_t m)
                : matrix_base({n, m}) {}

        template<typename T, typename std::enable_if_t<matrix_expr<typename T::base, typename T::data_t>::template is_derived_from_matrix_expr_v<T>, int> = 0>
        matrix_base(const T &rhs)
                : _size(rhs._size),
                  _storage(get_raw_size()) {
            auto[n, m] = rhs.get_size();
            for (size_t i = 0u; i < n; ++i) {
                for (size_t j = 0u; j < m; ++j) {
                    get_element(i, j) = rhs(i, j);
                }
            }
        }

        template<typename T, typename std::enable_if_t<std::is_convertible_v<T, matrix_expr<typename T::base, typename T::data_t>>, int> = 0>
        matrix_base(T
                    &&rhs)
                :
                _size(std::move(rhs._size)
                ),

                _storage(get_raw_size()) {
            auto[n, m] = _size;
            for (size_t i = 0u; i < n; ++i) {
                for (size_t j = 0u; j < m; ++j) {
                    get_element(i, j) = std::move(rhs(i, j));
                }
            }
        }

        template<typename T, typename std::enable_if_t<is_data_type_v<T>> = 0>
        matrix_base(matrix_base<T>
                    &&rhs)
                :
                _size(std::move(rhs._size)
                ),
                _storage(std::move(rhs._storage)
                ) {
        }

        virtual ~matrix_base() = default;

        virtual dim_t get_size() const override {
            return _size;
        }

        virtual size_t get_raw_size() const override {
            return _size.n * _size.m;
        }

        column_view get_collum(size_t j) {
            if (j >= _size.m)
                throw index_exception(__PRETTY_FUNCTION__, j, _size.n);
            return column_view(_size.n, std::shared_ptr<data_t[]>(_storage, _storage.get() + j), _size.m);
        }

        line_view get_line(size_t i) {
            if (i >= _size.n)
                throw index_exception(__PRETTY_FUNCTION__, i, _size.n);
            return line_view(_size.m, std::shared_ptr<data_t[]>(_storage, _storage.get() + i * _size.m));
        }

        matrix_transpose_expr<matrix_base<data_t>> transpose() const {
            return matrix_transpose_expr(*this);
        }

        virtual data_t operator()(size_t i, size_t j) const override {
            return get_element(i, j);
        }

        virtual data_t &operator()(size_t i, size_t j) override {
            return get_element(i, j);
        }

        bool operator==(const matrix_base &rhs) const {
            return _size == rhs._size &&
                   is_same_storage_data(rhs._storage);
        }

        bool operator!=(const matrix_base &rhs) const {
            return not(rhs == *this);
        }
    };

    template<typename __data_t>
    class vector_base : public matrix_base<__data_t> {
    public:
        using data_t  = __data_t;
        using base = typename matrix_base<data_t>::base;
        using dim_t = typename matrix_base<data_t>::dim_t;
    protected:
        using matrix_base<data_t>::matrix_base;

        template<typename T, typename std::enable_if_t<is_data_type_v<std::is_scalar_v<T>>> = 0>

        vector_base(std::initializer_list<std::initializer_list<T>> groovy_matrix);

        virtual dim_t get_size() const override {
            return matrix_base<data_t>::get_size();
        }

        virtual size_t get_raw_size() const override {
            return matrix_base<data_t>::get_raw_size();
        }

        virtual data_t &operator()(size_t i, size_t j) override {
            return matrix_base<data_t>::operator()(i, j);
        }

        vector_base() = default;

    public:

        size_t length() const {
            return this->get_raw_size();
        }

        size_t norm() const {
            return std::sqrt(std::inner_product(this->begin(), this->end(), this->begin(), 0));
        }

        explicit operator data_t() const {
            assert(this->get_raw_size() == 1u);
            return this->get_element(0u, 0u);
        }

        data_t operator()(size_t i) const {
            return get_element(i % this->_size.n, i % this->_size.m);
        }

        const data_t &operator[](size_t i) const {
            return get_element(i % this->_size.n, i % this->_size.m);
        }

        data_t &operator()(size_t i) {
            return get_element(i % this->_size.n, i % this->_size.m);
        }

        data_t &operator[](size_t i) {
            return get_element(i % this->_size.n, i % this->_size.m);
        }
    };
}
#endif //LINALG_BASIC_MATRIX_BASE_H
