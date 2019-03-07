//
// Created by joao on 03/03/2019.
//

#ifndef LINALG_BASIC_VECTOR_H
#define LINALG_BASIC_VECTOR_H

#include "matrix_base.h"

namespace linalg {
    template<typename __data_t>
    class vector : public vector_base<__data_t> {
    public:
        using data_t = __data_t;
        using vector_base<data_t>::vector_base;

        constexpr vector(size_t n)
                : vector_base<data_t>(n, 1u) {}

        vector(std::initializer_list<data_t> groovy_vector) {
            this->_size = {std::size(groovy_vector), 1u};
            this->_storage.reset(new data_t[this->_size.n]);

            std::move(std::begin(groovy_vector), std::end(groovy_vector), &this->_storage[0]);
        }

        transposed_vector <data_t> transpose() const {
            return vector_base<data_t>::transpose();
        }
    };

    template<typename __data_t>
    class transposed_vector : public vector_base<__data_t> {
    public:
        using data_t = __data_t;
        using vector_base<data_t>::vector_base;

        constexpr transposed_vector(size_t n)
                : vector_base<data_t>(1u, n) {}

        transposed_vector(std::initializer_list<data_t> groovy_vector) {
            this->_size = {1u, std::size(groovy_vector)};
            this->_storage.reset(new data_t[this->_size.m]);
            std::move(std::begin(groovy_vector), std::end(groovy_vector), &this->_storage[0]);
        }

        vector<data_t> transpose() const {
            vector<data_t> v = vector_base<data_t>::transpose();
            return vector_base<data_t>::transpose();
        }
    };

    template<typename __data_t>
    class vector_view : public vector_base<__data_t> {
    public:
        using data_t = __data_t;
    private:
        using vector_base<data_t>::vector_base;

        friend class matrix_base<data_t>;

        size_t _pitch;

        vector_view(size_t n, std::shared_ptr<data_t[]> ptr, size_t pitch = 0u) :
                _pitch(pitch) {
            this->_size = {n, 1u};
            this->_storage = std::move(ptr);
        }

        virtual data_t get_element(size_t i, size_t j) const override {
            return this->_storage[i * _pitch + j];
        }

        virtual data_t &get_element(size_t i, size_t j) override {
            return this->_storage[i * _pitch + j];
        }
    };

    template<typename __data_t>
    class transposed_vector_view : public vector_base<__data_t> {
    public:
        using data_t = __data_t;
    private:
        using vector_base<data_t>::vector_base;
        friend class matrix_base<data_t>;

        transposed_vector_view(size_t n, std::shared_ptr<data_t[]> ptr) {
            this->_size = {1u, n},
            this->_storage = std::move(ptr);
        }

    };

    template<typename T>
    inline constexpr bool is_column_vector_v = is_one_of_v<T, vector<typename T::data_t>,
            vector_view<typename T::data_t>>;

    template<typename T>
    inline constexpr bool is_line_vector_v = is_one_of_v<T, transposed_vector<typename T::data_t>,
            transposed_vector_view<typename T::data_t>>;

    template<typename vec_t, typename std::enable_if_t<is_column_vector_v<vec_t> > = 0>
    auto operator*(const vec_t &u, const vec_t &v) {
        return matrix_mult_expr(u.transpose(), v);
    }

    template<typename vec_t, typename std::enable_if_t<is_line_vector_v<vec_t> > = 0>
    auto operator*(const vec_t &u, const vec_t &v) {
        return matrix_mult_expr(u, v.transpose());
    }

}

#endif //LINALG_BASIC_VECTOR_H
