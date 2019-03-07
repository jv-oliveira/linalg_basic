//
// Created by joao on 07/03/2019.
//

#ifndef LINALG_BASIC_UTILS_H
#define LINALG_BASIC_UTILS_H

#include <type_traits>
#include <complex>

namespace linalg {

    template<typename T, typename ... Ts>
    using is_one_of = std::disjunction<std::is_same<T, Ts>...>;

    template<typename T, typename ... Ts>
    inline constexpr bool is_one_of_v = is_one_of<T, Ts ...>::value;

    template<class T>
    struct is_complex : std::false_type {
    };
    template<class T>
    struct is_complex<std::complex<T>> : std::true_type {
    };

    template<typename data_t>
    using is_data_type_t = std::disjunction<std::is_arithmetic<data_t>, is_complex<data_t>>;

    template<typename data_t>
    inline constexpr bool is_data_type_v = is_data_type_t<data_t>::value;




}

#endif //LINALG_BASIC_UTILS_H
