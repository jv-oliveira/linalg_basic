//
// Created by joao on 19/02/2019.
//

#ifndef LINALG_BASIC_INDEX_EXCEPTION_H
#define LINALG_BASIC_INDEX_EXCEPTION_H

#include <stdexcept>
#include <algorithm>

namespace numeric {
    struct index_exception : public std::out_of_range {
        index_exception(std::string func, size_t index, size_t upper_limit)
                : std::out_of_range(func +
                                    ": index (" +
                                    std::to_string(index) +
                                    ") not in range [0," +
                                    std::to_string(std::clamp(upper_limit, size_t(1), SIZE_MAX) - size_t(1)) +
                                    "]."
        ) {}


    };
}


#endif //LINALG_BASIC_INDEX_EXCEPTION_H
