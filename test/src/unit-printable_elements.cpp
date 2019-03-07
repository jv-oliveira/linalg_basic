//
// Created by joao on 06/03/2019.
//

#include "linalg_basic/linalg_elements.h"

#include <catch2/catch.hpp>

using namespace linalg;

TEST_CASE("printable_elements") {
    SECTION("matrix") {
        auto m = matrix<int>({
                                     {1, 2, 3, 4},
                                     {5, 6, 7, 8}
                             });
        std::ostringstream oss;
        oss << m;
        CHECK(oss.str() == "linalg::matrix<int> of size {2, 4}\n"
                           "{\n"
                           "\t{1, 2, 3, 4},\n"
                           "\t{5, 6, 7, 8}\n"
                           "}");
    }

    SECTION("square_matrix") {
        auto m = square_matrix<int>({
                                            {0, 1, 2, 3, 4},
                                            {5, 6, 7, 8, 9},
                                            {0, 1, 2, 3, 4},
                                            {5, 6, 7, 8, 9},
                                            {0, 1, 2, 3, 4},
                                    });
        std::ostringstream oss;
        oss << m;
        CHECK(oss.str() == "linalg::square_matrix<int> of size {5, 5}\n"
                           "{\n"
                           "\t{0, 1, 2, 3, 4},\n"
                           "\t{5, 6, 7, 8, 9},\n"
                           "\t{0, 1, 2, 3, 4},\n"
                           "\t{5, 6, 7, 8, 9},\n"
                           "\t{0, 1, 2, 3, 4}\n"
                           "}");
    }

    SECTION("vector") {
        auto v = vector<int>{0, 1, 2, 3, 4};
        std::ostringstream oss;
        oss << v;
        CHECK(oss.str() == "linalg::vector<int> of size {5, 1}\n"
                           "{\n"
                           "\t{0},\n"
                           "\t{1},\n"
                           "\t{2},\n"
                           "\t{3},\n"
                           "\t{4}\n"
                           "}");
    }

    SECTION("transposed_vector") {
        auto v = transposed_vector<int>{1, 2, 3, 4};
        std::ostringstream oss;
        oss << v;
        CHECK(oss.str() == "linalg::transposed_vector<int> of size {1, 4}\n"
                           "{\n"
                           "\t{1, 2, 3, 4}\n"
                           "}");
    }


    SECTION("vector_view") {
        auto m = square_matrix<int>({
                                            {0, 1, 2, 3, 4},
                                            {5, 6, 7, 8, 9},
                                            {0, 1, 2, 3, 4},
                                            {5, 6, 7, 8, 9},
                                            {0, 1, 2, 3, 4},
                                    });
        auto v = m.get_collum(0);
        std::ostringstream oss;
        oss << v;
        CHECK(oss.str() == "linalg::vector_view<int> of size {5, 1}\n"
                           "{\n"
                           "\t{0},\n"
                           "\t{5},\n"
                           "\t{0},\n"
                           "\t{5},\n"
                           "\t{0}\n"
                           "}");
    }

    SECTION("transposed_vector_view") {
        auto m = square_matrix<int>({
                                            {0, 1, 2, 3, 4},
                                            {5, 6, 7, 8, 9},
                                            {0, 1, 2, 3, 4},
                                            {5, 6, 7, 8, 9},
                                            {0, 1, 2, 3, 4},
                                    });
        auto v = m.get_line(0);
        std::ostringstream oss;
        oss << v;
        CHECK(oss.str() == "linalg::transposed_vector_view<int> of size {1, 5}\n"
                           "{\n"
                           "\t{0, 1, 2, 3, 4}\n"
                           "}");
    }

}
