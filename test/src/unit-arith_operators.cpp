//
// Created by joao on 07/03/2019.
//

#include "linalg_basic/linalg_elements.h"

#include <catch2/catch.hpp>

using namespace linalg;

TEST_CASE("arithmetic_operators") {
    SECTION("matrix_operators") {
        auto a = matrix<int>{
                {1, 2, 3, 4},
                {5, 6, 7, 8}
        };

        SECTION("addition") {
            auto b = matrix<int>{
                    {0,  -2, -3, -4},
                    {-5, -6, -7, -8}
            };
            CHECK((a + b) == matrix<int>{
                    {1, 0, 0, 0},
                    {0, 0, 0, 0}
            });
        }

        SECTION("subtraction") {
            auto b = matrix<int>{
                    {0, 2, 3, 4},
                    {5, 6, 7, 8}
            };
            CHECK((a - b) == matrix<int>{
                    {1, 0, 0, 0},
                    {0, 0, 0, 0}
            });
        }

        SECTION("multiplication") {
            SECTION("by_scalar") {
                CHECK(2 * a == matrix<int>{
                        {2, 4, 6, 8},
                        {10, 12, 14, 16}
                });
            }

            SECTION("by_matrix") {
                auto b = matrix<int>{
                        {0, 2},
                        {1, 0},
                        {0, 1},
                        {2, 0}
                };
                CHECK( a * b == matrix<int>{
                        {10, 5},
                        {22, 17}
                });
                auto c = matrix<int> {
                        {2, -1},
                        {2, 1}
                };
                CHECK( c * a == matrix<int>{
                        {-3, -2, -1, 0},
                        {7, 10, 13, 16}
                });
            }
        }

        SECTION("division_by_scalar") {
            CHECK( a / 2 == matrix<int>{
                   {0, 1, 1, 2},
                   {2, 3, 3, 4}
            });
        }
    }
}