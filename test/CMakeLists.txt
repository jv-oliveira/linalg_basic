option(LINALG_Sanitizer "Build test suite with Clang sanitizer" OFF)
option(LINALG_Valgrind "Execute test suite with Valgrind" OFF)
option(LINALG_NoExceptions "Build test suite without exceptions" OFF)
option(LINALG_Coverage "Build test suite with coverage information" OFF)

if (LINALG_Sanitizer)
    message(STATUS "Building test suite with Clang sanitizer")
    set(CMAKE_CXX_FLAGS "-g -O2 -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
endif ()

if (LINALG_Valgrind)
    find_program(CMAKE_MEMORYCHECK_COMMAND valgrind)
    message(STATUS "Executing test suite with Valgrind (${CMAKE_MEMORYCHECK_COMMAND})")
    set(memcheck_command "${CMAKE_MEMORYCHECK_COMMAND} ${CMAKE_MEMORYCHECK_COMMAND_OPTIONS} --error-exitcode=1 --leak-check=full")
    separate_arguments(memcheck_command)
endif ()

if (LINALG_NoExceptions)
    message(STATUS "Building test suite without exceptions")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLINALG_NOEXCEPTION")
    set(CATCH_TEST_FILTER -e)
endif ()

if (LINALG_Coverage)
    message(STATUS "Building test suite with coverage information")
    if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        message(FATAL_ERROR "LINALG_Coverage requires GCC.")
    endif ()

    # enable profiling
    set(CMAKE_CXX_FLAGS "--coverage -g -O0 -fprofile-arcs -ftest-coverage")

    # from https://github.com/RWTH-HPC/CMake-codecov/blob/master/cmake/FindGcov.cmake
    get_filename_component(COMPILER_PATH "${CMAKE_CXX_COMPILER}" PATH)
    string(REGEX MATCH "^[0-9]+" GCC_VERSION "${CMAKE_CXX_COMPILER_VERSION}")
    find_program(GCOV_BIN NAMES gcov-${GCC_VERSION} gcov HINTS ${COMPILER_PATH})

    # collect all source files from the chosen include dir
    file(GLOB_RECURSE SOURCE_FILES ${LINALG_INCLUDE_BUILD_DIR}*.h)

    # add target to collect coverage information and generate HTML file
    # (filter script from https://stackoverflow.com/a/43726240/266378)
    add_custom_target(lcov_html
            COMMAND lcov --directory . --capture --output-file LINALG.info --gcov-tool ${GCOV_BIN} --rc lcov_branch_coverage=1
            COMMAND lcov -e LINALG.info ${SOURCE_FILES} --output-file LINALG.info.filtered --rc lcov_branch_coverage=1
            COMMAND ${CMAKE_SOURCE_DIR}/test/thirdparty/imapdl/filterbr.py LINALG.info.filtered > LINALG.info.filtered.noexcept
            COMMAND genhtml --title "LINALG for Modern C++" --legend --demangle-cpp --output-directory html --show-details --branch-coverage LINALG.info.filtered.noexcept
            COMMENT "Generating HTML report test/html/index.html"
            )
endif ()

#############################################################################
# Catch library with the main function to speed up build
#############################################################################
find_package(Catch2 REQUIRED)

add_library(catch_main OBJECT
        "src/unit.cpp"
        )
target_compile_features(catch_main PUBLIC cxx_std_17)
target_link_libraries(catch_main Catch2::Catch2)

#############################################################################
# one executable for each unit test file
#############################################################################

file(GLOB files "src/unit-*.cpp")
foreach (file ${files})
    get_filename_component(file_basename ${file} NAME_WE)
    string(REGEX REPLACE "unit-([^$]+)" "test-\\1" testcase ${file_basename})

    add_executable(${testcase} $<TARGET_OBJECTS:catch_main> ${file})
    target_compile_definitions(${testcase} PRIVATE
            CATCH_CONFIG_FAST_COMPILE
            )
    target_compile_options(${testcase} PRIVATE
            $<$<CXX_COMPILER_ID:GNU>:-Wno-deprecated-declarations>
            )
    target_link_libraries(${testcase} ${LINALG_TARGET_NAME} Catch2::Catch2)

    add_test(NAME "${testcase}_default"
            COMMAND ${testcase} ${CATCH_TEST_FILTER}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            )
    set_tests_properties("${testcase}_default" PROPERTIES LABELS "default")

    add_test(NAME "${testcase}_all"
            COMMAND ${testcase} ${CATCH_TEST_FILTER} "*"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            )
    set_tests_properties("${testcase}_all" PROPERTIES LABELS "all")

    if (LINALG_Valgrind)
        add_test(NAME "${testcase}_valgrind"
                COMMAND ${memcheck_command} ${CMAKE_CURRENT_BINARY_DIR}/${testcase} ${CATCH_TEST_FILTER}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                )
        set_tests_properties("${testcase}_valgrind" PROPERTIES LABELS "valgrind")
    endif ()
endforeach ()