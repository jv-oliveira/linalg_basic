cmake_minimum_required(VERSION 3.13)

##
## PROJECT
## name and version
##
project(linalg_basic VERSION 0.0.1 LANGUAGES CXX)

## ===== -------  ------- ===== ##
##
## Tune CMake
## C++ standard and output CMAKE_CXX_FLAGS
##
set(CMAKE_CXX_STANDARD 17)
message(INFO "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

##
## INCLUDE
##
##
include(ExternalProject)

##
## OPTIONS
##
option(LINALG_BuildTests "Build the unit tests when BUILD_TESTING is enabled." ON)

##
## CONFIGURATION
##
set(LINALG_TARGET_NAME ${PROJECT_NAME})
set(LINALG_CONFIG_INSTALL_DIR "lib/cmake/${PROJECT_NAME}"
        CACHE INTERNAL "")
set(LINALG_INCLUDE_INSTALL_DIR "include")
set(LINALG_TARGETS_EXPORT_NAME "${PROJECT_NAME}Targets")
set(LINALG_CMAKE_CONFIG_TEMPLATE "cmake/config.cmake.in")
set(LINALG_CMAKE_CONFIG_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(LINALG_CMAKE_VERSION_CONFIG_FILE "${LINALG_CMAKE_CONFIG_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
set(LINALG_CMAKE_PROJECT_CONFIG_FILE "${LINALG_CMAKE_CONFIG_DIR}/${PROJECT_NAME}Config.cmake")
set(LINALG_CMAKE_PROJECT_TARGETS_FILE "${LINALG_CMAKE_CONFIG_DIR}/${PROJECT_NAME}Targets.cmake")

set(LINALG_INCLUDE_BUILD_DIR "${PROJECT_SOURCE_DIR}/include/")
message(STATUS "Using the multi-header code from ${LINALG_INCLUDE_BUILD_DIR}")

##
## TARGET
## create target and add include path
##
add_library(${LINALG_TARGET_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${LINALG_TARGET_NAME} ALIAS ${LINALG_TARGET_NAME})

target_include_directories(
        ${LINALG_TARGET_NAME}
        INTERFACE
        $<BUILD_INTERFACE:${LINALG_INCLUDE_BUILD_DIR}>
        $<INSTALL_INTERFACE:include>
)

##
## TESTS
## create and configure the unit test target
##
include(CTest) #adds option BUILD_TESTING (default ON)

if (BUILD_TESTING AND LINALG_BuildTests)
    enable_testing()
    add_subdirectory(test)
endif ()

##
## INSTALL
## install header files, generate and install cmake config files for find_package()
##
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${LINALG_CMAKE_VERSION_CONFIG_FILE} COMPATIBILITY SameMajorVersion
)
configure_file(
        ${LINALG_CMAKE_CONFIG_TEMPLATE}
        ${LINALG_CMAKE_PROJECT_CONFIG_FILE}
        @ONLY
)

install(
        DIRECTORY ${LINALG_INCLUDE_BUILD_DIR}
        DESTINATION ${LINALG_INCLUDE_INSTALL_DIR}
)
install(
        FILES ${LINALG_CMAKE_PROJECT_CONFIG_FILE} ${LINALG_CMAKE_VERSION_CONFIG_FILE}
        DESTINATION ${LINALG_CONFIG_INSTALL_DIR}
)
if (NLOHMANN_ADD_NATVIS)
    install(
            FILES ${NLOHMANN_NATVIS_FILE}
            DESTINATION .
    )
endif ()
export(
        TARGETS ${LINALG_TARGET_NAME}
        NAMESPACE ${PROJECT_NAME}::
        FILE ${LINALG_CMAKE_PROJECT_TARGETS_FILE}
)
install(
        TARGETS ${LINALG_TARGET_NAME}
        EXPORT ${LINALG_TARGETS_EXPORT_NAME}
        INCLUDES DESTINATION ${LINALG_INCLUDE_INSTALL_DIR}
)
install(
        EXPORT ${LINALG_TARGETS_EXPORT_NAME}
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${LINALG_CONFIG_INSTALL_DIR}
)