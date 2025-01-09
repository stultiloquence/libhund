# Build
cmake_minimum_required(VERSION 3.12)

project(libhund VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)

add_executable(runhund runhund.cpp)
target_include_directories(mt-kahypar/include/)


add_library(libhund INTERFACE)
target_link_libraries(libhund INTERFACE -Llib/mt-kahypar/build/lib/libmtkahypar.so)
target_link_libraries(libhund INTERFACE -lmtkahypar)

target_compile_features(libhund INTERFACE cxx_std_17)

# target_include_directories(
#     libhund
#     INTERFACE
#     $<BUILD_INTERFACE:${libhund_SOURCE_DIR}>
#     $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
# )

# Test
add_subdirectory(lib/Catch2)
add_executable(libhund_test test/test.cpp)

target_link_libraries(libhund_test PRIVATE -Llib/mt-kahypar/build/lib/libmtkahypar.so)
target_link_libraries(libhund_test PRIVATE Catch2::Catch2WithMain)
target_include_directories(libhund_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
include(CTest)
include(Catch)
catch_discover_tests(libhund_test)
enable_testing()