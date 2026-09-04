cmake_minimum_required(VERSION 3.15)
project(ParallelOfftargetSearch)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(extern/pybind11)
find_package(OpenMP REQUIRED)

# CORE LIB
add_library(core STATIC
    src/Sliding_window.cpp
    src/K_mers_naive.cpp
    src/Validation.cpp
    src/GenomeLoader.cpp
    src/K_mers_index_builder_serial.cpp
    src/K_mers_search.cpp
    src/K_mers_index_builder_parallel.cpp
    src/K_mers_core.cpp
    src/Rolling_hash.cpp )
    
set_target_properties(core PROPERTIES POSITION_INDEPENDENT_CODE ON)


# Add source to this project's executable.
add_executable(SlidingWindowTest "tests/sliding_window_test.cpp")
target_link_libraries(SlidingWindowTest PRIVATE core)

add_executable(KmerNaiveTest "tests/k_mer_naive_test.cpp")
target_link_libraries(KmerNaiveTest PRIVATE core)

add_executable(KmerIndexTest "tests/k_mer_index_builder_test.cpp")
target_link_libraries(KmerIndexTest PRIVATE core)

add_executable(KmerSearchTest "tests/k_mer_search_test.cpp")
target_link_libraries(KmerSearchTest PRIVATE core)

add_executable(GenomeLoaderTest "tests/genomeLoader_test.cpp")
target_link_libraries(GenomeLoaderTest PRIVATE core)

add_executable(SlidingwindowBenchmark  "benchmarks/SlidingWindow_benchmark.cpp")
target_link_libraries(SlidingwindowBenchmark  PRIVATE core)

add_executable(KmerBenchmark "benchmarks/Kmer_benchmark.cpp")
target_link_libraries(KmerBenchmark PRIVATE core)

target_include_directories(SlidingWindowTest PRIVATE include)
target_include_directories(KmerNaiveTest PRIVATE include)
target_include_directories(KmerIndexTest PRIVATE include)
target_include_directories(KmerSearchTest PRIVATE include)
target_include_directories(GenomeLoaderTest PRIVATE include)
target_include_directories(SlidingwindowBenchmark PRIVATE include)
target_include_directories(KmerBenchmark PRIVATE include)


target_include_directories(core PUBLIC include)
target_link_libraries(core PUBLIC OpenMP::OpenMP_CXX)
    
# PYTHON MODULE
pybind11_add_module(offtarget_search  python/bindings.cpp)

#target_compile_options(offtarget_search PRIVATE -fPIC)
target_link_libraries(offtarget_search PRIVATE core pybind11::module)
set_target_properties(offtarget_search PROPERTIES PREFIX "" SUFFIX ".so")
