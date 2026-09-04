#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>

std::vector < std::unordered_map < uint64_t, std::vector<int>>> build_index_parallel_hash(const std::string& genome, int k, int threads_number);

std::vector<std::vector<int>>build_index_parallel_array(const std::string& genome, int k, int threads_number);
