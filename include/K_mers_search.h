
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <cstdint>
#include "Validation.h"
#include "Hit.h"
#include "Rolling_hash.h"
#include <unordered_set>

std::vector<Hit> search_serial(const std::vector<int>& candidates, const std::string& genome, const std::string& guide, int max_mismatch);

std::pair<std::vector<int>, size_t> collect_candidates_serial(const std::vector<std::vector<int>>& array_index, const std::unordered_map<uint64_t, std::vector<int>>& hash_index, const std::string& guide, int max_mismatch, int k, bool use_array);

std::pair<std::vector<int>, size_t> collect_candidates_array(const std::vector<std::vector<int>>& all_index, const std::string& guide, int max_mismatch, int k, int threads);

std::pair<std::vector<int>, size_t> collect_candidates_hash(const std::vector<std::unordered_map<uint64_t, std::vector<int>>>& all_index, const std::string& guide, int max_mismatch, int k, int threads);
 
