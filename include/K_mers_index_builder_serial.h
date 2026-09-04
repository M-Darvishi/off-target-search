#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>

std::unordered_map< uint64_t, std::vector<int> > build_index_serial_hash(const std::string& genome, int k);

std::vector<std::vector<int>> build_index_serial_array(const std::string& genome, int k);
