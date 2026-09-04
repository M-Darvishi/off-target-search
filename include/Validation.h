#pragma once
#include <string>
#include <vector>

bool validateSearchInput( const std::string& genome,   const std::string& guide, int max_mismatch);

bool validateIndexInput(const std::string& genome, int k);

bool validateKmerInput(const std::string& genome, const std::string& guide, int max_mismatch, int k);