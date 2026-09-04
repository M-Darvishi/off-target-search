#include <iostream>
#include <vector>
#include "Hit.h"

KmerHits k_mers_naive_serial(const std::string& genome, const std::string& guide, int max_mismatch, int k);

KmerHits k_mers_naive_parallel(const std::string& genome, const std::string& guide, int max_mismatch, int k, int threads_number);
