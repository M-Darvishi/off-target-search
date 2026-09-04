#include <iostream>
#include <vector>
#include "Hit.h"

enum class Method {
	Array,
	Hash
};

enum class Mode {
	Serial,
    Parallel
};


 
KmerHits k_mers_serial_array(const std::string& genome, const std::string& guide, int max_mismatch);

KmerHits k_mers_serial_hash(const std::string& genome, const std::string& guide, int max_mismatch);


KmerHits k_mers_parallel_array(const std::string& genome, const std::string& guide, int max_mismatch, int threads_number);

KmerHits k_mers_parallel_hash(const std::string& genome, const std::string& guide, int max_mismatch, int threads_number);


void clear_kmer_cache();