#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <cstdint>
#include "K_mers_index_builder_serial.h"
#include "Validation.h"
#include "Rolling_hash.h"

std::unordered_map< uint64_t, std::vector<int> > build_index_serial_hash(const std::string& genome, int k)
{
    if (!validateIndexInput(genome, k)) return {};

    std::unordered_map<uint64_t, std::vector<int>> index;
    size_t n = genome.length() - k + 1;
    index.reserve(n); 

    if (n == 0) return index;

    uint64_t code = encode_kmer(genome, 0, k);
    index[code].emplace_back(0);
    for (int i = 1; i < (int)n; i++) {
        code = rolling_hash_next(code, genome[i + k - 1], k);uint64_t code = encode_kmer(genome , i , k);
        index[code].emplace_back(i);
    }

    return index;
}

std::vector<std::vector<int>> build_index_serial_array(const std::string& genome, int k)
{
    if (!validateIndexInput(genome, k)) return {};

    if ( k >= 11) {
        std::cerr << "\n\n ERROR: Array method cannot be used with k >= 11!\n";
        std::cerr << "   k = " << k << " (>= 11)\n";
        std::cerr << "   Please use Hash or Hybrid method for k >= 11.\n\n";
        return {};
    }
    
    size_t size = 1ULL << (2 * k);
    std::vector<std::vector<int>> index(size);
    size_t n = genome.size() - k + 1;

    if (n == 0) return index;

    uint64_t code = encode_kmer(genome , 0 , k);
    index[code].emplace_back(0);
    for (int i = 1; i < (int)n; i++) {
        code = rolling_hash_next(code, genome[i+k-1], k);
        index[code].push_back(i);
    } 
    return index;
}
