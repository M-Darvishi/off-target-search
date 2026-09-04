#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <cstdint>
#include "K_mers_index_builder_parallel.h"
#include "Validation.h"
#include "Rolling_hash.h"


std::vector < std::unordered_map < uint64_t, std::vector<int>>> build_index_parallel_hash(const std::string& genome, int k, int threads_number)
{
    if (!validateIndexInput(genome, k)) return {};

    size_t n = genome.size() - k + 1;
    uint64_t mask = (1ULL << (2 * (k - 1))) - 1;

    if (n == 0) return {};

    std::vector<std::unordered_map<uint64_t, std::vector<int>>> all_index(threads_number);

    omp_set_num_threads(threads_number);

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& local = all_index[tid];

        size_t chunk_size = (n + threads_number - 1) / threads_number;
        size_t start_i = tid * chunk_size;
        size_t end_i = std::min(start_i + chunk_size, n);

        if (start_i < n) {

            local.reserve(end_i - end_i);
            uint64_t code = encode_kmer(genome, start_i, k);
            local[code].emplace_back(start_i);
            for (size_t i = start_i+1 ; i < end_i; i++) {
                code = rolling_hash_next(code ,genome[i+k-1], k);
                local[code].emplace_back(i);
            }
        }
    }
    return all_index;
}


std::vector<std::vector<int>>build_index_parallel_array(const std::string& genome,int k, int threads_number)
{
    if (!validateIndexInput(genome, k)) return {};

    if (k >= 11) {
        std::cerr << "\n\n ERROR: Array method cannot be used with k >= 11!\n";
        std::cerr << "   k = " << k << " (>= 11)\n";
        std::cerr << "   Please use Hash or Hybrid method for k >= 11.\n\n";
        return {};
    }

    size_t table_size = 1ULL << (2 * k);
    std::vector<std::vector<int>> index(table_size);

    size_t n = genome.size() - k + 1;
    if (n == 0) return index;

    omp_set_num_threads(threads_number);

#pragma omp parallel
    {
        std::vector<std::vector<int>> local_index(table_size);

        int tid = omp_get_thread_num();

        size_t chunk_size = (n + threads_number - 1) / threads_number;
        size_t start = tid * chunk_size;
        size_t end   = std::min(start + chunk_size, n);

        if (start < end)
        { 
            uint64_t code = encode_kmer(genome, start, k);
            local_index[code].push_back((int)start);
             
            for (size_t i = start + 1; i < end; i++)
            {
                code = rolling_hash_next(code, genome[i + k - 1], k);
                local_index[code].push_back((int)i);
            }
        }

#pragma omp critical
        {
            for (size_t code = 0; code < table_size; code++)
            {
                auto& src = local_index[code];
                if (!src.empty())
                {
                    auto& dst = index[code];
                    dst.insert(dst.end(), src.begin(), src.end());
                }
            }
        }
    }

    return index;
}
