#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <chrono>
#include "K_mers_naive.h"
#include "Validation.h"

static std::unordered_map< std::string_view, std::vector<int> > build_index_serial(const std::string& genome, int k)
{
    if (!validateIndexInput(genome, k)) return {};

    std::unordered_map<std::string_view, std::vector<int>> index;

    for (int i = 0; i <= genome.length() - k; i++) {
        std::string_view kmer(&genome[i], k);
        index[kmer].push_back(i);
    }

    return index;
}

static std::unordered_map< std::string_view, std::vector<int> > build_index_parallel(const std::string& genome, int k, int threads_number)
{
    if (!validateIndexInput(genome, k)) return {};

    std::unordered_map<std::string_view, std::vector<int>> index;

    omp_set_num_threads(threads_number);

#pragma omp parallel
    {
        std::unordered_map<std::string_view, std::vector<int>> local_index;

#pragma omp for
        for (int i = 0; i <= (int)genome.length() - k; i++) {
            std::string_view kmer(&genome[i], k);
            local_index[kmer].push_back(i);
        }

#pragma omp critical
        {
            for (const auto& [kmer, positions] : local_index) {
                auto& target = index[kmer];
                target.insert(target.end(), positions.begin(), positions.end());
            }
        }
    }

    return index;
}

static std::vector<Hit> search_serial(const std::vector<int>& candidates, const std::string& genome, const std::string& guide, int max_mismatch) {
    std::vector<Hit> hits;
    int guide_len = guide.length();

    for (int pos : candidates) {
        if (pos + guide_len > (int)genome.length()) {
            continue;
        }

        int mismatches = 0;
        for (int j = 0; j < guide_len; j++) {
            if (genome[pos + j] != guide[j]) {
                mismatches++;
                if (mismatches > max_mismatch) {
                    break;
                }
            }
        }
        if (mismatches <= max_mismatch) {
            hits.push_back({ pos, mismatches });
        }
    }
    return hits;
} 


KmerHits k_mers_naive_serial(const std::string& genome,  const std::string& guide, int max_mismatch, int k)
{   
    KmerHits result = {{},0,0,0,0};

    if (!validateKmerInput(genome, guide, max_mismatch, k))
    {
        return {};
    }
    
    auto start_building_index = std::chrono::high_resolution_clock::now();

    std::unordered_map index = build_index_serial(genome, k);

    auto end_building_index = std::chrono::high_resolution_clock::now();
    auto duration_building_index = std::chrono::duration_cast<std::chrono::milliseconds>(end_building_index - start_building_index);
    result.index_time = duration_building_index.count();

    std::string first_part = guide.substr(0, k);
    auto it = index.find(first_part);
    if (it == index.end()) {
        return result;
    }

    const std::vector<int>& candidate_positions = it->second;
    result.candidate_count = candidate_positions.size();
    
    auto start_searching = std::chrono::high_resolution_clock::now();

    result.hits = search_serial(candidate_positions, genome, guide, max_mismatch);

    auto end_searching = std::chrono::high_resolution_clock::now();
    auto duration_searching = std::chrono::duration_cast<std::chrono::milliseconds>(end_searching - start_searching);
    result.search_time= duration_searching.count();

    return result;
}


KmerHits k_mers_naive_parallel(const std::string& genome, const std::string& guide, int max_mismatch, int k , int threads_number) 
{
    KmerHits result = {{},0,0,0,0};

    if (!validateKmerInput( genome, guide,  max_mismatch,  k))
    {
        return {};
    }
    
    auto start_building_index = std::chrono::high_resolution_clock::now();

    std::unordered_map index = build_index_parallel(genome, k , threads_number);

    auto end_building_index = std::chrono::high_resolution_clock::now();
    auto duration_building_index = std::chrono::duration_cast<std::chrono::milliseconds>(end_building_index - start_building_index);
    result.index_time = duration_building_index.count();

    std::string first_part = guide.substr(0, k);
    auto it = index.find(first_part);
    if (it == index.end()) {
        return result;
    }

    const std::vector<int>& candidate_positions = it->second;
    result.candidate_count = candidate_positions.size();
    
    auto start_searching = std::chrono::high_resolution_clock::now();

    result.hits = search_serial(candidate_positions, genome, guide, max_mismatch);

    auto end_searching = std::chrono::high_resolution_clock::now();
    auto duration_searching = std::chrono::duration_cast<std::chrono::milliseconds>(end_searching - start_searching);
    result.search_time= duration_searching.count();

    return result;
}
