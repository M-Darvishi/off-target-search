#include <iostream>
#include <vector>
#include <omp.h>
#include "Hit.h"
#include <chrono>
#include "Validation.h"
#include "Sliding_window.h"

SlidingwindowHits slidingWindow_serial(std::string genome, std::string gRNA, int max_mismatch)
{
    SlidingwindowHits result = { {},0 };
    if (!validateSearchInput(genome, gRNA, max_mismatch) ) return {};

    int genome_len = genome.length();
    int guide_len = gRNA.length();

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i <= genome_len - guide_len; i++) {
        int mismatches = 0;
        for (int j = 0; j < guide_len; j++) {
            if (genome[i + j] != gRNA[j]) mismatches++;
            if (mismatches > max_mismatch) break;
        }
        if (mismatches <= max_mismatch) {
            result.hits.push_back({ i, mismatches });
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    result.search_time = duration.count();

    return result;
}


SlidingwindowHits slidingWindow_parallel(std::string genome, std::string gRNA, int max_mismatch , int threads_number)
{
    if (!validateSearchInput(genome, gRNA, max_mismatch)) return {};

    SlidingwindowHits result = { {},0 };
    int genome_len = genome.length();
    int guide_len = gRNA.length();
    result.hits.reserve(genome_len / 100);

    omp_set_num_threads(threads_number);

    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel 
        {
            std::vector<Hit> local_hits;
            local_hits.reserve(genome_len / 100 / omp_get_num_threads());

    #pragma omp for
            for (int i = 0; i <= genome_len - guide_len; i++) {
                int mismatches = 0;

                for (int j = 0; j < guide_len; j++) {
                    if (genome[i + j] != gRNA[j]) mismatches++;
                    if (mismatches > max_mismatch) break;
                }
                if (mismatches <= max_mismatch) {
                    local_hits.push_back({ i, mismatches });
                }
            }

    #pragma omp critical
            {
                result.hits.insert(result.hits.end(), local_hits.begin(), local_hits.end());
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        result.search_time = duration.count();

        return result;
}