#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <cstdint>
#include "K_mers_search.h"
#include "Validation.h"
#include "Hit.h"
#include "Rolling_hash.h"
#include <unordered_set>

std::pair<std::vector<int>, size_t> collect_candidates_serial(const std::vector<std::vector<int>>& array_index,const std::unordered_map<uint64_t, std::vector<int>>& hash_index,const std::string& guide,int max_mismatch,int k, bool use_array)
{
    auto seeds = split_guide_into_seeds(guide, max_mismatch, k);

    std::unordered_set<int> candidate_set;

    for (const auto& seed : seeds) {
        uint64_t seed_code = encode_kmer(guide, seed.start, seed.length);

        if (use_array) {
            if (seed_code < array_index.size()) {
                for (int pos : array_index[seed_code]) {
                    candidate_set.insert(pos - seed.start);
                }
            }
        }
        else {
            auto it = hash_index.find(seed_code);
            if (it != hash_index.end()) {
                for (int pos : it->second) {
                    candidate_set.insert(pos - seed.start);
                }
            }
        }
    }
    std::vector<int> candidates(candidate_set.begin(), candidate_set.end());
    return { std::move(candidates), candidates.size() };
}

std::pair<std::vector<int>, size_t> collect_candidates_array(const std::vector<std::vector<int>>& all_index, const std::string& guide ,int max_mismatch, int k, int threads)
{
    auto seeds = split_guide_into_seeds(guide, max_mismatch, k);
    std::unordered_set<int> candidate_set;

    for (const auto& seed : seeds) {
        uint64_t seed_code = encode_kmer(guide, seed.start, seed.length);

        if (seed_code < all_index.size()) {
            const auto& positions = all_index[seed_code];
            for (int pos : positions) {
                candidate_set.insert(pos - seed.start);
            }
        }
    }
    std::vector<int> candidates(candidate_set.begin(), candidate_set.end());
    return { std::move(candidates), candidates.size() };
}


std::pair<std::vector<int>, size_t> collect_candidates_hash(const std::vector<std::unordered_map<uint64_t, std::vector<int>>>& all_index,const std::string& guide, int max_mismatch, int k, int threads)
{
    auto seeds = split_guide_into_seeds(guide, max_mismatch, k);

    std::unordered_set<int> candidate_set;

    for (const auto& seed : seeds) {
        uint64_t seed_code = encode_kmer(guide, seed.start, seed.length);

        for (int tid = 0; tid < threads; tid++) {
            auto it = all_index[tid].find(seed_code);
            if (it != all_index[tid].end()) {
                for (int pos : it->second) {
                    candidate_set.insert(pos - seed.start);
                }
            }
        }
    }
    std::vector<int> candidates(candidate_set.begin(), candidate_set.end());
    return { std::move(candidates), candidates.size() };
}



static std::vector<Hit> search_core(const std::vector<int>& candidates, const std::string& genome, const std::string& guide, int max_mismatch)
{
    if (!validateSearchInput(genome, guide, max_mismatch))return {};

    std::vector<Hit> hits;
    int guide_len = guide.length();
    hits.reserve(candidates.size());

    for (int pos : candidates) {
        if (pos + guide_len > (int)genome.length()) continue;

        int mismatches = 0;
        for (int j = 0; j < guide_len; j++) {
            if (genome[pos + j] != guide[j]) {
                mismatches++;
                if (mismatches > max_mismatch) break;
            }
        }
        if (mismatches <= max_mismatch) {
            hits.push_back({ pos, mismatches });
        }
    }
    return hits;
}



std::vector<Hit> search_serial(const std::vector<int>& candidates, const std::string& genome, const std::string& guide, int max_mismatch)
{
    if (candidates.empty())return {};
    return search_core(candidates, genome, guide, max_mismatch);
}