#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>

struct Seed {
    int start;
    int length;
};

uint64_t encode_kmer(const std::string& s, size_t start_i, int k);

uint64_t rolling_hash_next(uint64_t current, char next_base, int k);

std::vector<Seed> split_guide_into_seeds(const std::string& guide, int max_mismatch, int min_seed_len);