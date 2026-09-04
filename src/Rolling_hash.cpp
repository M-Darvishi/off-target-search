#include <iostream>
#include <string>
#include <vector>
#include "Rolling_hash.h"


uint64_t encode_kmer(const std::string& s, size_t start_i, int k) {
    uint64_t code = 0;
    for (int i = 0; i < k; i++) {
        code <<= 2;
        switch (s[start_i + i]) {
        case 'A': break;
        case 'C': code |= 1; break;
        case 'G': code |= 2; break;
        case 'T': code |= 3; break;
        }
    }
    return code;
}

static uint64_t base_to_bits(char c)
{
    switch (c)
    {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    }

    return 0;
}

uint64_t rolling_hash_next(uint64_t current, char next_base, int k)
{ 
    uint64_t mask = (1ULL << (2 * k)) - 1;

    current <<= 2;
    current &= mask;
    current |= base_to_bits(next_base);

    return current;
}

 std::vector<Seed> split_guide_into_seeds(const std::string& guide, int max_mismatch, int min_seed_len)
 {
     std::vector<Seed> seeds;
     int L = guide.length();
     int num_seeds = max_mismatch + 1;
      
     while (num_seeds > 1 && (L / num_seeds) < min_seed_len) {
         num_seeds--;
     }

     int seed_len = L / num_seeds;

     for (int i = 0; i < num_seeds; i++) {
         int start = i * seed_len;
         int end = (i == num_seeds - 1) ? L : start + seed_len;
         seeds.push_back({ start, end - start });
     } 
     return seeds;
 }