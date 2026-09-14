#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <cstdint>
#include <chrono>
#include <cmath>
#include "K_mers_core.h"
#include "K_mers_index_builder_serial.h"
#include "K_mers_index_builder_parallel.h"
#include "K_mers_search.h"
#include "Rolling_hash.h"


static std::vector<std::vector<int>> cached_array_index;
static std::vector<std::unordered_map<uint64_t, std::vector<int>>> cached_hash_index;

static std::vector<std::vector<int>>cached_array_index_serial;
static std::unordered_map< uint64_t, std::vector<int> > cached_hash_index_serial;

static std::string cached_genome;
static int cached_k = 0;
static int cached_threads = 0;
static Method cached_method = Method::Array;

void clear_kmer_cache() {
    cached_array_index.clear();
    cached_hash_index.clear();

    cached_array_index_serial.clear();
    cached_hash_index_serial.clear();

    cached_genome.clear();
    cached_k = 0;
    cached_threads = 0;

    cached_array_index.shrink_to_fit();
    cached_array_index_serial.shrink_to_fit();
    cached_genome.shrink_to_fit();

    cached_method = Method::Array;
}


static void build_cache_if_needed_serial(const std::string& genome, int k , Method method)
{
    if (cached_genome == genome && cached_k == k && cached_method == method) {
        return;
    }
    clear_kmer_cache();
    cached_method = method;

    if (method == Method::Array) 
    {
        if (!cached_array_index_serial.empty()) return;
        cached_array_index_serial = build_index_serial_array(genome, k);
    }
    else if(method == Method::Hash){
        if (!cached_hash_index_serial.empty()) return;
        cached_hash_index_serial = build_index_serial_hash(genome, k);
    }
    cached_genome = genome;
    cached_k = k;
}

static void build_cache_if_needed_parallel(const std::string& genome, int k, int threads , Method method)
{
    if (cached_genome == genome && cached_k == k && cached_threads == threads && cached_method == method) {
        return;
    }
    clear_kmer_cache();
    cached_method = method;

    switch (method) {
        case Method::Array: { 
            cached_array_index = build_index_parallel_array(genome, k, threads); 
            break;
        }
        case Method::Hash: { 
            cached_hash_index = build_index_parallel_hash(genome, k, threads); 
            break;
        }
    } 
    cached_genome = genome;
    cached_k = k;
    cached_threads = threads;
}

static  KmerHits kmer_switch(Method method,Mode mode, const std::string& genome, const std::string& guide, int max_mismatch, int threads)
{
    KmerHits result = {{},0,0,0,0};

    if (!validateSearchInput(genome, guide, max_mismatch)) return result ;
    int k = floor(guide.length() / float(max_mismatch + 1));
    if (k < 1) return result;

    if (method == Method::Array && k >= 11) {
        std::cerr << "\n\n ERROR: Array method cannot be used with k >= 11!\n";
        std::cerr << "   k = " << k << " (>= 11)\n";
        std::cerr << "   Please use Hash method for k >= 11.\n\n";
        return result;
    }

    std::pair<std::vector<int>, size_t> collected_candidates;


    auto start_building_index = std::chrono::high_resolution_clock::now();

    if (mode == Mode::Serial)
        build_cache_if_needed_serial(genome, k , method);  
    else
        build_cache_if_needed_parallel(genome, k, threads, method);
    
    auto end_building_index = std::chrono::high_resolution_clock::now();
    auto duration_building_index = std::chrono::duration_cast<std::chrono::milliseconds>(end_building_index - start_building_index);
    result.index_time = duration_building_index.count();


    auto start_collecting = std::chrono::high_resolution_clock::now();
    if (mode == Mode::Parallel) {
        switch (method) {
 
        case Method::Array: {
            collected_candidates = collect_candidates_array(cached_array_index, guide, max_mismatch, k, threads);
            break;
        }
        case Method::Hash: {
            collected_candidates = collect_candidates_hash(cached_hash_index, guide, max_mismatch, k, threads);
            break;
        }
        }
    }
    else {
        switch (method) {
        
        case Method::Array: {
            bool use_array = (method == Method::Array);
            collected_candidates = collect_candidates_serial(cached_array_index_serial, cached_hash_index_serial, guide, max_mismatch, k, use_array);
            break;
        }
        case Method::Hash: {
            bool use_array = (method == Method::Array);
            collected_candidates = collect_candidates_serial(cached_array_index_serial, cached_hash_index_serial, guide, max_mismatch, k, use_array);
            break;
        }
        }
    }

    auto end_collecting = std::chrono::high_resolution_clock::now();
    auto duration_collecting = std::chrono::duration_cast<std::chrono::milliseconds>(end_collecting - start_collecting);
    result.collect_time = duration_collecting.count();
    result.candidate_count = collected_candidates.second;


    auto start_searching = std::chrono::high_resolution_clock::now();

    result.hits = search_serial(collected_candidates.first , genome, guide, max_mismatch);

    auto end_searching = std::chrono::high_resolution_clock::now();
    auto duration_searching = std::chrono::duration_cast<std::chrono::milliseconds>(end_searching - start_searching);
    result.search_time= duration_searching.count();

    return result;
}



KmerHits k_mers_serial_hash(const std::string& genome, const std::string& guide, int max_mismatch)
{ 
    Mode mode = Mode::Serial;
    Method method = Method::Hash;
    return kmer_switch(method,mode, genome, guide, max_mismatch, 1);
}

KmerHits k_mers_serial_array(const std::string& genome, const std::string& guide, int max_mismatch)
{ 
    Mode mode = Mode::Serial;
    Method method = Method::Array;
    return kmer_switch(method, mode, genome, guide, max_mismatch, 1);
}


KmerHits k_mers_parallel_hash(const std::string& genome, const std::string& guide, int max_mismatch, int threads_number)
{ 
    Mode mode = Mode::Parallel;
    Method method = Method::Hash;
    return kmer_switch(method, mode, genome, guide, max_mismatch, threads_number);
}

KmerHits k_mers_parallel_array(const std::string& genome, const std::string& guide, int max_mismatch, int threads_number)
{
    Mode mode = Mode::Parallel;
    Method method = Method::Array;
    return kmer_switch(method, mode, genome, guide, max_mismatch, threads_number);
}

