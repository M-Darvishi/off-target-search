#include <iostream>
#include <cassert>
#include <algorithm>
#include "K_mers_index_builder_serial.h"
#include "K_mers_index_builder_parallel.h"
#include "Rolling_hash.h"

static void test_serial_hash();
static void test_serial_array(); 

static void test_parallel_hash();
static void test_parallel_array();

static void test_consistency_serial();
static void test_consistency_parallel();

static void test_edge_cases();

static void test_random_genome();


int main() {

    test_serial_hash();
    test_serial_array();
   
    test_parallel_hash();
    test_parallel_array();

    test_consistency_serial();
    test_consistency_parallel();

    test_edge_cases();

    test_random_genome();

    std::cout << "\n All index builder tests passed!\n";
    return 0;
}


static void test_serial_hash() {
    std::string genome = "GCTAGCTAGC";
    int k = 2;
    std::vector<std::string> expected_kmers = { "GC", "CT", "TA", "AG", "GC", "CT", "TA", "AG", "GC" };
    std::vector<int> expected_positions = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };


    auto index = build_index_serial_hash(genome, k);

    int filled_count = 0;
    for (const auto& [code, positions] : index) {
        if (!positions.empty()) {
            filled_count++;
        }
    }
    assert(filled_count == 4 && "Serial Hash : un useable index filled!");

    for (size_t i = 0; i < expected_kmers.size(); i++) {
        uint64_t code = encode_kmer(expected_kmers[i], 0, expected_kmers[i].length());
        auto it = index.find(code);
        assert(it != index.end() && "Serial Hash : k-mer not found!");

        bool found = false;
        for (int pos : it->second) {
            if (pos == expected_positions[i]) {
                found = true;
                break;
            }
        }
        assert(found && "Serial Hash : Position mismatch!");
    }

    std::cout << "test_serial_hash passed!\n";
}


static void test_serial_array() {
    std::string genome = "GCTAGCTAGC";
    int k = 2;
    std::vector<std::string> expected_kmers = { "GC", "CT", "TA", "AG", "GC", "CT", "TA", "AG", "GC" };
    std::vector<int> expected_positions = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };

    auto index = build_index_serial_array(genome, k);   

    int filled_count = 0;
    for (size_t code = 0; code < index.size(); code++) {
        if (!index[code].empty()) {
            filled_count++;
        }
    }
    assert(filled_count == 4 && "Serial Array : un useable array index filled!");

    for (size_t i = 0; i < expected_kmers.size(); i++) {
        uint64_t code = encode_kmer(expected_kmers[i], 0, expected_kmers[i].length());
        assert(code < index.size() && !index[code].empty() && "Serial Array : k-mer not found!");

        bool found = false;
        for (int pos : index[code]) {
            if (pos == expected_positions[i]) {
                found = true;
                break;
            }
        }
        assert(found && "Serial Array : Position mismatch!");
    }

    std::cout << "test_serial_array passed!\n";
}


static void test_parallel_hash() {
    std::string genome = "ATCGATCG";
    int k = 2;
    std::vector<std::string> expected_kmers = { "AT", "TC", "CG", "GA", "AT", "TC", "CG" };
    std::vector<int> expected_positions = { 0, 1, 2, 3, 4, 5, 6 };

    int threads = 4;
    auto all_index = build_index_parallel_hash(genome, k, threads);

    std::unordered_map<uint64_t, std::vector<int>> merged;
    for (int tid = 0; tid < threads; tid++) {
        for (const auto& [code, positions] : all_index[tid]) {
            auto& vec = merged[code];
            vec.insert(vec.end(), positions.begin(), positions.end());
        }
    }

    int filled_count = 0;
    for (const auto& [code, positions] : merged) {
        if (!positions.empty()) {
            filled_count++;
        }
    }
    assert(filled_count == 4 && "Parallel Hash : un useable index filled!");

    for (size_t i = 0; i < expected_kmers.size(); i++) {
        uint64_t code = encode_kmer(expected_kmers[i], 0, expected_kmers[i].length());

        auto it = merged.find(code);
        assert(it != merged.end() && "Parallel Hash : k-mer not found in hash!");

        bool found = false;
        for (int pos : it->second) {
            if (pos == expected_positions[i]) {
                found = true;
                break;
            }
        }
        assert(found && "Parallel Hash : Position mismatch!");
    }
    std::cout << "test_parallel_hash passed!\n";
}


static void test_parallel_array() {
    std::string genome = "ATCGATCG";
    int k = 2;
    std::vector<std::string> expected_kmers = { "AT", "TC", "CG", "GA", "AT", "TC", "CG" };
    std::vector<int> expected_positions = { 0, 1, 2, 3, 4, 5, 6 };

    int threads = 4;
    auto index = build_index_parallel_array(genome, k, threads);

    int filled_count = 0;
    for (size_t code = 0; code < index.size(); code++) {
        if (!index[code].empty()) {
            filled_count++;
        }
    }
    assert(filled_count == 4 && "Parallel Array : un useable index filled!");
   
    for (size_t i = 0; i < expected_kmers.size(); i++) {
        uint64_t code = encode_kmer(expected_kmers[i], 0, expected_kmers[i].length());

        assert(code < index.size() && "Parallel Array : code out of range!");
        assert(!index[code].empty() && "Parallel Array : k-mer not found in parallel array!");

        bool found = false;
        for (int pos : index[code]) {
            if (pos == expected_positions[i]) {
                found = true;
                break;
            }
        }
        assert(found && "Parallel Array : Position mismatch!");
    }
    std::cout << "test_parallel_array passed!\n";
}


static void test_consistency_serial() {
    std::string genome = "ATCGATCG";
    int k = 2;

    auto array_index = build_index_serial_array(genome, k);
    auto hash_index = build_index_serial_hash(genome, k);

    for (size_t code = 0; code < array_index.size(); code++) {
        if (!array_index[code].empty()) {
            auto it = hash_index.find(code);
            assert(it != hash_index.end());
            assert(array_index[code] == it->second);
        }
    }

    std::cout << "test_consistency passed!\n";
}


static void test_consistency_parallel() {
    std::string genome = "ATCGATCG";
    int k = 2;
    int threads = 4;

    auto serial_array = build_index_serial_array(genome, k);

    auto parallel_array = build_index_parallel_array(genome, k, threads);

    auto parallel_hash = build_index_parallel_hash(genome, k, threads);
    std::unordered_map<uint64_t, std::vector<int>> hash_merged;
    for (int tid = 0; tid < threads; tid++) {
        for (const auto& [code, positions] : parallel_hash[tid]) {
            auto& vec = hash_merged[code];
            vec.insert(vec.end(), positions.begin(), positions.end());
        }
    }

    for (size_t code = 0; code < serial_array.size(); code++) {
        if (!serial_array[code].empty()) {

            assert(code < parallel_array.size() && "Parallel Array: code out of range!");
            assert(!parallel_array[code].empty() && "Parallel Array missing a key from serial!");
            assert(serial_array[code] == parallel_array[code] && "Parallel Array positions mismatch!");


            auto it_hash = hash_merged.find(code);
            assert(it_hash != hash_merged.end() && "Parallel Hash missing a key from serial!");
            assert(serial_array[code] == it_hash->second && "Parallel Hash positions mismatch!");
 
        }
    }

    std::cout << "test_consistency_parallel passed!\n";
}


static void test_random_genome() {

    std::string genome = "";
    for (int i = 0; i < 100; i++) {
        genome += "ATCG"[rand() % 4];
    }
    int k = 3;
    int threads = 4;

    auto idx1 = build_index_serial_array(genome, k);
    assert(!idx1.empty() && "serial_array : random genome failed!");

    auto idx2 = build_index_serial_hash(genome, k);
    assert(!idx2.empty() && "serial_hash : random genome failed!");

    auto idx3 = build_index_parallel_array(genome, k, threads);
    assert(!idx3.empty() && "parallel_array : random genome failed!");

    auto idx4 = build_index_parallel_hash(genome, k, threads);
    assert(!idx4.empty() && "parallel_hash : random genome failed!");
 

    std::cout << "test_random_genome passed!\n";
}


static void test_edge_cases() {

    auto idx1 = build_index_serial_array("", 2);
    assert(idx1.empty() && "serial_array : empty genome failed!");
    auto idx2 = build_index_serial_array("ATCG", 5);
    assert(idx2.empty() && "serial_array : unvalid k failed!");
    auto idx3 = build_index_serial_array("ATCG", 0);
    assert(idx3.empty() && "serial_array : empty k failed!");

    auto idx4 = build_index_serial_hash("", 2);
    assert(idx4.empty() && "serial_hash : empty genome failed!");
    auto idx5 = build_index_serial_hash("ATCG", 5);
    assert(idx5.empty() && "serial_hash : unvalid k failed!");
    auto idx6 = build_index_serial_hash("ATCG", 0);
    assert(idx6.empty() && "serial_array : empty k failed!");

    int threads = 4;
    auto idx7 = build_index_parallel_array("", 2, threads);
    assert(idx7.empty() && "parallel_array : empty genome failed!");
    auto idx8 = build_index_parallel_array("ATCG", 5, threads);
    assert(idx8.empty() && "parallel_array : unvalid k failed!");
    auto idx9 = build_index_parallel_array("ATCG", 0, threads);
    assert(idx9.empty() && "parallel_array : empty k failed!");

    auto idx10 = build_index_parallel_hash("", 2, threads);
    assert(idx10.empty() && "parallel_hash : empty genome failed!");
    auto idx11 = build_index_parallel_hash("ATCG", 5, threads);
    assert(idx11.empty() && "parallel_hash : unvalid k failed!");
    auto idx12 = build_index_parallel_hash("ATCG", 0, threads);
    assert(idx12.empty() && "parallel_hash : empty k failed!");

    std::cout << "test_edge_cases passed!\n";
}
