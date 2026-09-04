#include <iostream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <string>
#include "K_mers_core.h"
#include "K_mers_search.h"
#include "Sliding_window.h"
#include "Rolling_hash.h"


static void test_search_serial_array();
static void test_search_serial_hash();

static void  test_search_parallel_array();
static void test_search_parallel_hash();

static void test_search_no_match();

static void test_empty_genome();
static void test_empty_guide();
static void test_exact_match();
static void test_large_mismatch();

static void test_seeding_theorem();


int main() {

    test_search_serial_array();
    test_search_serial_hash();

    test_search_parallel_array();
    test_search_parallel_hash();

    test_search_no_match();

    test_empty_genome();
    test_empty_guide();
    test_exact_match();
    test_large_mismatch();


    test_seeding_theorem();

    std::cout << "\n All search tests passed!\n";
	return 0;
}



static void test_search_serial_array() {
    std::string genome = "ATCGATCCACTG";
    std::string guide = "ATCG";
    int max_mismatch = 1; // k=2

    auto hits = k_mers_serial_array(genome, guide, max_mismatch).hits;
    auto check_hits = slidingWindow_serial(genome, guide, max_mismatch).hits;

    assert(!hits.empty() && "Serial Array: no hits found!");
    assert(hits.size() == check_hits.size() && "Serial Array: wrong number of hits!");

    for (const auto& hit : hits) {
        assert(hit.position >= 0 && hit.position < (int)genome.length() - guide.length() + 1);
        assert(hit.mismatches <= max_mismatch);
    }

    genome = "ATCGATCCACTGATCGAACCACTGATCGATGTACTGATCGATCCACTGATCGAACCACTGATCGAACCACTG";
    guide = "ATCGATCCACTGATCGAACCACTG";
    max_mismatch = 1; //k=12

    hits = k_mers_serial_array(genome, guide, max_mismatch).hits;
    check_hits = slidingWindow_serial(genome, guide, max_mismatch).hits;

    assert(hits.empty() && "Serial array: no hits found!"); // k=12 is not aloud

    std::cout << "test_search_serial_array passed!\n";
}


static void test_search_serial_hash() {
    std::string genome = "ATCGATCCACTGATCGAACCACTGATCGATGTACTGATCGATCCACTGATCGAACCACTGATCGAACCACTG";
    std::string guide = "ATCGATCCACTGATCGAACCACTG";
    int max_mismatch = 1; //k=12

    auto hits = k_mers_serial_hash(genome, guide, max_mismatch).hits;
    auto check_hits = slidingWindow_serial(genome,guide , max_mismatch).hits;
    
    assert(!hits.empty() && "Serial Hash: no hits found!");
    assert(hits.size() == check_hits.size() && "Serial Hash: wrong number of hits!");
    for (const auto& hit : hits) {
        assert(hit.position >= 0 && hit.position < (int)genome.length() - guide.length() + 1);
        assert(hit.mismatches <= max_mismatch);
    }
    std::cout << "test_search_serial passed!\n";
}


static void test_search_parallel_array() {
    std::string genome = "ATCGATCGATCG";
    std::string guide = "ATCG";
    int max_mismatch = 0; //k=3
    int threads = 4;

    auto hits = k_mers_parallel_array(genome, guide, max_mismatch, threads).hits;
    auto check_hits = slidingWindow_serial(genome, guide, max_mismatch).hits;

    assert(!hits.empty() && "Parallel Array: no hits found!");
    assert(hits.size() == check_hits.size() && "Parallel Array: wrong number of hits!");

    for (const auto& hit : hits) {
        assert(hit.position >= 0 && hit.position < (int)genome.length() - guide.length() + 1);
        assert(hit.mismatches <= max_mismatch);
    }

    std::cout << "test_search_parallel_array passed!\n";
}


static void test_search_parallel_hash() {
    std::string genome = "ATCGATCGACCG";
    std::string guide = "ATCG";
    int max_mismatch = 1; //k=2
    int threads = 4;

    auto hits = k_mers_parallel_hash(genome, guide, max_mismatch, threads).hits;
    auto check_hits = slidingWindow_serial(genome, guide, max_mismatch).hits;

    assert(!hits.empty() && "Parallel Hash: no hits found!");
    assert(hits.size() == check_hits.size() && "Parallel Hash: wrong number of hits!");

    for (const auto& hit : hits) {
        assert(hit.position >= 0 && hit.position < (int)genome.length() - guide.length() + 1);
        assert(hit.mismatches <= max_mismatch);
    }

    std::cout << "test_search_parallel_hash passed!\n";
}


static void test_search_no_match() {
    std::string genome = "ATCGATCGATCG";
    std::string guide = "GGGG"; 
    int max_mismatch = 0;
    int threads = 4;

    auto hits = k_mers_parallel_array(genome, guide, max_mismatch, threads).hits;

    assert(hits.empty() && "No Match: should return empty!");

    std::cout << "test_search_no_match passed!\n";
}

static void test_empty_genome() {
    std::string genome = "";
    std::string guide = "ATCG";
    int max_mismatch = 0;
    int threads = 4;

    auto hits1 = k_mers_serial_array(genome, guide, max_mismatch).hits;
    assert(hits1.empty() && "Serial: empty genome!");

    auto hits2 = k_mers_parallel_array(genome, guide, max_mismatch, threads).hits;
    assert(hits2.empty() && "Parallel Array: empty genome!");

    auto hits3 = k_mers_parallel_hash(genome, guide, max_mismatch, threads).hits;
    assert(hits3.empty() && "Parallel Hash: empty genome!");

    std::cout << "test_empty_genome passed!\n";
}

static void test_empty_guide() {
    std::string genome = "ATCGATCG";
    std::string guide = "";
    int max_mismatch = 0;
    int threads = 4;

    auto hits1 = k_mers_serial_array(genome, guide, max_mismatch).hits;
    assert(hits1.empty() && "Serial: empty guide!");

    auto hits2 = k_mers_parallel_array(genome, guide, max_mismatch, threads).hits;
    assert(hits2.empty() && "Parallel Array: empty guide!");

    auto hits3 = k_mers_parallel_hash(genome, guide, max_mismatch, threads).hits;
    assert(hits3.empty() && "Parallel Hash: empty guide!"); 

    std::cout << "test_empty_guide passed!\n";
}

static void test_exact_match() {
    std::string genome = "ATCGATCG";
    std::string guide = "ATCG";
    int max_mismatch = 0;
    int threads = 4;

    auto hits1 = k_mers_serial_array(genome, guide, max_mismatch).hits;
    assert(!hits1.empty() && "Serial: exact match should find hits!");

    auto hits2 = k_mers_parallel_array(genome, guide, max_mismatch, threads).hits;
    assert(!hits2.empty() && "Parallel Array: exact match should find hits!");

    auto hits3 = k_mers_parallel_hash(genome, guide, max_mismatch, threads).hits;
    assert(!hits3.empty() && "Parallel Hash: exact match should find hits!"); 

    std::cout << "test_exact_match passed!\n";
}

static void test_large_mismatch() {
    std::string genome = "ATCGATCG";
    std::string guide = "ATCG";
    int max_mismatch = 100;
    int threads = 4;

    auto hits1 = k_mers_serial_array(genome, guide, max_mismatch).hits;
    assert(hits1.empty() && "Serial: large mismatch should be empty!");

    auto hits2 = k_mers_parallel_array(genome, guide, max_mismatch, threads).hits;
    assert(hits2.empty() && "Parallel Array: large mismatch should be empty!");

    auto hits3 = k_mers_parallel_hash(genome, guide, max_mismatch, threads).hits;
    assert(hits3.empty() && "Parallel Hash: large mismatch should be empty!"); 

    std::cout << "test_large_mismatch passed!\n";
}

static void test_seeding_theorem() {
    std::string genome = "ACTGA";
    std::string guide = "ATCGA";
    int max_mismatch = 2;
    int threads = 4;

    auto hits1 = k_mers_serial_array(genome, guide, max_mismatch).hits;
   // assert(hits1.empty() && "Serial: large mismatch should be empty!");

    auto hits2 = k_mers_parallel_array(genome, guide, max_mismatch, threads).hits;
  //  assert(hits2.empty() && "Parallel Array: large mismatch should be empty!");

    auto hits3 = k_mers_parallel_hash(genome, guide, max_mismatch, threads).hits;
  //  assert(hits3.empty() && "Parallel Hash: large mismatch should be empty!"); 

    std::cout << "test_large_mismatch passed!\n";
}

