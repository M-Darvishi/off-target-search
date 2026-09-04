#include <iostream>
#include <cassert>
#include <algorithm>
#include "K_mers_naive.h"


void simpleTest();
void completeMatchTest();
void noMatchTest();
void oneMatchTest();
void oneMismatchTest();
void invalidKTest();
void emptyGenomeTest();
void emptyGRNATest();

static void test_seeding_theorem() {
    std::string genome = "ACTGA";
    std::string gRNA = "ATCGA";
    int max_mismatch = 2;
    int threads = 4;


    auto serial = k_mers_naive_serial(genome, gRNA, max_mismatch, 4).hits;
    auto index_parallel = k_mers_naive_parallel(genome, gRNA, max_mismatch, 5, 16).hits;


    std::cout << serial.size() << " _ " << index_parallel.size() << " ;";

    std::cout << "test_large_mismatch passed!\n";
}


int main() {
    
    simpleTest();
    completeMatchTest();
    noMatchTest();
    oneMatchTest();
    oneMismatchTest();
    invalidKTest();
    emptyGenomeTest();
    emptyGRNATest();

    test_seeding_theorem();

	return 0;
}

void simpleTest() {
    std::string genome = "ATGCGTACGTAGCTAGCTAGCATGCAATGCGTAGCTAGCATGCA";
    std::string gRNA = "ATGCGTACGTAGCTAGCTA";
    int max_mismatch = 2;

    auto serial = k_mers_naive_serial(genome , gRNA , max_mismatch , 5).hits;
    auto index_parallel = k_mers_naive_parallel(genome, gRNA, max_mismatch, 5 , 16 ).hits;

    assert(serial.size() == index_parallel.size());

    for (size_t i = 0; i < serial.size(); i++) {
        assert(serial[i].position ==index_parallel[i].position);
        assert(serial[i].mismatches ==index_parallel[i].mismatches);

        assert(index_parallel[i].position == index_parallel[i].position );
        assert(index_parallel[i].mismatches == index_parallel[i].mismatches);
    }
    std::cout << "Sipmle Test Passed!\n";
}

void completeMatchTest() {

    std::string genome = "ATCGCATCG";
    std::string guide = "ATCG";

    int max_mismatch = 0;
    int k = 2;

    auto serial = k_mers_naive_serial( genome, guide, max_mismatch, k).hits;
    auto index_parallel = k_mers_naive_parallel(genome, guide, max_mismatch, k, 16).hits;
    auto full_parallel = k_mers_naive_parallel(genome, guide, max_mismatch, k, 16).hits;


    std::sort(index_parallel.begin(), index_parallel.end(),[](const Hit& a, const Hit& b) {
            return a.position < b.position;
     });


    std::sort(serial.begin(), serial.end(), [](const Hit& a, const Hit& b) {
                return a.position < b.position;
     });


    assert(serial.size() == 2);

    assert(serial[0].position == 0);
    assert(serial[0].mismatches == 0);

    assert(serial[1].position == 5);
    assert(serial[1].mismatches == 0);

    assert(index_parallel.size() == 2);

    assert(index_parallel[0].position == 0);
    assert(index_parallel[0].mismatches == 0);

    assert(index_parallel[1].position == 5);
    assert(index_parallel[1].mismatches == 0);

    std::cout << "Complete Match Test Passed!\n";
}

void noMatchTest() {
    std::string genome = "AAAAAAA";
    std::string guide = "TTTT";
    int max_mismatch = 0;
    int k = 2;

    auto serial = k_mers_naive_serial(genome, guide, max_mismatch, k).hits;
    auto parallel = k_mers_naive_parallel(genome, guide, max_mismatch, k , 16).hits;

    assert(serial.empty());
    assert(parallel.empty());

    std::cout << "No Match Test Passed!\n";
}

void oneMatchTest() {
    std::string genome = "GGATCGCC";
    std::string guide = "ATCG";
    int max_mismatch = 0;
    int k = 2;

    auto serial = k_mers_naive_serial(genome, guide, max_mismatch, k).hits;
    auto parallel = k_mers_naive_parallel(genome, guide, max_mismatch, k , 16).hits;

    assert(serial.size() == 1);
    assert(parallel.size() == 1);

    assert(serial[0].position == 2);
    assert(serial[0].mismatches == 0);

    assert(parallel[0].position == 2);
    assert(parallel[0].mismatches == 0);

    std::cout << "One Match Test Passed!\n";
}

void oneMismatchTest() {
    std::string genome = "GGATGGCC";
    std::string guide = "ATCG";
    int max_mismatch = 1;
    int k = 2;

    auto serial = k_mers_naive_serial(genome, guide, max_mismatch, k).hits;
    auto parallel = k_mers_naive_parallel(genome, guide, max_mismatch, k , 16).hits;

    assert(serial.size() == 1);
    assert(parallel.size() == 1);

    assert(serial[0].position == 2);
    assert(serial[0].mismatches == 1);

    assert(parallel[0].position == 2);
    assert(parallel[0].mismatches == 1);

    std::cout << "One Mismatch Test Passed!\n";
}

void invalidKTest() {
    std::string genome = "GGATGGCC";
    std::string guide = "ATCG";
    int max_mismatch = 1;
    int k = 10;

    auto serial = k_mers_naive_serial(genome, guide, max_mismatch, k).hits;
    auto parallel = k_mers_naive_parallel(genome, guide, max_mismatch, k , 16).hits;

    assert(serial.empty());
    assert(parallel.empty());;

    std::cout << "Invalid K Test Passed!\n";
}

void emptyGenomeTest() {
    std::string genome = "";
    std::string guide = "ATCG";
    int max_mismatch = 1;
    int k = 10;

    auto serial = k_mers_naive_serial(genome, guide, max_mismatch, k).hits;
    auto parallel = k_mers_naive_parallel(genome, guide, max_mismatch, k , 16).hits;

    assert(serial.empty());
    assert(parallel.empty());;

    std::cout << "Empty Genome Test Passed!\n";
}

void emptyGRNATest() {
    std::string genome = "GGATGGCC";
    std::string guide = "";
    int max_mismatch = 1;
    int k = 10;

    auto serial = k_mers_naive_serial(genome, guide, max_mismatch, k).hits;
    auto parallel = k_mers_naive_parallel(genome, guide, max_mismatch, k , 16).hits;

    assert(serial.empty());
    assert(parallel.empty());;

    std::cout << "Empty gRNA Test Passed!\n";
}
