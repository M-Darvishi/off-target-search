#include <iostream>
#include <cassert>
#include <algorithm>
#include "Sliding_window.h"

void simpleTest();
void compeletMatchTest();
void noMatchTest();
void oneMismatchTest();
void emptyGenomeTest();
void emptyGRNATest();
void gRNA_longerThanGenome();

int main()
{
    simpleTest();
    compeletMatchTest();
    noMatchTest();
    oneMismatchTest();
    emptyGenomeTest();
    emptyGRNATest();
    gRNA_longerThanGenome();

    return 0;
}

void simpleTest() {
    std::string genome = "ATGCGTACGTAGCTAGCTAGCATGCAATGCGTAGCTAGCATGCA";
    std::string gRNA = "ATGCGTACGTAGCTAGCTA";

    int max_mismatch = 2;

    auto serial =
        slidingWindow_serial(genome, gRNA, max_mismatch).hits;

    auto parallel =
        slidingWindow_parallel(genome, gRNA, max_mismatch, 5).hits;

    assert(serial.size() == parallel.size());

    std::cout << " Simple Test passed!" << std::endl;
}

void compeletMatchTest() {
    std::string genome = "ATCGCATCG";
    std::string gRNA = "ATCG";
    int max_mismatch = 2;

    auto serial =
        slidingWindow_serial(genome, gRNA, max_mismatch).hits;

    auto parallel =
        slidingWindow_parallel(genome, gRNA, max_mismatch, 5).hits;

    std::sort(parallel.begin(), parallel.end(), [](const Hit& a, const Hit& b) {
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

    assert(parallel.size() == 2);
    assert(parallel[0].position == 0);
    assert(parallel[0].mismatches == 0);

    assert(parallel[1].position == 5);
    assert(parallel[1].mismatches == 0);

    std::cout << " Compelet Match Test passed!" << std::endl;
}

void noMatchTest() {
    std::string genome = "AAAAAA";
    std::string gRNA = "TTTT";

    int max_mismatch = 2;

    auto serial =
        slidingWindow_serial(genome, gRNA, max_mismatch).hits;

    auto parallel =
        slidingWindow_parallel(genome, gRNA, max_mismatch, 5).hits;

    assert(serial.size() == 0);
    assert(parallel.size() == 0);

    std::cout << " No Match Test passed!" << std::endl;
}

void oneMismatchTest(){
    std::string genome = "CTATGGCATGCT";
    std::string gRNA = "ATCG";
    int max_mismatch = 1;

    auto serial =
        slidingWindow_serial(genome, gRNA, max_mismatch).hits;

    auto parallel =
        slidingWindow_parallel(genome, gRNA, max_mismatch, 5).hits;

    assert(serial.size() == 1);
    assert(serial[0].position == 2);
    assert(serial[0].mismatches == 1);

    assert(parallel.size() == 1);
    assert(parallel[0].position == 2);
    assert(parallel[0].mismatches == 1);

    std::cout << " One Mismatch Test passed!" << std::endl;
}

void emptyGenomeTest() {
    std::string genome = "";
    std::string gRNA = "ATGC";
    int max_mismatch = 2;

    auto serial =
        slidingWindow_serial(genome, gRNA, max_mismatch).hits;

    auto parallel =
        slidingWindow_parallel(genome, gRNA, max_mismatch, 5).hits;

    assert(serial.size() == 0);
    assert(parallel.size() == 0);

    std::cout << " Empty Genome Test passed!" << std::endl;
}

void emptyGRNATest() {
    std::string genome = "ATCGATCG";
    std::string gRNA = "";
    int max_mismatch = 2;

    auto serial =
        slidingWindow_serial(genome, gRNA, max_mismatch).hits;

    auto parallel =
        slidingWindow_parallel(genome, gRNA, max_mismatch, 5).hits;

    assert(serial.size() == 0);
    assert(parallel.size() == 0);

    std::cout << " Empty gRNA Test passed!" << std::endl;
}

void gRNA_longerThanGenome() {
    std::string genome = "ATCG";
    std::string gRNA = "ATCGTCG";
    int max_mismatch = 2;

    auto serial =
        slidingWindow_serial(genome, gRNA, max_mismatch).hits;

    auto parallel =
        slidingWindow_parallel(genome, gRNA, max_mismatch, 5).hits;

    assert(serial.size() == 0);
    assert(parallel.size() == 0);

    std::cout << " gRNA Longer Than Genome Test passed!" << std::endl;

}
