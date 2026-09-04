#include <iostream>
#include <chrono>
#include <omp.h>
#include "GenomeLoader.h"
#include "Sliding_window.h"

void threadsTest(std::string genome, std::string guide, int threads_number);

int main()
{
    std::string genome = loadFasta("data/GCF_000005845.2_ASM584v2_genomic.fna");
    std::string gRNA = "ATTAGCGGCCAGGATG";
    std::cout << "Genome size: "<< genome.size()<< '\n';

    threadsTest(genome, gRNA, 1);
    threadsTest(genome, gRNA, 2);
    threadsTest(genome, gRNA, 4);
    threadsTest(genome, gRNA, 8);
    threadsTest(genome, gRNA, 16);

}


void threadsTest( std::string genome , std::string guide , int threads_number) {

    std::cout << "\n\n Threads number: " << threads_number << std::endl;


    auto serial_time = 0;
    auto parallel_time = 0;
    std::cout << "Serial - ";
    std::cout << "Hits: \n";

    for (int i=0 ; i<5 ; i++){
    auto start = std::chrono::high_resolution_clock::now();
    auto hits = slidingWindow_serial(genome, guide , 2).hits;
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    serial_time += duration.count();

    std::cout  << duration.count() << " ms -  ";
    std::cout << hits.size() << " , ";
    }

    std::cout << "\n\nParallel - ";
    std::cout << "Hits: \n";
    for (int i = 0; i < 5; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        auto hits = slidingWindow_parallel(genome, guide, 2, threads_number).hits;
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        parallel_time += duration.count();
        std::cout << duration.count() << " ms - ";
        std::cout  << hits.size() << " , ";
    }

    std::cout << "\nSpeedup: " << (double)serial_time / parallel_time << '\n';
}
