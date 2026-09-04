#include <iostream>
#include <chrono>
#include <omp.h>
#include <cmath>
#include "GenomeLoader.h"
#include "K_mers_core.h"
#include "K_mers_naive.h"
#include "K_mers_index_builder_parallel.h"

void threadsTest(std::string genome, std::string guide,int  max_mismatch);

void naive_benchmark(std::string genome, std::string guide, int  max_mismatch, int threads_number);
void hash_benchmark(std::string genome, std::string guide, int  max_mismatch, int threads_number);
void array_benchmark(std::string genome, std::string guide, int  max_mismatch, int threads_number);

int main()
{
    std::string genome = loadFasta("data/GCF_000005845.2_ASM584v2_genomic.fna");
    std::string gRNA = "AGCTTTTCATTCTGACTGCA";
    //std::string gRNA = "AGCTTTTCATTCTGACTGCAACGGGCAATA";
    int max_mismatch = 2;
    std::cout << "Genome size: " << genome.size() << '\n';

   threadsTest(genome, gRNA, max_mismatch);

}

void threadsTest(std::string genome, std::string guide, int max_mismatch) 
{ 
    int threads[5] = {1,2,4,8,16};

    std::cout << "\n\n\n\t===Naive== " << std::endl;
    for (int i : threads)
        naive_benchmark(genome, guide, max_mismatch, i);
    
    std::cout << "\n\n\n\t===Hash=== " << std::endl;
    for (int i : threads)
        hash_benchmark(genome, guide, max_mismatch, i);

    std::cout << "\n\n\n\t===Array== " << std::endl;
    for (int i : threads)  
        array_benchmark(genome, guide, max_mismatch, i);  
}


void naive_benchmark(std::string genome, std::string guide, int  max_mismatch, int threads_number) 
{
    int k = floor(guide.length() / float(2 + 1));

    auto naive_serial_time = 0;
    auto naive_parallel_time = 0;

    std::cout << "\nNaive Serial - ";
    std::cout << "Hits: \n\t";
    for (int i = 0; i < 1; i++) {
        auto naive_start = std::chrono::high_resolution_clock::now();
        auto naive_hits = k_mers_naive_serial(genome, guide, max_mismatch, k).hits;
        auto naive_end = std::chrono::high_resolution_clock::now();

        auto naive_duration = std::chrono::duration_cast<std::chrono::milliseconds>(naive_end - naive_start);
        naive_serial_time += naive_duration.count();

        std::cout << naive_duration.count() << " ms -  ";
        std::cout << naive_hits.size() << " , ";
    }
    clear_kmer_cache();

    std::cout << "\nNaive Parallel - ";
    std::cout << "Hits: \n\t";
    for (int i = 0; i < 1; i++) {
        auto naive_start = std::chrono::high_resolution_clock::now();
        auto naive_hits = k_mers_naive_parallel(genome, guide, max_mismatch, k, threads_number).hits;
        auto naive_end = std::chrono::high_resolution_clock::now();

        auto naive_duration = std::chrono::duration_cast<std::chrono::milliseconds>(naive_end - naive_start);
        naive_parallel_time += naive_duration.count();

        std::cout << naive_duration.count() << " ms -  ";
        std::cout << naive_hits.size() << " , ";
    }
    clear_kmer_cache();

    std::cout << "\n\n~Threads number: " << threads_number;
    std::cout << " ~Speedup: " << (double)naive_serial_time / naive_parallel_time<<"\n\n";
}


void hash_benchmark(std::string genome, std::string guide, int  max_mismatch, int threads_number)
{
    auto hash_serial_time = 0;
    auto hash_parallel_time = 0;

    std::cout << "\nHash serial - ";
    std::cout << "Hits: \n\t";
    for (int i = 0; i < 5; i++) {
        auto kmer_hits = k_mers_serial_hash(genome, guide, max_mismatch);
        auto duration = kmer_hits.index_time + kmer_hits.collect_time + kmer_hits.search_time;;
        hash_serial_time += duration;

        std::cout << duration << " ms -  ";
        std::cout << kmer_hits.hits.size() << " , ";
    }
    clear_kmer_cache();

    std::cout << "\nHash Parallel - ";
    std::cout << "Hits: \n\t";
    for (int i = 0; i < 5; i++) {
        auto kmer_hits = k_mers_parallel_hash(genome, guide, max_mismatch, threads_number);
        auto duration = kmer_hits.index_time + kmer_hits.collect_time + kmer_hits.search_time;;
        hash_parallel_time += duration;

        std::cout << duration << " ms -  ";
        std::cout << kmer_hits.hits.size() << " , ";
    }
    clear_kmer_cache();

    std::cout << "\n\n~Threads number: " << threads_number;
    std::cout << " ~Speedup: " << (double)hash_serial_time / hash_parallel_time <<"\n\n";
}

void array_benchmark(std::string genome, std::string guide, int  max_mismatch,  int threads_number) 
{
    auto array_serial_time = 0;
    auto array_parallel_time = 0;

    std::cout << "\nArray Serial - ";
    std::cout << "Hits: \n\t";
    for (int i = 0; i < 5; i++) {
        auto kmer_hits = k_mers_serial_array(genome, guide, max_mismatch);
        auto duration = kmer_hits.index_time + kmer_hits.collect_time + kmer_hits.search_time;;
        array_serial_time += duration;

        std::cout << duration << " ms -  ";
        std::cout << kmer_hits.hits.size() << " , ";
    }
    clear_kmer_cache();

    std::cout << "\nArray Parallel - ";
    std::cout << "Hits: \n\t";
    for (int i = 0; i < 5; i++) {
        auto kmer_hits = k_mers_parallel_array(genome, guide, max_mismatch, threads_number);
        auto duration = kmer_hits.index_time + kmer_hits.collect_time + kmer_hits.search_time;;
        array_parallel_time += duration;

        std::cout << duration << " ms -  ";
        std::cout << kmer_hits.hits.size() << " , ";
    }
    clear_kmer_cache();

    std::cout << "\n\n~Threads number: " << threads_number;
    std::cout << " ~Speedup: " << (double)array_serial_time / array_parallel_time << "\n\n";
}

