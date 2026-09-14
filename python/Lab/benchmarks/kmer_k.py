import os 
import sys 
import csv

PROJECT_ROOT = os.path.abspath("../../..")
sys.path.append(PROJECT_ROOT)
sys.path.append(os.path.join(PROJECT_ROOT, "build"))

import offtarget_search as ot
from python.Lab.modules.kmer import (benchmark_naive_serial, benchmark_naive_parallel, benchmark_array_serial, benchmark_array_parallel, benchmark_hash_serial, benchmark_hash_parallel)


def k_test_runner( GENOME_PATH , GUIDE_LIST , MAX_MISMATCH ,  K_LIST , THREADS , REPEATS , WARMUPS , clear_cache):

    genome = ot.load_fasta(GENOME_PATH)

    print("Module loaded") 
    print(f"Genome length: {len(genome)}") 
    print(f"Guide: {GUIDE_LIST}") 
    print(f"Maximum mismatches: {MAX_MISMATCH}") 
    print(f"Threads: {THREADS}") 
    print(f"Repeats: {REPEATS}") 
    print(f"Warmups: {WARMUPS}") 
    print()

    results = {}


    for i in range(5):
        k = K_LIST[i]
        guide = GUIDE_LIST[i]
        print( f"Running Naive serial " f"with {k}-k..." ) 
        serial_result = benchmark_naive_serial( ot, genome, guide , MAX_MISMATCH, k , REPEATS, WARMUPS )
        serial_result["speedup"] = 1.0
        results[f"naive_serial_{k}"] = serial_result 
        print( f"Naive Serial {k:2d}-k → " f"Mean: {serial_result['index_mean']:.2f} ms | " f"Median: {serial_result['index_median']:.2f} ms | " f"Candidates: {serial_result['candidates']} | " f"Hits: {serial_result['hits']}" )
        print()
        
    for i in range(5):
        k = K_LIST[i]
        guide = GUIDE_LIST[i]
        print( f"Running Naive parallel " f"with {k}-k..." ) 
        result = benchmark_naive_parallel( ot, genome, guide, MAX_MISMATCH, k , THREADS, REPEATS, WARMUPS ) 
        speedup = ( results[f"naive_serial_{k}"]["median"] / result["median"] ) 
        result["speedup"] = speedup 
        results[f"naive_parallel_{k}"] = result 
        print( f"Naive {k:2d}-k → " f"Mean: {result['index_mean']:.2f} ms | " f"Median: {result['index_median']:.2f} ms | " f"Speedup: {speedup:.2f}x | " f"Candidates: {result['candidates']} | " f"Hits: {result['hits']}" ) 
        print()

  

    for i in range(5):
        k = K_LIST[i]
        guide = GUIDE_LIST[i]
        print( f"Running Array serial " f"with {k}-k..." ) 
        serial_result = benchmark_array_serial( ot, genome, guide, MAX_MISMATCH, REPEATS, WARMUPS, clear_cache )
        serial_result["speedup"] = 1.0
        results[f"array_serial_{k}"] = serial_result 
        print( f"Array Serial {k:2d}-k → " f"Mean: {serial_result['index_mean']:.2f} ms | " f"Median: {serial_result['index_median']:.2f} ms | " f"Candidates: {serial_result['candidates']} | " f"Hits: {serial_result['hits']}" )
        print()
        
    for i in range(5):
        k = K_LIST[i]
        guide = GUIDE_LIST[i]
        print( f"Running Array parallel " f"with {k}-k..." ) 
        result =benchmark_array_parallel( ot, genome, guide, MAX_MISMATCH, THREADS, REPEATS, WARMUPS, clear_cache ) 
        speedup = ( results[f"array_serial_{k}"]["median"] / result["median"] ) 
        result["speedup"] = speedup 
        results[f"array_parallel_{k}"] = result 
        print( f"Array  {k:2d}-k → " f"Mean: {result['index_mean']:.2f} ms | " f"Median: {result['index_median']:.2f} ms | " f"Candidates: {result['candidates']} | " f"Hits: {result['hits']}" )
        print()

    

    for i in range(5):
        k = K_LIST[i]
        guide = GUIDE_LIST[i]
        print( f"Running Hash serial " f"with {k}-k..." ) 
        serial_result = benchmark_hash_serial( ot, genome, guide, MAX_MISMATCH, REPEATS, WARMUPS, clear_cache )
        serial_result["speedup"] = 1.0
        results[f"hash_serial_{k}"] = serial_result 
        print( f"Hash Serial {k:2d}-k → " f"Mean: {serial_result['index_mean']:.2f} ms | " f"Median: {serial_result['index_median']:.2f} ms | " f"Candidates: {serial_result['candidates']} | " f"Hits: {serial_result['hits']}" )
        print()
        
    for i in range(5):
        k = K_LIST[i]
        guide = GUIDE_LIST[i]
        print( f"Running Hash  parallel " f"with {k}-k..." ) 
        result =benchmark_hash_parallel( ot, genome, guide, MAX_MISMATCH, THREADS, REPEATS, WARMUPS, clear_cache ) 
        speedup = ( results[f"hash_serial_{k}"]["median"] / result["median"] ) 
        result["speedup"] = speedup 
        results[f"hash_parallel_{k}"] = result 
        print( f"Hash {k:2d}-k → " f"Mean: {result['index_mean']:.2f} ms | " f"Median: {result['index_median']:.2f} ms | " f"Speedup: {speedup:.2f}x | " f"Candidates: {result['candidates']} | " f"Hits: {result['hits']}" ) 
        print()

    return results



def save_k_results_csv(results, CSV_PATH):

    os.makedirs(os.path.dirname(CSV_PATH), exist_ok=True)

    with open(CSV_PATH, "w", newline="") as file:

        writer = csv.writer(file)

        writer.writerow(["algorithm", "mode", "k", "threads",
            "index_mean_ms", "search_mean_ms", "collect_mean_ms", "total_mean_ms",
            "index_median_ms", "search_median_ms", "collect_median_ms", "total_median_ms",
            "speedup", "candidates", "hits"])

        for key, result in results.items():

            parts = key.split("_")

            algorithm = parts[0]

            if algorithm not in ["naive", "array", "hash"]:
                continue

            mode = parts[1]
            k = int(parts[2])

            if mode == "serial":
                threads = 1

            elif mode == "parallel":
                threads = None

            else:
                continue

            writer.writerow([ algorithm, mode, k, threads,
                             result["index_mean"], result["search_mean"], result.get("collect_mean"), result["mean"],
                             result["index_median"], result["search_median"], result.get("collect_median"), result["median"],
                             result["speedup"], result["candidates"], result["hits"]])

    print(f"Results saved to: {CSV_PATH}")




THREADS_LIST = [1, 2, 4, 8, 16]
REPEATS = 10 
WARMUPS = 2
GENOME= os.path.join( PROJECT_ROOT, "data", "GCF_000005845.2_ASM584v2_genomic.fna" )



guide_list= ["ACCATTA" ,
             "ACCATTACCCCC",
             "ACCATTACCCCCATCGCC",
             "ACCATTACCCCCATCGCCCAGTTC", 
             "ACCATTACCCCCATCGCCCAGTTCCAGATCC"]
k_list=[2,4,6,8,10]
k_test_maxmismatch=2
k_test_cvs_path= os.path.join( PROJECT_ROOT, "python/Lab/results" , "k_results.csv")
k_test_cache = True #کش ‍اک میشه
k_test_result = k_test_runner( GENOME , guide_list , k_test_maxmismatch , k_list , 2 , REPEATS , WARMUPS , k_test_cache)
save_k_results_csv(k_test_result,k_test_cvs_path)
