import os 
import sys 
import csv

PROJECT_ROOT = os.path.abspath("../../..")
sys.path.append(PROJECT_ROOT)
sys.path.append(os.path.join(PROJECT_ROOT, "build"))

import offtarget_search as ot
from python.Lab.modules.kmer import (benchmark_naive_serial, benchmark_naive_parallel, benchmark_array_serial, benchmark_array_parallel, benchmark_hash_serial, benchmark_hash_parallel)


def threads_test_runner( GENOME_PATH , GUIDE , MAX_MISMATCH ,  K , THREADS_LIST , REPEATS , WARMUPS , clear_cache):

    genome = ot.load_fasta(GENOME_PATH)

    print("Module loaded") 
    print(f"Genome length: {len(genome)}") 
    print(f"Guide: {GUIDE}") 
    print(f"Maximum mismatches: {MAX_MISMATCH}") 
    print(f"K: {K}") 
    print(f"Repeats: {REPEATS}") 
    print(f"Warmups: {WARMUPS}") 
    print()

    results = {}


    print("Running Naive serial benchmark...")

    result = benchmark_naive_serial( ot, genome, GUIDE, MAX_MISMATCH, K, REPEATS, WARMUPS )
    result["speedup"] = 1.0
    results["naive_serial"] = result
    print( f"Naive Serial → " f"Mean: {result['mean']:.2f} ms | " f"Median: {result['median']:.2f} ms | " f"Candidates: {result['candidates']} | " f"Hits: {result['hits']}" )
    print()


    for threads in THREADS_LIST[1:]:
        print( f"Running Naive parallel " f"with {threads} threads..." ) 
        result = benchmark_naive_parallel( ot, genome, GUIDE, MAX_MISMATCH, K, threads, REPEATS, WARMUPS ) 
        speedup = ( results["naive_serial"]["median"] / result["median"] ) 
        result["speedup"] = speedup 
        results[f"naive_{threads}"] = result 
        print( f"Naive {threads:2d} threads → " f"Mean: {result['mean']:.2f} ms | " f"Median: {result['median']:.2f} ms | " f"Speedup: {speedup:.2f}x | " f"Candidates: {result['candidates']} | " f"Hits: {result['hits']}" ) 
    print()


    print("Running Array serial benchmark...")

    result = benchmark_array_serial( ot, genome, GUIDE, MAX_MISMATCH, REPEATS, WARMUPS, clear_cache )
    result["speedup"] = 1.0
    results["array_serial"] = result
    print( f"Array Serial → " f"Mean: {result['mean']:.4f} ms | " f"Median: {result['median']:.4f} ms | " f"Candidates: {result['candidates']} | " f"Hits: {result['hits']}" )
    print()
    

    for threads in THREADS_LIST[1:]:
        print( f"Running Array parallel " f"with {threads} threads..." ) 
        result = benchmark_array_parallel( ot, genome, GUIDE, MAX_MISMATCH, threads, REPEATS, WARMUPS, clear_cache ) 
        speedup = ( results["array_serial"]["median"] / result["median"] ) 
        result["speedup"] = speedup  
        results[f"array_{threads}"] = result 
        print( f"Array {threads:2d} threads → " f"Mean: {result['mean']:.4f} ms | " f"Median: {result['median']:.4f} ms | " f"Speedup: {speedup:.4f}x | " f"Candidates: {result['candidates']} | " f"Hits: {result['hits']}" ) 
    print()
    

    print("Running Hash serial benchmark...")

    result = benchmark_hash_serial( ot, genome, GUIDE, MAX_MISMATCH, REPEATS, WARMUPS , clear_cache)
    result["speedup"] = 1.0
    results["hash_serial"] = result
    print( f"Hash Serial → " f"Mean: {result['mean']:.4f} ms | " f"Median: {result['median']:.4f} ms | " f"Candidates: {result['candidates']} | " f"Hits: {result['hits']}" )
    print()
    

    for threads in THREADS_LIST[1:]:
        print( f"Running Hash parallel " f"with {threads} threads..." ) 
        result = benchmark_hash_parallel( ot, genome, GUIDE, MAX_MISMATCH, threads, REPEATS, WARMUPS, clear_cache ) 
        speedup = ( results["hash_serial"]["median"] / result["median"] ) 
        result["speedup"] = speedup 
        results[f"hash_{threads}"] = result 
        print( f"Hash {threads:2d} threads → " f"Mean: {result['mean']:.4f} ms | " f"Median: {result['median']:.4f} ms | " f"Speedup: {speedup:.4f}x | " f"Candidates: {result['candidates']} | " f"Hits: {result['hits']}" ) 
    print()

    return results



def save_threads_csv( results , CSV_PATH ):

    os.makedirs( os.path.dirname(CSV_PATH), exist_ok=True )
    with open(CSV_PATH, "w", newline="") as file:

        writer = csv.writer(file) 
        writer.writerow([ "algorithm", "threads", 
                         "index_mean_ms", "search_mean_ms", "collect_mean_ms", "total_mean_ms", 
                         "index_median_ms", "search_median_ms", "collect_median_ms", "total_median_ms",
                         "speedup", "candidates", "hits"])

        writer.writerow(["naive", 1, 
                         results["naive_serial"]["index_mean"], results["naive_serial"]["search_mean"] ,None, results["naive_serial"]["mean"] ,
                         results["naive_serial"]["index_median"], results["naive_serial"]["search_median"],None, results["naive_serial"]["median"] ,
                         1.0,results["naive_serial"]["candidates"],results["naive_serial"]["hits"]])
        for threads in THREADS_LIST[1:]: 
            result = results[f"naive_{threads}"] 
            writer.writerow(["naive", threads, 
                             result["index_mean"], result["search_mean"], None, result["mean"], 
                             result["index_median"], result["search_median"], None, result["median"],  
                             result["speedup"], result["candidates"],result["hits"]])

        writer.writerow([ "array", 1, 
                         results["array_serial"]["index_mean"], results["array_serial"]["search_mean"],results["array_serial"]["collect_mean"], results["array_serial"]["mean"],
                         results["array_serial"]["index_median"], results["array_serial"]["search_median"],results["array_serial"]["collect_median"],  results["array_serial"]["median"],
                         1.0, results["array_serial"]["candidates"], results["array_serial"]["hits"]])
        for threads in THREADS_LIST[1:]: 
            result = results[f"array_{threads}"] 
            writer.writerow([ "array", threads, 
                            result["index_mean"], result["search_mean"],  result["collect_mean"], result["mean"],
                            result["index_median"], result["search_median"], result["collect_median"], result["median"],
                            result["speedup"], result["candidates"], result["hits"]])
                  
        writer.writerow(["hash", 1, results["hash_serial"]["index_mean"], results["hash_serial"]["search_mean"], results["hash_serial"]["collect_mean"],results["hash_serial"]["mean"],
                          results["hash_serial"]["index_median"], results["hash_serial"]["search_median"], results["hash_serial"]["collect_median"], results["hash_serial"]["median"], 
                          1.0, results["hash_serial"]["candidates"], results["hash_serial"]["hits"]])
        for threads in THREADS_LIST[1:]: 
            result = results[f"hash_{threads}"] 
            writer.writerow(["hash", threads, 
                             result["index_mean"], result["search_mean"], result["collect_mean"], result["mean"],
                             result["index_median"], result["search_median"], result["collect_median"], result["median"],
                             result["speedup"], result["candidates"], result["hits"]])
     
    print(f"Results saved to: {CSV_PATH}")




THREADS_LIST = [1, 2, 4, 8, 16]
REPEATS = 10 
WARMUPS = 2
GENOME= os.path.join( PROJECT_ROOT, "data", "GCF_000005845.2_ASM584v2_genomic.fna" )


cache_test_guide= "ACGCGCCGATTGTTGCGAGA" 
cache_test_maxmismatch= 2
cache_test_k = 6
cache_test_cvs_path= os.path.join( PROJECT_ROOT, "python/Lab/results", "threads_results.csv" )
cache_test_clear_cache = True # کش رو حذف کنم ؟ بله 
cache_test_result = threads_test_runner( GENOME , cache_test_guide , cache_test_maxmismatch , cache_test_k  , THREADS_LIST , REPEATS , WARMUPS  , cache_test_clear_cache)
save_threads_csv(cache_test_result, cache_test_cvs_path )

