import os 
import sys 
import math
import csv

PROJECT_ROOT = os.path.abspath("../../..")
sys.path.append(PROJECT_ROOT)
sys.path.append(os.path.join(PROJECT_ROOT, "build"))

import offtarget_search as ot
from python.Lab.modules.kmer import (benchmark_naive_serial, benchmark_naive_parallel, benchmark_array_serial, benchmark_array_parallel, benchmark_hash_serial, benchmark_hash_parallel)


def kmer_cache_test( GENOME_PATH , GUIDE_LIST , MAX_MISMATCH , THREADS ):

    genome = ot.load_fasta(GENOME_PATH)

    print("Module loaded") 
    print(f"Genome length: {len(genome)}") 
    print(f"Guide: {GUIDE_LIST}") 
    print(f"Maximum mismatches: {MAX_MISMATCH}") 
    print(f"Threads: {THREADS}") 
    print()

    results = []
  

    for g in GUIDE_LIST:
        k = math.floor(len(g)  / float(MAX_MISMATCH+ 1));
        print( f"Running Array serial " f"with {k}-k..." ) 
        result =ot.kmer_serial_array(genome, g,  MAX_MISMATCH )
        print( f"Array Serial {k:2d}-k → " f"Index: { result.index_time :.2f} ms | " f"Serach: {result.search_time :.2f} ms | " f"Collect: {result.collect_time :.2f} ms | " f"Candidates: {result.candidate_count} | "  )
        results.append({
            "algorithm": "array",
            "mode": "serial",
            "guide": g,
            "guide_length": len(g),
            "k": k,
            "threads": THREADS,
            "index_time_ms": result.index_time,
            "collect_time_ms": result.collect_time,
            "search_time_ms": result.search_time,
            "total_time_ms": result.index_time + result.search_time + result.collect_time,
            "candidates": result.candidate_count,
            "hits": len(result.hits)
        })
    print()

      
    for g in GUIDE_LIST:
        k = math.floor(len(g)  / float(MAX_MISMATCH+ 1));
        print( f"Running Array parallel " f"with {k}-k..." ) 
        result =ot.kmer_parallel_array(genome, g,  MAX_MISMATCH , THREADS)
        print( f"Arrary parallel {k:2d}-k → " f"Index: { result.index_time :.2f} ms | " f"Serach: {result.search_time :.2f} ms | " f"Collect: {result.collect_time :.2f} ms | "  f"Candidates: {result.candidate_count} | " )
        results.append({
            "algorithm": "array",
            "mode": "parallel",
            "guide": g,
            "guide_length": len(g),
            "k": k,
            "threads": THREADS,
            "index_time_ms": result.index_time,
            "collect_time_ms": result.collect_time,
            "search_time_ms": result.search_time,
            "total_time_ms": result.index_time + result.search_time + result.collect_time,
            "candidates": result.candidate_count,
            "hits": len(result.hits)
        })
    print()

    

    for g in GUIDE_LIST:
        k = math.floor(len(g)  / float(MAX_MISMATCH+ 1));
        print( f"Running Hash serial " f"with {k}-k..." ) 
        result =ot.kmer_serial_hash(genome, g,  MAX_MISMATCH )
        print( f"Hash Serial {k:2d}-k → " f"Index: { result.index_time :.2f} ms | " f"Serach: {result.search_time :.2f} ms | " f"Collect: {result.collect_time :.2f} ms |  " f"Candidates: {result.candidate_count} | " )
        results.append({
            "algorithm": "hash",
            "mode": "serial",
            "guide": g,
            "guide_length": len(g),
            "k": k,
            "threads": THREADS,
            "index_time_ms": result.index_time,
            "collect_time_ms": result.collect_time,
            "search_time_ms": result.search_time,
            "total_time_ms": result.index_time + result.search_time +result.collect_time,
            "candidates": result.candidate_count,
            "hits": len(result.hits)
        })
    print()
       
    for g in GUIDE_LIST:
        k = math.floor(len(g)  / float(MAX_MISMATCH+ 1));
        print( f"Running Hash  parallel " f"with {k}-k..." ) 
        result =ot.kmer_parallel_hash(genome, g,  MAX_MISMATCH , THREADS)
        print( f"Hash parallel {k:2d}-k → " f"Index: { result.index_time :.2f} ms | " f"Serach: {result.search_time :.2f} ms | " f"Collect: {result.collect_time :.2f} ms | " f"Candidates: {result.candidate_count} | "  )
        results.append({
            "algorithm": "hash",
            "mode": "parallel",
            "guide": g,
            "guide_length": len(g),
            "k": k,
            "threads": THREADS,
            "index_time_ms": result.index_time,
            "collect_time_ms": result.collect_time,
            "search_time_ms": result.search_time,
            "total_time_ms": result.index_time + result.search_time + result.collect_time,
            "candidates": result.candidate_count,
            "hits": len(result.hits)
        })
    print()

    return results


def save_kmer_cache_csv(results, csv_path):
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)

    with open(csv_path, "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "algorithm",
            "mode",
            "guide",
            "guide_length",
            "k",
            "threads",
            "index_time_ms",
            "collect_time_ms",
            "search_time_ms",
            "total_time_ms",
            "candidates",
            "hits"
        ])

        for item in results:
            writer.writerow([
                item["algorithm"],
                item["mode"],
                item["guide"],
                item["guide_length"],
                item["k"],
                item["threads"],
                item["index_time_ms"],
                item["collect_time_ms"],
                item["search_time_ms"],
                item["total_time_ms"],
                item["candidates"],
                item["hits"]
            ])

    print(f"Results saved to: {csv_path}")



GENOME= os.path.join( PROJECT_ROOT, "data", "GCF_000005845.2_ASM584v2_genomic.fna" )
MAX_MISMATCH = 2

csv_path= os.path.join( PROJECT_ROOT, "python/Lab/results", "cache_results.csv" )

guide_list= [ "ACCATTACCCCCATCGCC", "AGGGATGGTCGGCATGGC","AATTGATCAACACAGCGAC","ATGAGATGCTGA", "AGTTGTACATGG" ]
results = kmer_cache_test( GENOME , guide_list , MAX_MISMATCH , 2  )
save_kmer_cache_csv(results, csv_path)
