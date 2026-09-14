import os
import sys
import csv

PROJECT_ROOT = os.path.abspath("../../..")
sys.path.append(PROJECT_ROOT)
sys.path.append(os.path.join(PROJECT_ROOT, "build"))

import offtarget_search as ot

from python.Lab.modules.sliding_window import ( benchmark_serial, benchmark_parallel)

GENOME_PATH = os.path.join( PROJECT_ROOT, "data", "GCF_000005845.2_ASM584v2_genomic.fna")
GUIDE = "AGCTTTTCATTCTGACTGCA"
MAX_MISMATCH = 2
THREADS_LIST = [1, 2, 4, 8, 16]

REPEATS = 10
WARMUPS = 2

CSV_PATH = os.path.join( PROJECT_ROOT, "python/Lab/results", "sliding_window_results.csv")
 
genome = ot.load_fasta(GENOME_PATH)

print("Module loaded") 
print(f"Genome length: {len(genome)}") 
print(f"Guide: {GUIDE}") 
print(f"Maximum mismatches: {MAX_MISMATCH}") 
print(f"Repeats: {REPEATS}") 
print()

 
print("Running serial benchmark...")
serial_result = benchmark_serial( ot, genome, GUIDE, MAX_MISMATCH, REPEATS, WARMUPS )
print( f"Serial → " f"Mean: {serial_result['mean']:.2f} ms | " f"Median: {serial_result['median']:.2f} ms | " f"Hits: {serial_result['hits']}" )
print()

results = { "serial": { **serial_result, "speedup": 1.0 } }

for threads in THREADS_LIST[1:]:
    print(f"Running {threads} threads...") 
    result = benchmark_parallel( ot, genome, GUIDE, MAX_MISMATCH, threads, REPEATS, WARMUPS ) 
    speedup = serial_result["median"] / result["median"] 
    result["speedup"] = speedup 
    results[threads] = result 
    print( f"Threads: {threads:2d} → " f"Mean: {result['mean']:.2f} ms | " f"Median: {result['median']:.2f} ms | " f"Speedup: {speedup:.2f}x | " f"Hits: {result['hits']}" ) 
print()


#Saving CSV

with open(CSV_PATH, "w", newline="") as file:

    writer = csv.writer(file)
    writer.writerow([
        "threads",
        "mean_ms",
        "median_ms",
        "speedup",
        "hits"
    ])

    writer.writerow([
        1,
        serial_result["mean"],
        serial_result["median"],
        1.0,
        serial_result["hits"]
    ])

    for threads in THREADS_LIST[1:]:

        result = results[threads]

        writer.writerow([
            threads,
            result["mean"],
            result["median"],
            result["speedup"],
            result["hits"]
        ])

print(f"Results saved to: {CSV_PATH}")
