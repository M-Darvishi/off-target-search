import os 
import sys 

PROJECT_ROOT = os.path.abspath("../../..")
sys.path.append(PROJECT_ROOT)
sys.path.append(os.path.join(PROJECT_ROOT, "build"))

import offtarget_search as ot
from python.Lab.modules.kmer import (benchmark_hash_serial)


def kmer_hash_memory( GENOME_PATH , GUIDE , MAX_MISMATCH ,  K, REPEATS , WARMUPS , clear_cache):

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

    print("Running Hash serial benchmark...")

    result = benchmark_hash_serial( ot, genome, GUIDE, MAX_MISMATCH, K, REPEATS, WARMUPS )
    results["hash_serial"] = result
    print( f"Hash Serial → " f"Time: {result['index_mean']:.2f} ms | " f"Candidates: {result['candidates']} | " f"Hits: {result['hits']}" )
    print()



THREADS_LIST = [1, 2, 4, 8, 16]
REPEATS = 10 
WARMUPS = 2
GENOME = os.path.join( PROJECT_ROOT, "data", "GCF_000005845.2_ASM584v2_genomic.fna" )
GUIDE_LIST = ["ACCATTA" ,
             "ACCATTACCCCC",
             "ACCATTACCCCCATCGCC",
             "ACCATTACCCCCATCGCCCAGTTC", 
             "ACCATTACCCCCATCGCCCAGTTCCAGATCC"]
K_LIST = [2,4,6,8,10]
MAX_MISSMATCH =2
CSV_PATH = os.path.join( PROJECT_ROOT, "CSV", "kmer_k_results.csv")
CLEAR_CACHE = False # با استفادخ از کش

while(True):
    k = int(input("Enter '0' to exite. \nChose betwean [2, 4, 6, 8, 10]: "))
    if k in K_LIST:
        i = K_LIST.index(k)
        result = kmer_hash_memory( GENOME , GUIDE_LIST[i] , MAX_MISSMATCH , k  , REPEATS , WARMUPS , CLEAR_CACHE)
        break
    elif(k==0):
        break
    else:
        continue
