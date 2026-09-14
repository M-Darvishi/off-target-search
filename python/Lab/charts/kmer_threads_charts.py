import csv
import matplotlib.pyplot as plt
import os
import sys

PROJECT_ROOT = os.path.abspath("../../..")
sys.path.append(PROJECT_ROOT)


def load_kmer_cache_results(csv_path):

    results = {
        "naive": {},
        "array": {},
        "hash": {}
    }

    with open(csv_path, "r", newline="") as file:

        reader = csv.DictReader(file)

        for row in reader:

            algorithm = row["algorithm"]
            threads = int(row["threads"])
            median = float(row["total_median_ms"])
            speedup = float(row["speedup"])

            if algorithm in results:
                results[algorithm][threads] = { "median": median, "speedup": speedup}
    return results



def plot_execution_time(results, threads):

    plt.figure(figsize=(8, 5))

    algorithms = ["naive", "array", "hash"]

    for algorithm in algorithms:

        serial_median = results[algorithm][1]["median"]

        parallel_medians = [  results[algorithm][t]["median"]  for t in threads ]

        line = plt.plot( threads, parallel_medians, marker="o", linewidth=2, label=f"{algorithm.capitalize()} Parallel")[0]

        plt.axhline( y=serial_median, linestyle="--", linewidth=1.5, color=line.get_color(), label=f"{algorithm.capitalize()} Serial")

    plt.xlabel("Number of Threads")
    plt.ylabel("Total Median Execution Time (ms)")
    plt.title("K-mer Search: Serial vs Parallel Execution Time")

    plt.xticks(threads)
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()

def plot_speedup(results, threads):

    plt.figure(figsize=(8, 5))

    for algorithm in ["naive", "array", "hash"]:

        speedups = [ results[algorithm][t]["speedup"]  for t in threads ]

        plt.plot( threads, speedups, marker="o", linewidth=2, label=f"{algorithm.capitalize()} Speedup" )

        for t, speedup in zip(threads, speedups):
            plt.text( t, speedup + 0.05, f"{speedup:.2f}x", ha="center")

    plt.plot( threads, threads, linestyle="--", linewidth=1.5, label="Ideal Speedup")

    plt.xlabel("Number of Threads")
    plt.ylabel("Speedup")
    plt.title("K-mer Search: Speedup vs Number of Threads")

    plt.xticks(threads)
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()
    
    
CSV_PATH = os.path.join(
    PROJECT_ROOT,
    "python/Lab/results",
    "threads_results.csv"
)
results = load_kmer_cache_results(CSV_PATH)
threads= [2 , 4, 8, 16]
plot_execution_time(results, threads)
plot_speedup(results, threads)







