import csv
import os
import sys

import matplotlib.pyplot as plt

PROJECT_ROOT = os.path.abspath("../../..")
sys.path.append(PROJECT_ROOT)


def load_results_from_csv(csv_path):

    results = {}
    with open(csv_path, "r", newline="") as file:

        reader = csv.DictReader(file)
        for row in reader:

            threads = int(row["threads"])

            result = {
                "mean": float(row["mean_ms"]),
                "median": float(row["median_ms"]),
                "speedup": float(row["speedup"]),
                "hits": int(row["hits"])
            }

            if threads == 1:
                results["serial"] = result
            else:
                results[threads] = result

    return results


def plot_execution_time(results, threads):

    serial_mean = results["serial"]["mean"]

    parallel_mean_times = [
        results[t]["mean"]
        for t in threads
    ]

    plt.figure(figsize=(8, 5))

    plt.plot(
        threads,
        parallel_mean_times,
        marker="o",
        linewidth=2,
        label="Parallel"
    )

    plt.axhline(
        y=serial_mean,
        linewidth=2,
        label="Serial"
    )

    plt.xlabel("Number of Threads")
    plt.ylabel("Execution Time (ms)")
    plt.title("Sliding Window: Execution Time")

    plt.xticks(threads)
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


def plot_speedup(results, threads):

    speedups = [
        results[t]["speedup"]
        for t in threads
    ]

    plt.figure(figsize=(8, 5))

    plt.plot(
        threads,
        speedups,
        marker="o",
        linewidth=2,
        label="Measured Speedup"
    )

    plt.plot(
        threads,
        threads,
        linestyle="--",
        linewidth=1.5,
        label="Ideal Speedup"
    )

    plt.xlabel("Number of Threads")
    plt.ylabel("Speedup")
    plt.title("Sliding Window: Speedup vs Number of Threads")

    plt.xticks(threads)
    plt.grid(True)
    plt.legend()

    for t, speedup in zip(threads, speedups):
        plt.text(
            t,
            speedup + 0.05,
            f"{speedup:.2f}x",
            ha="center"
        )

    plt.tight_layout()
    plt.show()
    
 
 
CSV_PATH = os.path.join(
    PROJECT_ROOT,
    "python/Lab/results",
    "sliding_window_results.csv"
)
   
    
results = load_results_from_csv(CSV_PATH)
threads= [2 , 4, 8, 16]
plot_execution_time(results, threads)
plot_speedup(results, threads)




    
