import statistics


def benchmark_serial(ot, genome, guide, max_mismatch, repeats=10, warmups=2):

    times = []
    hit_count = None

    for _ in range(warmups):
        ot.sliding_window_serial(
            genome,
            guide,
            max_mismatch
        )

    for _ in range(repeats):

        result = ot.sliding_window_serial(
            genome,
            guide,
            max_mismatch
        )

        times.append(result.search_time )

        current_hit_count = len(result.hits)

        if hit_count is None:
            hit_count = current_hit_count
        elif current_hit_count != hit_count:
            raise RuntimeError("Different hit count detected!")

    return {
        "times": times,
        "mean": statistics.mean(times),
        "median": statistics.median(times),
        "hits": hit_count
    }


def benchmark_parallel(
    ot,
    genome,
    guide,
    max_mismatch,
    threads,
    repeats=10,
    warmups=2
):

    times = []
    hit_count = None

    for _ in range(warmups):
        ot.sliding_window_parallel(
            genome,
            guide,
            max_mismatch,
            threads
        )

    for _ in range(repeats):

        result = ot.sliding_window_parallel(
            genome,
            guide,
            max_mismatch,
            threads
        )

        times.append(result.search_time)

        current_hit_count = len(result.hits)

        if hit_count is None:
            hit_count = current_hit_count
        elif current_hit_count != hit_count:
            raise RuntimeError(
                f"Different hit count detected for {threads} threads!"
            )

    return {
        "times": times,
        "mean": statistics.mean(times),
        "median": statistics.median(times),
        "hits": hit_count
    }
