import statistics

def benchmark_naive_serial(ot, genome, guide, max_mismatch, k, repeats=10, warmups=2):

    index_times = []
    search_times = []
    total_times = []
    hit_count = None
    candidate_count = None

    for _ in range(warmups):
        ot.kmer_naive_serial(  genome, guide,  max_mismatch,  k )

    for _ in range(repeats):

        result = ot.kmer_naive_serial( genome, guide, max_mismatch, k )

        index_times.append(result.index_time)
        search_times.append(result.search_time)
        total_times.append(result.search_time + result.index_time + result.collect_time)

        current_hit_count = len(result.hits)
        current_candidate_count = result.candidate_count

        if hit_count is None:
            hit_count = current_hit_count
            candidate_count = current_candidate_count

        elif current_hit_count != hit_count:
            raise RuntimeError(
                "Different hit count detected for Naive serial!"
            )

        elif current_candidate_count != candidate_count:
            raise RuntimeError(
                "Different candidate count detected for Naive serial!"
            )

    return {
        "index_times": index_times,
        "search_times": search_times,
        "index_mean": statistics.mean(index_times),
        "search_mean": statistics.mean(search_times),
        "index_median": statistics.median(index_times),
        "search_median": statistics.median(search_times),
        "mean": statistics.mean(total_times),
        "median": statistics.median(total_times),
        "hits": hit_count,
        "candidates": candidate_count
    }


def benchmark_naive_parallel(ot, genome, guide, max_mismatch, k, threads, repeats=10, warmups=2):

    index_times = []
    search_times = []
    total_times = []
    hit_count = None
    candidate_count = None

    for _ in range(warmups):
        ot.kmer_naive_parallel(  genome, guide,  max_mismatch, k,  threads)

    for _ in range(repeats):

        result = ot.kmer_naive_parallel( genome, guide, max_mismatch, k, threads )

        index_times.append(result.index_time)
        search_times.append(result.search_time)
        total_times.append(result.search_time + result.index_time + result.collect_time)

        current_hit_count = len(result.hits)
        current_candidate_count = result.candidate_count

        if hit_count is None:
            hit_count = current_hit_count
            candidate_count = current_candidate_count

        elif current_hit_count != hit_count:
            raise RuntimeError(
                f"Different hit count detected for Naive "
                f"parallel with {threads} threads!"
            )

        elif current_candidate_count != candidate_count:
            raise RuntimeError(
                f"Different candidate count detected for Naive "
                f"parallel with {threads} threads!"
            )

    return {
        "index_times": index_times,
        "search_times": search_times,
        "index_mean": statistics.mean(index_times),
        "search_mean": statistics.mean(search_times),
        "index_median": statistics.median(index_times),
        "search_median": statistics.median(search_times),
        "mean": statistics.mean(total_times),
        "median": statistics.median(total_times),
        "hits": hit_count,
        "candidates": candidate_count
    }


def benchmark_array_serial(ot, genome, guide, max_mismatch, k, repeats=10, warmups=2, clear_cache=False):

    index_times = []
    search_times = []
    collect_times =[]
    total_times = []
    hit_count = None
    candidate_count = None

    for _ in range(warmups):
         
        ot.kmer_serial_array( genome, guide,  max_mismatch )
        if(clear_cache):
            ot.clear_cache()

    ot.clear_cache()
    for _ in range(repeats):

        if(clear_cache):
            ot.clear_cache()
        result = ot.kmer_serial_array( genome, guide, max_mismatch )
        
        index_times.append(result.index_time)
        search_times.append(result.search_time)
        collect_times.append(result.collect_time)
        total_times.append(result.search_time + result.index_time + result.collect_time)

        current_hit_count = len(result.hits)
        current_candidate_count = result.candidate_count

        if hit_count is None:
            hit_count = current_hit_count
            candidate_count = current_candidate_count

        elif current_hit_count != hit_count:
            raise RuntimeError(
                "Different hit count detected for Array serial!"
            )

        elif current_candidate_count != candidate_count:
            raise RuntimeError(
                "Different candidate count detected for Array serial!"
            )

        if(clear_cache):
            ot.clear_cache()

    return {
        "index_times": index_times,
        "search_times": search_times,
        "collect_times": collect_times,
        "index_mean": statistics.mean(index_times),
        "search_mean": statistics.mean(search_times),
        "collect_mean": statistics.mean(collect_times),
        "index_median": statistics.median(index_times),
        "search_median": statistics.median(search_times),
        "collect_median": statistics.median(collect_times),
        "mean": statistics.mean(total_times),
        "median": statistics.median(total_times),
        "hits": hit_count,
        "candidates": candidate_count
    }
 

def benchmark_array_parallel(ot, genome, guide, max_mismatch, k, threads, repeats=10, warmups=2, clear_cache=False):

    index_times = []
    search_times = []
    collect_times =[]
    total_times = []
    hit_count = None
    candidate_count = None

    for _ in range(warmups): 
 
        ot.kmer_parallel_array( genome,  guide,  max_mismatch,  threads )
        if(clear_cache):
            ot.clear_cache()
        
    ot.clear_cache()
    for _ in range(repeats):
 
        result = ot.kmer_parallel_array( genome, guide, max_mismatch, threads )

        index_times.append(result.index_time)
        search_times.append(result.search_time)
        collect_times.append(result.collect_time)
        total_times.append(result.search_time + result.index_time + result.collect_time)

        current_hit_count = len(result.hits)
        current_candidate_count = result.candidate_count

        if hit_count is None:
            hit_count = current_hit_count
            candidate_count = current_candidate_count

        elif current_hit_count != hit_count:
            raise RuntimeError(
                f"Different hit count detected for Array "
                f"parallel with {threads} threads!"
            )

        elif current_candidate_count != candidate_count:
            raise RuntimeError(
                f"Different candidate count detected for Array "
                f"parallel with {threads} threads!"
            )
        
        if(clear_cache):
            ot.clear_cache()

    return {
        "index_times": index_times,
        "search_times": search_times,
        "collect_times": collect_times,
        "index_mean": statistics.mean(index_times),
        "search_mean": statistics.mean(search_times),
        "collect_mean": statistics.mean(collect_times),
        "index_median": statistics.median(index_times),
        "search_median": statistics.median(search_times),
        "collect_median": statistics.median(collect_times),
        "mean": statistics.mean(total_times),
        "median": statistics.median(total_times),
        "hits": hit_count,
        "candidates": candidate_count
    }
 

def benchmark_hash_serial(ot, genome, guide, max_mismatch, k, repeats=10, warmups=2, clear_cache=False):

    index_times = []
    search_times = []
    collect_times =[]
    total_times = []
    hit_count = None
    candidate_count = None

    for _ in range(warmups):
 
        ot.kmer_serial_hash( genome, guide, max_mismatch )
        if(clear_cache):
            ot.clear_cache()

    ot.clear_cache()
    for _ in range(repeats):
        
        result = ot.kmer_serial_hash( genome, guide, max_mismatch )

        index_times.append(result.index_time)
        search_times.append(result.search_time)
        collect_times.append(result.collect_time)
        total_times.append(result.search_time + result.index_time + result.collect_time)

        current_hit_count = len(result.hits)
        current_candidate_count = result.candidate_count

        if hit_count is None:
            hit_count = current_hit_count
            candidate_count = current_candidate_count

        elif current_hit_count != hit_count:
            raise RuntimeError(
                "Different hit count detected for Hash serial!"
            )

        elif current_candidate_count != candidate_count:
            raise RuntimeError(
                "Different candidate count detected for Hash serial!"
            )
        
        if(clear_cache):
            ot.clear_cache()

    return {
        "index_times": index_times,
        "search_times": search_times,
        "collect_times": collect_times,
        "index_mean": statistics.mean(index_times),
        "search_mean": statistics.mean(search_times),
        "collect_mean": statistics.mean(collect_times),
        "index_median": statistics.median(index_times),
        "search_median": statistics.median(search_times),
        "collect_median": statistics.median(collect_times),
        "mean": statistics.mean(total_times),
        "median": statistics.median(total_times),
        "hits": hit_count,
        "candidates": candidate_count
    }
 

def benchmark_hash_parallel(ot, genome, guide, max_mismatch, k, threads, repeats=10, warmups=2, clear_cache=False):

    index_times = []
    search_times = []
    collect_times =[]
    total_times = []
    hit_count = None
    candidate_count = None

    for _ in range(warmups):
 
        ot.kmer_parallel_hash( genome, guide,  max_mismatch, threads )
        if(clear_cache):
            ot.clear_cache()       
        
    ot.clear_cache()
    for _ in range(repeats):
         
        result = ot.kmer_parallel_hash( genome,  guide, max_mismatch, threads )

        index_times.append(result.index_time)
        search_times.append(result.search_time)
        collect_times.append(result.collect_time)
        total_times.append(result.search_time + result.index_time + result.collect_time)

        current_hit_count = len(result.hits)
        current_candidate_count = result.candidate_count

        if hit_count is None:
            hit_count = current_hit_count
            candidate_count = current_candidate_count

        elif current_hit_count != hit_count:
            raise RuntimeError(
                f"Different hit count detected for Hash "
                f"parallel with {threads} threads!"
            )

        elif current_candidate_count != candidate_count:
            raise RuntimeError(
                f"Different candidate count detected for Hash "
                f"parallel with {threads} threads!"
            )
        
        if(clear_cache):
            ot.clear_cache()

    return {
        "index_times": index_times,
        "search_times": search_times,
        "collect_times": collect_times,
        "index_mean": statistics.mean(index_times),
        "search_mean": statistics.mean(search_times),
        "collect_mean": statistics.mean(collect_times),
        "index_median": statistics.median(index_times),
        "search_median": statistics.median(search_times),
        "collect_median": statistics.median(collect_times),
        "mean": statistics.mean(total_times),
        "median": statistics.median(total_times),
        "hits": hit_count,
        "candidates": candidate_count
    }
