#pragma once

struct Hit {
    int position;
    int mismatches;
};


struct KmerHits {
    std::vector<Hit> hits;

    size_t candidate_count;

    double index_time;

    double collect_time;

    double search_time;
};

struct SlidingwindowHits {
    std::vector<Hit> hits;

    double search_time;
};