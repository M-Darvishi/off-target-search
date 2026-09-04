#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "K_mers_core.h"
#include "K_mers_naive.h"
#include "K_mers_index_builder_serial.h"
#include "K_mers_index_builder_parallel.h"
#include "Sliding_window.h" 
#include "Hit.h"
#include "GenomeLoader.h"  

namespace py = pybind11;

PYBIND11_MODULE(offtarget_search, m) {
    m.doc() = "Off-target search";

    py::class_<Hit>(m, "Hit")
        .def_readonly("position", &Hit::position)
        .def_readonly("mismatches", &Hit::mismatches)
        .def("__repr__", [](const Hit& h) {
        return "<Hit position=" + std::to_string(h.position) +
            ", mismatches=" + std::to_string(h.mismatches) + ">";
            });

    py::class_<SlidingwindowHits>(m, "SlidingwindowHits")
    .def_readonly("hits", &SlidingwindowHits::hits)
    .def_readonly("search_time", &SlidingwindowHits::search_time)
    .def("repr", [](const SlidingwindowHits& r) {
        return "<SlidingwindowHits hits=" +
               std::to_string(r.hits.size()) +
               ", search_time=" +
               std::to_string(r.search_time) +
               "ms>";
    });

    py::class_<KmerHits>(m, "KmerHits")
        .def_readonly("hits", &KmerHits::hits)
        .def_readonly("index_time", &KmerHits::index_time)
        .def_readonly("collect_time", &KmerHits::collect_time)
        .def_readonly("search_time", &KmerHits::search_time)
        .def_readonly("candidate_count", &KmerHits::candidate_count)
        .def("__repr__", [](const KmerHits& r) {
        return "<KmerHits hits=" + std::to_string(r.hits.size()) +
            ", index_time=" + std::to_string(r.index_time) +
            "ms, collect_time=" + std::to_string(r.collect_time) +
            "ms, search_time=" + std::to_string(r.search_time) +
            "ms, candidates=" + std::to_string(r.candidate_count) + ">";
            });


    m.def("load_fasta", &loadFasta,
        "Load genome from FASTA file",
        py::arg("file_path"));

    m.def("sliding_window_serial", &slidingWindow_serial,
        "Serial sliding window search (baseline)",
        py::arg("genome"), py::arg("guide"), py::arg("max_mismatch"));

    m.def("sliding_window_parallel", &slidingWindow_parallel,
        "Parallel sliding window search",
        py::arg("genome"), py::arg("guide"), py::arg("max_mismatch"), py::arg("threads"));


    m.def("kmer_serial_array", &k_mers_serial_array,
        "Serial search for off-target sites",
        py::arg("genome"), py::arg("guide"), py::arg("max_mismatch"));

    m.def("kmer_serial_hash", &k_mers_serial_hash,
        "Serial search for off-target sites",
        py::arg("genome"), py::arg("guide"), py::arg("max_mismatch"));


    m.def("kmer_parallel_array", &k_mers_parallel_array,
        "Parallel search using array index",
        py::arg("genome"), py::arg("guide"), py::arg("max_mismatch"), py::arg("threads"));

    m.def("kmer_parallel_hash", &k_mers_parallel_hash,
        "Parallel search using hash index",
        py::arg("genome"), py::arg("guide"), py::arg("max_mismatch"), py::arg("threads"));


    m.def("kmer_build_index_array", &build_index_serial_array,
        "Build index array serial",
        py::arg("genome"), py::arg("k"));

    m.def("kmer_build_index_hash", &build_index_serial_hash,
        "Build index hash serial",
        py::arg("genome"), py::arg("k"));

    
    m.def("kmer_naive_serial", &k_mers_naive_serial,
        "Serial K-mer naive search",
        py::arg("genome"),
        py::arg("guide"),
        py::arg("max_mismatch"),
        py::arg("k"));

    m.def("kmer_naive_parallel", &k_mers_naive_parallel,
        "Parallel K-mer naive search",
        py::arg("genome"),
        py::arg("guide"),
        py::arg("max_mismatch"),
        py::arg("k"),
        py::arg("threads"));


    m.def("clear_cache", &clear_kmer_cache,
        "Clear the internal cache");
}
