#include "motif_finder_cpu_external.h"

#include <thread>

#include <safe_counter.h>
#include <timer.h>
#include <sequences_to_hashes.h>

#include <cpu_range_motif_finder.h>

using namespace std;

namespace {
void finder_thread_external(std::vector<uint16_t> &out_motif_weights,
                            const std::vector<uint32_t> &motif_hashes,
                            const CpuRangeMotifFinder &finder,
                            SafeCounter &motifs_counter)
{
    int chunk_size = MOT_PER_CHUNK / 10;

    std::vector<uint32_t> motif_hashes_range(chunk_size, 0);
    std::vector<uint16_t> motif_occurrences(chunk_size, 0);

    while (true) {
        auto range_info = motifs_counter.get_and_increment_range_info(chunk_size);
        if (range_info.start >= range_info.end) {
            break;
        }

        uint32_t count = range_info.end - range_info.start;
        motif_hashes_range.resize(count, 0);
        motif_occurrences.resize(count, 0);

        std::copy(motif_hashes.begin() + range_info.start, motif_hashes.begin() + range_info.end, motif_hashes_range.begin());
        finder.get_occurrence_in_range(motif_hashes_range, motif_occurrences);
        std::copy(motif_occurrences.begin(), motif_occurrences.end(), out_motif_weights.begin() + range_info.start);
    }
}
} // namespace

void find_motifs_external_cpu(std::vector<uint16_t> &out_motif_weights, const CpuExternalParams &params)
{
    out_motif_weights.resize(params.motif_hashes.size(), 0);
    SequenceHashes hashes_info = sequences_to_hashes(params.sequences, params.search_complementary);
    SequenceHashesParams finder_params;
    finder_params.count = params.sequences.size();
    finder_params.lengths = hashes_info.lengths;
    finder_params.seq_begins = hashes_info.seq_begins;
    finder_params.hashes = hashes_info.hashes;

    CpuRangeMotifFinder finder(finder_params);
    SafeCounter motifs_counter(params.motif_hashes.size());
    vector<thread> threads;
    int thread_count = params.thread_count ? params.thread_count : thread::hardware_concurrency();
    for (int i = 0; i < thread_count; i++) {
        threads.push_back(std::thread(finder_thread_external,
                                      ref(out_motif_weights),
                                      ref(params.motif_hashes),
                                      ref(finder),
                                      ref(motifs_counter)));
    }

    for (int i = 0; i < params.thread_count; i++) {
        threads[i].join();
    }
}
