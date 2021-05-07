#include "motif_finder_gpu_internal.h"
#include <thread>

#include <config.h>
#include <safe_counter.h>
#include <timer.h>
#include <sequences_to_hashes.h>
#include <hash_conversions.h>

#include <gpu_info.h>
#include <gpu_range_motif_finder_internal.h>

using namespace std;


namespace {
void finder_thread_internal(std::vector<uint16_t> &out_motif_weights,
                            const GpuRangeMotifFinderInternal &finder,
                            uint32_t motif_range_size,
                            SafeCounter &motifs_counter)
{
    std::vector<uint16_t> motif_occurrences(motif_range_size, 0);

    while (true) {
        auto range_info = motifs_counter.get_and_increment_range_info(motif_range_size);
        if (range_info.start >= range_info.end) {
            break;
        }

        uint32_t count = range_info.end - range_info.start;
        motif_occurrences.resize(count, 0);
        finder.get_occurrence_in_range(range_info.start, range_info.end, motif_occurrences);
        std::copy(motif_occurrences.begin(), motif_occurrences.end(), out_motif_weights.begin() + range_info.start);
    }
}
} // namespace

void find_motifs_internal_gpu(std::vector<uint16_t> &out_motif_weights, const GpuInternalParams &params)
{
    out_motif_weights.resize(TOTAL_MOT, 0);

    // Parse sequences
    SequenceHashes hashes_info = sequences_to_hashes(params.sequences, params.search_complementary);

    // Prepare motif finder
    SequenceHashesParams finder_params;
    finder_params.count = params.sequences.size();
    finder_params.lengths = hashes_info.lengths;
    finder_params.seq_begins = hashes_info.seq_begins;
    finder_params.hashes = hashes_info.hashes;

    GpuCudaParams cuda_params;
    cuda_params.unified_memory = params.unified_memory;
    cuda_params.motif_range_size = params.motif_range_size;
    cuda_params.threads_per_block = params.threads_per_block;

    GpuRangeMotifFinderInternal finder(finder_params, cuda_params);

    // Run in parallel
    SafeCounter motifs_counter(TOTAL_MOT);
    vector<thread> threads;
    int thread_count = params.gpu_count ? params.gpu_count : gpu_count();
    for (int i = 0; i < thread_count; i++) {
        threads.push_back(std::thread(finder_thread_internal,
                                      ref(out_motif_weights),
                                      ref(finder),
                                      cuda_params.motif_range_size,
                                      ref(motifs_counter)));
    }

    // Wait for finished
    for (int i = 0; i < thread_count; i++) {
        threads[i].join();
    }
}
