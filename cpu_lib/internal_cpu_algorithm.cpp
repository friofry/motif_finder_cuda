#include "internal_cpu_algorithm.h"

#include <thread>

#include <config.h>
#include <hash_conversions.h>
#include <run_parallel.h>
#include <safe_counter.h>

#include "cpu_common_algorithm.h"

void internal_cpu_algorithm(const SequenceHashes &sequence_hashes,
                            std::vector<uint16_t> &out_motif_weights,
                            uint32_t threads)
{
    out_motif_weights.resize(TOTAL_MOT, 0);
    threads = (threads > 0) ? threads : std::thread::hardware_concurrency();

    SafeCounter motifs_counter(TOTAL_MOT);

    run_parallel(threads, [&] (uint32_t thread_id) {
        while (true) {
            auto range = motifs_counter.get_and_increment_range_info(MOT_PER_CHUNK/100);
            if (range.count() == 0) {
                break;
            }
            printf("%d %d\n", thread_id, range.start/1000000);
            for (uint32_t i = range.start; i < range.end; i++) {
                out_motif_weights[i] = get_motif_occurrence(idx_to_hash(i), sequence_hashes);
            }
        }
    });
}

