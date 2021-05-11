#ifndef MOTIF_FINDER_EXTERNAL_CPU_ALGORITHM_H
#define MOTIF_FINDER_EXTERNAL_CPU_ALGORITHM_H

#include <SequenceHashes.h>
#include <config.h>
#include <cstdint>

// Get motif occurrences and store to out_motif_weights
// @param motif_hashes - motif hashes
// @param sequence_hashes - sequences hashes
// @param out_motif_weights - [output] array for motif occurrences
// @param thread_count - number of threads. if 0 - then use maximum available threads
void external_cpu_algorithm(const std::vector<uint32_t> &motif_hashes,
                            const SequenceHashes &sequence_hashes,
                            std::vector<uint16_t> &out_motif_weights,
                            uint32_t thread_count = 0);

#endif //MOTIF_FINDER_EXTERNAL_CPU_ALGORITHM_H
