#ifndef MOTIF_FINDER_CPU_COMMON_ALGORITHM_H
#define MOTIF_FINDER_CPU_COMMON_ALGORITHM_H

#include <SequenceHashes.h>
#include <cstdint>

uint32_t get_motif_occurrence(uint32_t motif_hash, const SequenceHashes &sequence_hashes);

#endif //MOTIF_FINDER_CPU_COMMON_ALGORITHM_H
