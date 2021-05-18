#ifndef MOTIF_FINDER_SEQUENCE_HASHES_UTILS_H
#define MOTIF_FINDER_SEQUENCE_HASHES_UTILS_H

#include <sequence_hashes.h>

int nonempty_sequence_count(const SequenceHashes &sequences_hashes);
void remove_motif_hashes(SequenceHashes &sequences_hashes, uint32_t latest_motif, bool complementary);
void mask_motif_hashes(SequenceHashes &sequences_hashes, uint32_t latest_motif, bool complementary);

#endif // MOTIF_FINDER_SEQUENCE_HASHES_UTILS_H
