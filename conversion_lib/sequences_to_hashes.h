#ifndef MOTIF_FINDER_SEQUENCES_TO_HASHES_H
#define MOTIF_FINDER_SEQUENCES_TO_HASHES_H

#include <vector>
#include <string>
#include <SequenceHashes.h>

SequenceHashes sequences_to_hashes(const std::vector<std::string> &sequences, bool complementary);
std::vector<std::string> hashes_to_sequences(const SequenceHashes &seq_hashes, bool complementary);

#endif // MOTIF_FINDER_SEQUENCES_TO_HASHES_H
