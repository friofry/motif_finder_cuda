#ifndef MOTIF_FINDER_SEQUENCESPARAMS_H
#define MOTIF_FINDER_SEQUENCESPARAMS_H

#include <vector>

struct SequenceHashesParams {
    uint32_t count; // number of sequences
    std::vector<uint32_t> hashes; // subsequence hashes
    std::vector<uint32_t> lengths; // sequence lengths array
    std::vector<uint32_t> seq_begins; // sequence start indexes array
};
#endif // MOTIF_FINDER_SEQUENCESPARAMS_H
