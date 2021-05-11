#ifndef MOTIF_FINDER_SEQUENCE_HASHES_H
#define MOTIF_FINDER_SEQUENCE_HASHES_H
#include <vector>
#include <cstdint>

struct SequenceHashes {
    uint32_t count {0};
    std::vector<uint32_t> hashes;
    std::vector<uint32_t> lengths;
    std::vector<uint32_t> seq_begins;
};

#endif // MOTIF_FINDER_SEQUENCE_HASHES_H
