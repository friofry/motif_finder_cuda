#ifndef MOTIF_FINDER_SEQUENCE_NUMS_H
#define MOTIF_FINDER_SEQUENCE_NUMS_H

#include <vector>
#include <cstdint>

struct SequenceNums {
    uint32_t count {0};
    std::vector<uint8_t> nums;
    std::vector<uint32_t> lengths;
    std::vector<uint32_t> seq_begins;
};

#endif // MOTIF_FINDER_SEQUENCE_NUMS_H
