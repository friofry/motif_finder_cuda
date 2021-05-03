#ifndef MOTIF_FINDER_SEQUENCEHASHES_H
#define MOTIF_FINDER_SEQUENCEHASHES_H

struct SequenceHashes {
    int count {0};
    std::vector<uint32_t> hashes;
    std::vector<uint32_t> lengths;
    std::vector<uint32_t> seq_begins;
};

#endif //MOTIF_FINDER_SEQUENCEHASHES_H
