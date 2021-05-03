#ifndef MOTIF_FINDER_MOTIFFINDERCPUEXTERNAL_H
#define MOTIF_FINDER_MOTIFFINDERCPUEXTERNAL_H

#include <vector>
#include <stdint.h>


struct CpuExternalParams {
    CpuExternalParams(std::vector<std::string> sequences, bool search_complementary,
                      const std::vector<uint32_t> &motif_hashes, int thread_count)
        : sequences(sequences)
        , search_complementary(search_complementary)
        , motif_hashes(motif_hashes)
        , thread_count(thread_count){}

    std::vector<std::string> sequences;
    bool search_complementary;
    const std::vector<uint32_t> &motif_hashes;
    int thread_count;
};

void find_motifs_external_cpu(std::vector<uint16_t> &out_motif_weights, const CpuExternalParams &params);

#endif //MOTIF_FINDER_MOTIFFINDERCPUEXTERNAL_H
