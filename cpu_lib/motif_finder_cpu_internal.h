#ifndef MOTIF_FINDER_MOTIF_FINDER_CPU_INTERNAL_H
#define MOTIF_FINDER_MOTIF_FINDER_CPU_INTERNAL_H

#include <string>
#include <vector>

struct CpuInternalParams {
    std::vector<std::string> sequences;
    bool search_complementary { false };
    int thread_count { 0 };
};

void find_motifs_internal_cpu(std::vector<uint16_t> &out_motif_weights, const CpuInternalParams &params);

#endif // MOTIF_FINDER_MOTIF_FINDER_CPU_INTERNAL_H
