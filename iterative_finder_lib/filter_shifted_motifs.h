#ifndef MOTIF_FINDER_CPU_FILTER_SHIFTED_MOTIFS_H
#define MOTIF_FINDER_CPU_FILTER_SHIFTED_MOTIFS_H

#include <vector>
#include <cstdint>


/* Test if motifs that are equal to one of the motifs with a shift */
bool hash_matches(uint32_t hash, const std::vector<uint32_t> &hashes, bool complementary);

#endif // MOTIF_FINDER_CPU_FILTER_SHIFTED_MOTIFS_H
