#ifndef MOTIF_FINDER_MATCH_FUNCTIONS_H
#define MOTIF_FINDER_MATCH_FUNCTIONS_H

#include <cstdint>


bool hash_matches_motif(uint32_t seq_hash, uint32_t motif_hash);
bool hash_matches_motif(uint32_t seq_hash, uint32_t motif_hash, bool complementary);

#endif // MOTIF_FINDER_MATCH_FUNCTIONS_H
