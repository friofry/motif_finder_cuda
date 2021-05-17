#include "match_functions.h"

#include <hash_conversions.h>

bool hash_matches_motif(uint32_t seq_hash, uint32_t motif_hash)
{
    return (seq_hash & motif_hash) == seq_hash;
}

bool hash_matches_motif(uint32_t seq_hash, uint32_t motif_hash, bool complementary)
{
    if (hash_matches_motif(seq_hash, motif_hash)) {
        return true;
    }

    if (complementary) {
        uint32_t compl_seq_hash = to_compl_hash_reverse(seq_hash);
        if (hash_matches_motif(compl_seq_hash, motif_hash)) {
            return true;
        }
    }
    return false;
}
