//
// Created by andrey on 13.06.2021.
//

#include "filter_shifted_motifs.h"
#include "hash_conversions.h"

namespace {

bool matches(uint32_t hash, uint32_t motif)
{
    auto left_shift = 0xF | (motif << 4);
    auto right_shift = 0xF0000000 | (motif >> 4);
    return hash == left_shift || hash == right_shift;
}

bool matches(uint32_t hash, uint32_t motif, bool complementary)
{
    return matches(hash, motif) || (complementary && matches(hash, to_compl_hash_reverse(motif)));
}

} // namespace

bool hash_matches(uint32_t hash, const std::vector<uint32_t> &hashes, bool complementary)
{
    for (const auto &motif : hashes) {
        if (hash == motif) {
            return false;
        }
        if (matches(hash, motif, complementary)) {
            return true;
        }
    }
    return false;
}