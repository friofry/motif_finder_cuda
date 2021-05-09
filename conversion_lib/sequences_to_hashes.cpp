#include "sequences_to_hashes.h"

#include <config.h>

#include <hash_conversions.h>

using namespace  std;

SequenceHashes sequencesToHashes(const std::vector<std::string> &sequences, bool complementary)
{
    SequenceHashes ret;
    ret.count = sequences.size();
    for (const auto &s: sequences) {
        int hashLen = s.size() - MOTIV_LEN + 1;

        ret.seq_begins.push_back(ret.hashes.size());

        // Forward strand hashes
        for (int j = 0; j < hashLen; j++) {
            uint32_t hash = string_to_hash(&s[j]);
            ret.hashes.push_back(hash);
        }

        // Complementary strand hashes
        if (complementary) {
            for (int j = 0; j < hashLen; j++) {
                uint32_t hash = string_to_hash_compl(&s[j]);
                ret.hashes.push_back(hash);
            }
        }
        ret.lengths.push_back(complementary ? hashLen * 2 : hashLen);
    }
    return ret;
}

std::vector<std::string> hashesToSequences(const SequenceHashes &seqHashes, bool complementary)
{
    std::vector<std::string> ret;

    for (uint32_t i = 0; i < seqHashes.count; i++) {
        string seq;
        uint32_t start = seqHashes.seq_begins[i];
        uint32_t length = complementary ? seqHashes.lengths[i] / 2 : seqHashes.lengths[i];

        for (uint32_t j = 0; j < length; j++) {
            auto s = hash_to_string(seqHashes.hashes[start + j]);
            if (j < length - 1) {
                seq += s[0];
            } else {
                seq += s;
                if (complementary) {
                    seq += " ";
                }
            }
        }

        if (complementary) {
            for (uint32_t j = 0; j < length; j++) {
                auto s = compl_hash_to_string(seqHashes.hashes[start + j + length]);
                if (j < length - 1) {
                    seq += s[0];
                } else {
                    seq += s;
                }
            }
        }
        ret.push_back(seq);
    }
    return ret;
}