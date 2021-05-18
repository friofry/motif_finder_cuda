#include "sequences_to_hashes.h"

#include <config.h>

#include <hash_conversions.h>

using namespace std;

SequenceHashes sequences_to_hashes(const std::vector<std::string> &sequences, bool complementary)
{
    SequenceHashes ret;
    ret.count = sequences.size();
    for (const auto &s : sequences) {
        int hash_len = s.size() - MOTIV_LEN + 1;

        ret.seq_begins.push_back(ret.hashes.size());

        // Forward strand hashes
        for (int j = 0; j < hash_len; j++) {
            uint32_t hash = string_to_hash(&s[j]);
            ret.hashes.push_back(hash);
        }

        // Complementary strand hashes
        if (complementary) {
            for (int j = 0; j < hash_len; j++) {
                uint32_t hash = string_to_hash_compl(&s[j]);
                ret.hashes.push_back(hash);
            }
        }
        ret.lengths.push_back(complementary ? hash_len * 2 : hash_len);
    }
    return ret;
}

std::vector<std::string> hashes_to_sequences(const SequenceHashes &seq_hashes, bool complementary)
{
    std::vector<std::string> ret;

    for (uint32_t i = 0; i < seq_hashes.count; i++) {
        string seq;
        uint32_t start = seq_hashes.seq_begins[i];
        uint32_t length = complementary ? seq_hashes.lengths[i] / 2 : seq_hashes.lengths[i];

        for (uint32_t j = 0; j < length; j++) {
            auto s = hash_to_string(seq_hashes.hashes[start + j]);
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
                auto s = compl_hash_to_string(seq_hashes.hashes[start + j + length]);
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