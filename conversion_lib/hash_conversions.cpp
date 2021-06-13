#include "hash_conversions.h"

#include <iostream>

#include <config.h>
#include <letter_conversions.h>

using namespace std;

uint32_t idx_to_hash(uint32_t idx)
{
    uint32_t mult = 1;
    uint32_t ret = 0;
    for (int j = MOTIV_LEN - 1; j >= 0; j--) {
        ret += (idx % ALPH_SIZE + 1) * mult;
        mult *= HASH_BASE;
        idx /= ALPH_SIZE;
    }
    return ret;
}

uint32_t string_to_hash(const char *s)
{
    uint32_t hash = 0;
    uint32_t mult = 1;
    for (int j = MOTIV_LEN - 1; j >= 0; j--) {
        uint8_t numeric_letter = convert_to_code(s[j]);
        hash += numeric_letter * mult;
        mult *= HASH_BASE;
    }
    return hash;
}

uint32_t string_to_hash_compl(const char *s)
{
    uint32_t hash = 0;
    uint32_t mult = 1;
    for (int j = 0; j < MOTIV_LEN; j++) {
        uint8_t numeric_letter = convert_to_code(s[j]);
        hash += numeric_letter * mult;
        mult *= HASH_BASE;
    }
    return to_compl_hash(hash);
}

uint32_t to_compl_hash(uint32_t hash)
{
    // swap AT and GC
    return ((hash & 0x55555555) << 1) | ((hash & 0xAAAAAAAA) >> 1);
}

uint32_t to_compl_hash_reverse(uint32_t v)
{
    // swap odd and even bits
    v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
    // swap consecutive pairs
    v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
    // swap nibbles ...
    v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
    // swap bytes
    v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
    // swap 2-byte long pairs
    v = (v >> 16) | (v << 16);

    return ((v & 0x33333333) << 2) | ((v & 0xCCCCCCCC) >> 2);
}

// hash to readable string
std::string hash_to_string(uint32_t hash)
{
    uint8_t mm[MOTIV_LEN];
    // from end to start
    for (int j = MOTIV_LEN - 1; j >= 0; j--) {
        mm[j] = hash % HASH_BASE;
        hash /= HASH_BASE;
    }

    std::string result;
    for (int i = 0; i < MOTIV_LEN; i++) {
        result += hash_letter_to_symbol(mm[i]);
    }
    return result;
}

std::string compl_hash_to_string(uint32_t hash)
{
    hash = to_compl_hash(hash);
    uint8_t mm[MOTIV_LEN];

    for (int j = 0; j < MOTIV_LEN; j++) {
        mm[j] = hash % HASH_BASE;
        hash /= HASH_BASE;
    }

    std::string result;
    for (int i = 0; i < MOTIV_LEN; i++) {
        result += hash_letter_to_symbol(mm[i]);
    }
    return result;
}

uint32_t idx_to_motif_prefix_hash(uint32_t prefix_idx, uint32_t prefix_size)
{
    uint32_t current_motif_hash = 0;
    uint32_t i = 0;
    while (prefix_idx > 0) {
        if (i > prefix_size) {
            cout << "error in idx_to_motif_prefix_hash" << endl;
            break;
        }

        uint32_t symbol_hash = prefix_idx % 16;
        if (symbol_hash == 0) {
            return 0;
        }

        current_motif_hash |= symbol_hash << (i * 4);
        prefix_idx /= 16;
        i++;
    }
    return current_motif_hash << (32 - prefix_size * 4);
}

uint32_t olig_id_to_prefix_hash(uint32_t olig_idx, uint32_t prefix_size)
{
    uint32_t current_hash_mask = 0;
    for (uint32_t i = 0; i < prefix_size; i++) {
        uint32_t nucl_hash = 1 << (olig_idx % 4);
        current_hash_mask |= nucl_hash << (i * 4);
        olig_idx /= 4;
    }

    return current_hash_mask << (32 - prefix_size * 4);
}