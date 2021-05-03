//
// Created by Андрей Бочарников on 03.05.2021.
//

#include "hash_conversions.h"

#include <config.h>
#include <letter_conversions.h>

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

uint32_t string_to_hash(const char *s) {
    uint32_t hash = 0;
    uint32_t mult = 1;
    for (int j = MOTIV_LEN - 1; j >= 0; j--) {
        uint8_t numeric_letter = to_num(s[j]);
        hash += (1 << numeric_letter) * mult;
        mult *= HASH_BASE;
    }
    return hash;
}

uint32_t string_to_hash_compl(const char *s) {
    uint32_t hash = 0;
    uint32_t mult = 1;
    for (int j = 0; j < MOTIV_LEN; j++) {
        uint8_t numeric_letter = to_num(s[j]);
        hash += (1 << numeric_letter) * mult;
        mult *= HASH_BASE;
    }
    return to_compl_hash(hash);
}

uint32_t to_compl_hash(uint32_t hash)
{
    // swap AT and GC
    return ((hash & 0x55555555) << 1) | ((hash & 0xAAAAAAAA) >> 1);
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
