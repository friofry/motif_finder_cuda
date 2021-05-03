#ifndef MOTIF_FINDER_HASH_CONVERSIONS_H
#define MOTIF_FINDER_HASH_CONVERSIONS_H

#include <stdint.h>
#include <string>

uint32_t idx_to_hash(uint32_t idx);

uint32_t string_to_hash(const char *s);
uint32_t string_to_hash_compl(const char *s);

std::string hash_to_string(uint32_t hash);
std::string compl_hash_to_string(uint32_t hash);

uint32_t to_compl_hash(uint32_t hash);

#endif //MOTIF_FINDER_HASH_CONVERSIONS_H

