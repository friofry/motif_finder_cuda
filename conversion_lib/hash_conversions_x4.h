#ifndef MOTIF_FINDER_HASH_CONVERSIONS_X4_H
#define MOTIF_FINDER_HASH_CONVERSIONS_X4_H

#include <vector>
#include <cstdint>

/// предрассчитать хэши для всех возможны 4-буквенных индексов
std::vector<uint32_t> calc_hash_x4();

// предрассчитать индекс для всех возможных 4-буквенных хэшей
std::vector<uint32_t> calc_idx_x4();

// посчитать хэш по индексу
uint32_t idx_to_hash_x4(uint32_t idx, const std::vector<uint32_t> &hash_x4);

// посчитать индексу по хэшу
uint32_t hash_to_idx_x4(uint32_t hash, const std::vector<uint32_t> &idx_x4);


#endif // MOTIF_FINDER_HASH_CONVERSIONS_X4_H
