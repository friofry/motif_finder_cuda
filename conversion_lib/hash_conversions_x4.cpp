#include "hash_conversions_x4.h"

#include <math.h>

#include "hash_conversions.h"

#include <config.h>

using namespace std;

namespace {

const uint32_t HI_PART_15_MIN = 15*15*15*15;

uint32_t hash_to_idx_with_check(uint32_t hash, bool &out_ok)
{
    uint32_t idx = 0;
    uint32_t mult = 1;
    out_ok = true;
    for (int j = MOTIV_LEN - 1; j >= 0; j--) {
        if (hash % HASH_BASE == 0) {
            out_ok = false;
            return 0;
        }
        idx += (hash % HASH_BASE - 1) * mult;
        mult *= ALPH_SIZE;
        hash /= HASH_BASE;
    }
    return idx;
}

} // namespace

// Преобразование индекса в хэш
vector<uint32_t> calc_hash_x4()
{
    uint32_t total_idx_4 = 15 * 15 * 15 * 15;
    vector<uint32_t> result_hashes(total_idx_4, 0);

    for (uint32_t idx = 0; idx < total_idx_4; idx++) {
        uint32_t hash = idx_to_hash(idx) & 0xFFFF;
        result_hashes[idx] = hash;
    }

    return result_hashes;
}

vector<uint32_t> calc_idx_x4()
{
    uint32_t total_hashes_4 = 1 << 16;
    vector<uint32_t> result_indexes(total_hashes_4, 0);

    bool ok = true;
    for (uint32_t hash = 0; hash < total_hashes_4; hash++) {
        uint32_t idx = hash_to_idx_with_check(0x11110000 | hash, ok);
        if (!ok) {
            continue;
        }

        result_indexes[hash] = idx;
    }

    return result_indexes;
}

// посчитать хэш по индексу
uint32_t idx_to_hash_x4(uint32_t idx, const vector<uint32_t> &hash_x4)
{
    return hash_x4[idx % HI_PART_15_MIN] + (hash_x4[idx / HI_PART_15_MIN] << 16);
}

// посчитать индексу по хэшу
uint32_t hash_to_idx_x4(uint32_t hash, const vector<uint32_t> &idx_x4)
{
    return idx_x4[hash & 0xFFFF] + (idx_x4[hash >> 16] * HI_PART_15_MIN);
}
