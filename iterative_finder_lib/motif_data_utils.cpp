//
// Created by andrey on 13.06.2021.
//

#include "motif_data_utils.h"

std::vector<uint32_t> get_hashes(const std::vector<MotifData> &motif_data)
{
    std::vector<uint32_t> results;
    results.reserve(motif_data.size());
    for (const auto &d: motif_data) {
        results.push_back(d.hash);
    }
    return results;
}
