#ifndef MOTIF_FINDER_CPU_MOTIF_DATA_UTILS_H
#define MOTIF_FINDER_CPU_MOTIF_DATA_UTILS_H

#include <vector>
#include <cstdint>
#include "motif_data.h"

std::vector<uint32_t> get_hashes(const std::vector<MotifData> &motif_data);

#endif // MOTIF_FINDER_CPU_MOTIF_DATA_UTILS_H
