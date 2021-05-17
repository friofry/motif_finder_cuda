#pragma once
#include <vector>
#include <cstdint>

#include <stat_model.h>
#include "motif_data.h"

// найти мотив с самым высоким значением хи2 (не меньше чем min_chi2)
// @param weights - веса мотивов
// return индекс самого значимого мотива, либо -1
MotifData find_most_important_motif(
    const std::vector<uint16_t> &motif_weights,
    const std::vector<uint32_t> &motif_hashes,
    const StatModel &stat_model,
    double min_presence,
    bool use_binom_instead_chi2,
    const std::vector<uint16_t> &prev_weights,
    std::vector<double> &binomial_prob);
