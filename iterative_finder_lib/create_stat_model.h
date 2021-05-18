#ifndef MOTIF_FINDER_CREATE_STAT_MODEL_H
#define MOTIF_FINDER_CREATE_STAT_MODEL_H

#include <string>
#include <vector>

#include <stat_model.h>

StatModelPtr create_stat_model(bool use_real_nucl_frequences,
                               uint32_t markov_level,
                               const std::vector<std::string> &sequences,
                               bool complementary,
                               bool use_binom_instead_of_chi2);

#endif // MOTIF_FINDER_CREATE_STAT_MODEL_H
