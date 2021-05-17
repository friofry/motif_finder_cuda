#ifndef MOTIF_FINDER_CREATE_STAT_MODEL_H
#define MOTIF_FINDER_CREATE_STAT_MODEL_H

#include <vector>
#include <string>

#include <stat_model.h>

StatModelPtr create_stat_model(bool use_real_nucl_frequences, uint32_t markov_level, const std::vector<std::string> &sequences, bool complementary);

#endif // MOTIF_FINDER_CREATE_STAT_MODEL_H
