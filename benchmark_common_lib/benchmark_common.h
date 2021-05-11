#ifndef MOTIF_FINDER_BENCHMARK_COMMON_H
#define MOTIF_FINDER_BENCHMARK_COMMON_H

#include <string>
#include <vector>
#include <stat_model.h>

void save_to_file(const std::vector<std::string> &sequences, const char *filename);

// Prepares stat model
// @param model - stat model. one of bernulli, ideal, markov1, morkov2, markov3
// @param sequences - nucleotide sequences
// @param complementary - true - use both fwd and compl strands, false - use only fwd
StatModelPtr create_stat_model(const std::string &model, const std::vector<std::string> &sequences, bool complementary);

#endif //MOTIF_FINDER_BENCHMARK_COMMON_H

