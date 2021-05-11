#ifndef MOTIF_FINDER_SEQUENCES_TO_NUMS_H
#define MOTIF_FINDER_SEQUENCES_TO_NUMS_H

#include <sequence_nums.h>
#include <string>
#include <vector>

SequenceNums sequences_to_nums(const std::vector<std::string> &sequences);
std::vector<std::string> nums_to_sequences(const SequenceNums &seq_nums);

#endif // MOTIF_FINDER_SEQUENCES_TO_NUMS_H
