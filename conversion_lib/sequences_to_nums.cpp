#include "sequences_to_nums.h"
#include <letter_conversions.h>

using namespace  std;

SequenceNums sequences_to_nums(const std::vector<std::string> &sequences)
{
    SequenceNums ret;
    ret.count = sequences.size();

    for (const auto &s: sequences) {
        ret.seq_begins.push_back(ret.nums.size());
        for (uint32_t j = 0; j < s.size(); j++) {
            ret.nums.push_back(to_num(s[j]));
        }
        ret.lengths.push_back(s.size());
    }
    return ret;
}

std::vector<std::string> nums_to_sequences(const SequenceNums &seq_nums)
{
    vector<string> result;

    for (uint32_t i = 0; i < seq_nums.count; i++) {
        string s;
        for (uint32_t j = 0; j < seq_nums.lengths[i]; j++) {
            s += to_symbol(seq_nums.seq_begins[i] + j);
        }
        result.push_back(s);
    }
    return result;
}

