#include "ideal_stat_model.h"

#include <probability.h>

using namespace std;

IdealStatModel::IdealStatModel(const std::vector<std::string> &sequences,
                               bool complementary,
                               bool use_binom_instead_of_chi2)
    : StatModel(sequences, complementary, use_binom_instead_of_chi2)
    , _frequencies_ratio({ 0.25, 0.25, 0.25, 0.25 })
{
    _probability_x4 = ::calc_probability_x4(_frequencies_ratio);
}

double IdealStatModel::motif_probability(uint32_t hash) const
{
    return ::motif_probability(hash, _frequencies_ratio);
}

double IdealStatModel::motif_probability_x4(uint32_t hash) const
{
    return ::motif_probability_x4(hash, _probability_x4);
}
