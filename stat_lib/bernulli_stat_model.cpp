#include "bernulli_stat_model.h"
#include "probability.h"

using namespace std;

BernulliStatModel::BernulliStatModel(const std::vector<std::string> &sequences,
                                     bool complementary,
                                     bool use_binom_instead_of_chi2)
    : StatModel(sequences, complementary, use_binom_instead_of_chi2)
{
    // TODO: определиться как считать вероятность для комплементарных последовательностей.
    _frequencies_ratio = calc_frequencies_ratio(sequences, false);
    _probability_x4 = ::calc_probability_x4(_frequencies_ratio);
}

double BernulliStatModel::motif_probability(uint32_t hash) const
{
    return ::motif_probability(hash, _frequencies_ratio);
}

double BernulliStatModel::motif_probability_x4(uint32_t hash) const
{
    return ::motif_probability_x4(hash, _probability_x4);
}
