#include "create_stat_model.h"
#include <memory>

#include <bernulli_stat_model.h>
#include <ideal_stat_model.h>
#include <markov_stat_model.h>

using namespace std;

StatModelPtr create_stat_model(bool use_real_nucl_frequences, uint32_t markov_level, const std::vector<std::string> &sequences, bool complementary)
{
    if (!use_real_nucl_frequences) {
        return make_shared<IdealStatModel>(sequences, complementary);
    } else {
        if (markov_level == 0) {
            return make_shared<BernulliStatModel>(sequences, complementary);
        } else if (markov_level < 5) {
            return make_shared<MarkovStatModel>(sequences, complementary, markov_level);
        }
        throw invalid_argument("Only Markov levels 0-4 are supported");
    }
}