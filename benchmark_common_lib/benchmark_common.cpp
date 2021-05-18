#include "benchmark_common.h"

#include <cstdint>
#include <fstream>
#include <memory>

#include <bernulli_stat_model.h>
#include <ideal_stat_model.h>
#include <markov_stat_model.h>

using namespace std;

void save_to_file(const std::vector<std::string> &sequences, const char *filename)
{
    ofstream fs(filename);
    for (uint32_t i = 0; i < sequences.size(); i++) {
        fs << ">" << i << endl;
        fs << sequences[i] << endl;
    }
    fs << endl;
}

StatModelPtr create_stat_model(const std::string &model, const std::vector<std::string> &sequences, bool complementary)
{
    if (model == "bernulli") {
        return make_shared<BernulliStatModel>(sequences, complementary, false);
    } else if (model == "markov1") {
        return make_shared<MarkovStatModel>(sequences, complementary, 1, false);
    } else if (model == "markov2") {
        return make_shared<MarkovStatModel>(sequences, complementary, 2, false);
    } else if (model == "markov3") {
        return make_shared<MarkovStatModel>(sequences, complementary, 3, false);
    } else {
        // "ideal"
        return std::make_shared<IdealStatModel>(sequences, complementary, false);
    }
}