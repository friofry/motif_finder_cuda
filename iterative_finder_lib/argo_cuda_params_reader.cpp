#include "argo_cuda_params_reader.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

using namespace std;

namespace {
vector<string> readlines(const char *fillename)
{
    vector<string> result;
    ifstream input(fillename);

    if (input.eof() || input.bad()) {
        return result;
    }

    string line;

    while (getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        result.push_back(line);
    }
    return result;
}

string get_first_word(const string &s)
{
    regex rgx("[^\\s]+");
    sregex_iterator it(s.begin(), s.end(), rgx);
    if (it != sregex_iterator()) {
        smatch match = *it;
        return match.str();
    }
    return "";
}

} // namespace

ArgoCudaParams read_ini_file(const char *inifile)
{
    auto lines = readlines(inifile);
    ArgoCudaParams result;
    if (lines.size() < 13) {
        cout << "Ini file should contain exactly 15 parameters. " << inifile << lines.size() << endl;
        return result;
    }

    result.complementary = atoi(lines[0].c_str());
    result.min_motif_score = stod(lines[1]);
    result.max_motif_prob_by_chance = stod(lines[2]) / 100.0;
    result.min_motif_presence = stod(lines[3]) / 100.0;
    result.positive_sequences = get_first_word(lines[4]);
    result.use_real_nucl_frequences = atoi(lines[5].c_str());
    result.use_binom_instead_chi2 = atoi(lines[6].c_str());
    result.motif_to_find = atoi(lines[7].c_str());
    if (result.motif_to_find == 0) {
        result.motif_to_find = 1000000;
    }
    result.markov_level = atoi(lines[8].c_str());
    result.use_old_motifs_file = atoi(lines[9].c_str());
    result.chip_seq_percentage = stod(lines[10]) / 100.0;
    result.contrast_sequences = get_first_word(lines[11]);
    result.max_motif_score_contrast = stod(lines[12]);
    result.min_motif_chi2 = stod(lines[13]);
    result.bonferroni_correction = atoi(lines[14].c_str());
    return result;
}

void print_argo_cuda_params(const ArgoCudaParams &params)
{
    cout << (int)params.complementary << "\t\tComplementarity. 0 - forward strand. 1 - forward + reverse strand."
         << endl;
    cout << params.min_motif_score << "\t\tMinimum score" << endl;
    cout << params.max_motif_prob_by_chance * 100
         << "\t\tMaximum presence of motif for random reasons in the positive set of sequences [0-100]" << endl;
    cout << params.min_motif_presence * 100
         << "\t\tMinimum presence of motif for random reasons in a positive set of sequences [0-100]" << endl;
    cout << params.positive_sequences << "\tFile with positive set of sequences" << endl;
    cout << (int)params.use_real_nucl_frequences
         << "\t\t0 - neutral frequencies, 1 - real nucleotide frequencies in the set of sequences [0, 1]" << endl;
    cout << (int)params.use_binom_instead_chi2 << "\t\tScore type: 0 - chi-squared or 1 - binomial [0, 1]" << endl;
    cout << params.motif_to_find << "\t\tMaximum number of result motifs. 0 - reveal all significant motifs [0, ]"
         << endl;
    cout << params.markov_level
         << "\t\tMarkov chain order (0-Bernulli, 1-dinucleotide, 2-trinucleotide), when using real nucleotide "
            "frequencies [0-3]"
         << endl;
    cout << (int)params.use_old_motifs_file << "\t\t[deprecated]" << endl;
    cout << params.chip_seq_percentage * 100 << "\t\t[deprecated]" << endl;
    cout << params.contrast_sequences << "\tFile with contrast set of sequences" << endl;
    cout << params.max_motif_score_contrast << "\t\tMaximum score in a contrast set of sequences" << endl;
    cout << params.min_motif_chi2 << "\t\t[deprecated]" << endl;
    cout << params.bonferroni_correction << "\t\tOutput results with Bonferroni correction [0, 1]" << endl;
}