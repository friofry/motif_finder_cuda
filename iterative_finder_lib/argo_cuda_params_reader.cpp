#include "argo_cuda_params_reader.h"
#include <fstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <regex>

using namespace std;

namespace
{
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
        cout << "Ini file should contain exactly 13 parameters. " << inifile << lines.size() << endl;
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

    return result;
}

void print_argo_cuda_params(const ArgoCudaParams &params)
{
    cout << (int)params.complementary << "\t\tComplementarity" << endl;
    cout << params.min_motif_score << "\t\tMin score value" << endl;
    cout << params.max_motif_prob_by_chance << "\t\tMax presence of motif by chance" << endl;
    cout << params.min_motif_presence << "\t\tMin presence of motif in the seq set" << endl;
    cout << params.positive_sequences << "\tSequences set" << endl;
    cout << (int)params.use_real_nucl_frequences << "\t\tUsing of real frequencies. 0-not/1-use" << endl;
    cout << (int)params.use_binom_instead_chi2 << "\t\tusing of the Chi2 (0) or binom (1) score" << endl;
    cout << params.motif_to_find  << "\t\tNumber of motifs. 0->1000000" << endl;
    cout << params.markov_level << "\t\tMarkov level of frequencies, if RealFreqs==1. 0- bernulli/1- dinucleotides/2- trinucleotides." << endl;
    cout << (int)params.use_old_motifs_file << "\t\tUse (1) or not (0) of earlier found motifs from oldmotifs.mot file." << endl;
    cout << params.chip_seq_percentage << "\t\tUse chip seq (minimum percentage of sequences)" << endl;
    cout << params.contrast_sequences << "\tContrast sequences set" << endl;
    cout << params.max_motif_score_contrast << "\t\tMax score value in contrast seq set" << endl;
    cout << params.min_motif_chi2 << "\t\tMinimum motif chi2 score" << endl;
}