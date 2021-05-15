#pragma once
#include <vector>
#include <string>

struct ArgoCudaParams
{
    bool complementary { true };
    double min_motif_score { 100.0 };
    double max_motif_prob_by_chance { 0.3 };
    double min_motif_presence { 0.0 };
    std::string positive_sequences;
    bool use_real_nucl_frequences { false };
    bool use_binom_instead_chi2 { false };
    int motif_to_find { 10000000 };
    int markov_level { 0 }; // 0 - bernulli, 1 - 2-nucleotides, 2 - 3-nucleotides (if use_real_nucl_frequences is true);
    bool use_old_motifs_file { false };
    double chip_seq_percentage { 0.0 };
    std::string contrast_sequences;
    double max_motif_score_contrast { 0.0 };
    double min_motif_chi2 { 0.0 };
};
