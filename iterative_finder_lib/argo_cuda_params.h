#pragma once

#include <string>

// clang-format off
struct ArgoCudaParams {
    bool complementary{ true };                 // Complementarity. 0 - forward strand. 1 - forward + reverse strand.
    double min_motif_score{ 100.0 };            // Minimum score
    double max_motif_prob_by_chance{ 0.3 };     // Maximum presence of motif for random reasons in the positive set of sequences [0-1]
    double min_motif_presence{ 0.0 };           // Minimum presence of motif for random reasons in a positive set of sequences [0-1]
    std::string positive_sequences;             // File with positive set of sequences
    bool use_real_nucl_frequences{ false };     // 0 - neutral frequencies, 1 - real nucleotide frequencies in the set of sequences [0, 1]
    bool use_binom_instead_chi2{ true };       // [deprecated] Score type: 0 - chi-squared or 1 - binomial [0, 1]
    int motif_to_find{ 10000000 };              // Maximum number of result motifs. 0 - reveal all significant motifs [0, ]
    int markov_level{ 0 };                      // Markov chain order (0-Bernulli, 1-dinucleotide, 2-trinucleotide), when using real nucleotide frequencies [0-3]
    bool use_old_motifs_file{ false };          // [deprecated]
    double chip_seq_percentage{ 0.0 };          // [deprecated]
    std::string contrast_sequences;             // File with contrast set of sequences
    double max_motif_score_contrast{ 0.0 };     // Maximum score in a contrast set of sequences
    double min_motif_chi2{ 0.0 };               // [deprecated]
    bool bonferroni_correction{ false };        // Output results with Bonferroni correction [0, 1]
    std::string output_file;                    // Output file
    bool int_results { true };                  // Output integer scores
};
// clang-format on