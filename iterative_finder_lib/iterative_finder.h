#pragma once

#include <string>
#include <vector>

#include <stat_model.h>

#include "argo_cuda_params.h"
#include "motif_data.h"
#include "sequence_hashes.h"

using FindOccurrencesAlgorithm = std::function<void(const std::vector<uint32_t> &motif_hashes,
                                                    const SequenceHashes &sequence_hashes,
                                                    std::vector<uint16_t> &out_motif_weights)>;

class ImportantMotifFinder {
public:
    ImportantMotifFinder(const ArgoCudaParams &params,
                         const FindOccurrencesAlgorithm &external_algorithm,
                         const char *output_file = "a.txt");
    std::vector<uint32_t> find();

private:
    void find_motifs_iterative(std::vector<uint32_t> &motif_hashes, bool exclude_sequence = false);
    void exclude_motifs_by_score(std::vector<uint32_t> &motif_hashes, std::vector<uint16_t> &weights) const;
    void write_results_old();

private:
    ArgoCudaParams _params;
    const FindOccurrencesAlgorithm &_external_algorithm;
    StatModelPtr _stat_model;
    std::vector<uint32_t> _found_motifs;
    std::vector<MotifData> _found_motifs_data;
    std::string _output_file;
    SequenceHashes _sequence_hashes;
};

std::vector<uint32_t> find_important_motifs(const ArgoCudaParams &params,
                                            const FindOccurrencesAlgorithm &external_algorithm,
                                            const char *output_file = "a.txt");