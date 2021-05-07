#ifndef MOTIF_FINDER_GPURANGEFINDER_EXTERNAL_H
#define MOTIF_FINDER_GPURANGEFINDER_EXTERNAL_H

#include <sequence_hashes_params.h>
#include <gpu_cuda_params.h>
#include <gpu_memory.h>
#include <vector>

class GpuRangeMotifFinderExternal {
public:
    GpuRangeMotifFinderExternal(const SequenceHashesParams &params, const GpuCudaParams &gpu_params, const std::vector<uint32_t> &motif_hashes);
    virtual ~GpuRangeMotifFinderExternal();

    // Get motif occurrences and store to out_motif_weights
    // @params begin - first motif index
    // @params end - last motif index + 1 (points after the last motif)
    // @params out_motif_weights - array for results (of at least (end - begin) size)
    void get_occurrence_in_range(uint32_t begin, uint32_t end, std::vector<uint16_t> &out_motif_weights) const;

private:
    void prepare_memory();
    void free_memory();

private:
    const SequenceHashesParams &_params;
    const GpuCudaParams &_gpu_params;
    const std::vector<uint32_t> &_motif_hashes;
    GpuMemoryPtr _mem;
    uint32_t _shared_memory_size {0};
    uint16_t *_weights_dev {nullptr};
    uint32_t *_seq_hashes_dev {nullptr};
    uint32_t *_seq_lengths_dev {nullptr};
    uint32_t *_seq_begins_dev {nullptr};
    uint32_t *_motif_hashes_dev {nullptr};
};

#endif // MOTIF_FINDER_GPURANGEFINDER_EXTERNAL_H
