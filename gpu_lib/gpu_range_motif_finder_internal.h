#ifndef MOTIF_FINDER_GPURANGEFINDER_INTERNAL_H
#define MOTIF_FINDER_GPURANGEFINDER_INTERNAL_H

#include <sequence_hashes_params.h>
#include <gpu_cuda_params.h>
#include <gpu_memory.h>
#include <vector>

class GpuRangeMotifFinderInternal {
public:
    GpuRangeMotifFinderInternal(const SequenceHashesParams &params, const GpuCudaParams &gpu_params);
    virtual ~GpuRangeMotifFinderInternal();

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
    GpuMemoryPtr _mem;
    uint32_t _shared_memory_size {0};
    uint16_t *_weights_out {nullptr};
    uint32_t *_seq_hashes {nullptr};
    uint32_t *_seq_lengths {nullptr};
    uint32_t *_seq_begins {nullptr};
};

#endif // MOTIF_FINDER_GPURANGEFINDER_INTERNAL_H
