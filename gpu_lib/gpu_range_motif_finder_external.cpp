#include "gpu_range_motif_finder_external.h"

#include <algorithm>

#include "finder_kernel.cuh"

GpuRangeMotifFinderExternal::GpuRangeMotifFinderExternal(const SequenceHashesParams &params, const GpuCudaParams &gpu_params, const std::vector<uint32_t> &_motif_hashes)
    : _params(params)
    , _gpu_params(gpu_params)
    , _motif_hashes(motif_hashes)
{
    prepare_memory();
    _shared_memory_size = std::max_element(params._hash_lengths.begin(), params._hash_lengths.end()) * sizeof(uint32_t);
}

GpuRangeMotifFinderExternal::~GpuRangeMotifFinderExternal()
{
    free_memory();
}

void GpuRangeMotifFinderExternal::get_occurrence_in_range(uint32_t begin, uint32_t end, std::vector<uint16_t> &out_motif_weights) const
{
    auto mots_to_copy = end - begin;

    mem->MEMCPY_TO_DEVICE(_motif_hashes_dev, &_motif_hashes[begin], mots_to_copy * sizeof(uint32_t));
    motif_finder_gpu_external(_weights_out,
                              _seq_hashes_dev,
                              _params.count,
                              _seq_lengths_dev,
                              _seq_begins_dev,
                              _motif_hashes_dev,
                              mots_to_copy,
                              _gpu_params.motif_range_size
                              _shared_memory_size,
                              begin,
                              _gpu_params.threads_per_block);
    _mem->MEMCPY_TO_HOST(&(out_motif_weights[0]), _weights_out, mots_to_copy * sizeof(uint16_t));
}

void GpuRangeMotifFinderExternal::prepare_memory()
{
    _mem = create_memory_allocator(_gpu_params.unified_memory);

    _weights_dev = allocate_on_device_and_init<uint16_t>(mem.get(), _gpu_params.motif_range_size, nullptr);
    _seq_hashes_dev = allocate_on_device_and_init<uint32_t>(mem.get(), _params.hashes.size(), &(_params.hashes[0]));
    _seq_lengths_dev = allocate_on_device_and_init<uint32_t>(mem.get(), _params.lengths.size(), &(_params.lengths[0]));
    _seq_begins_dev = allocate_on_device_and_init<uint32_t>(mem.get(), _params.seq_begins.size(), &(_params.seq_begins[0]));
    _motif_hashes_dev =  allocate_on_device_and_init<uint32_t>(mem.get(), _gpu_params.motif_range_size, nullptr);
}

void GpuRangeMotifFinderExternal::free_memory()
{
    cudaFree(_weights_dev);
    cudaFree(_seq_hashes_dev);
    cudaFree(_seq_lengths_dev);
    cudaFree(_seq_begins_dev);
    cudaFree(_motif_hashes_dev);
}
