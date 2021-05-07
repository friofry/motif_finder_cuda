#include "gpu_range_motif_finder_internal.h"

#include <algorithm>

#include "finder_kernel.cuh"

GpuRangeMotifFinderInternal::GpuRangeMotifFinderInternal(const SequenceHashesParams &params, const GpuCudaParams &gpu_params)
    : _params(params)
    , _gpu_parmas(gpu_params)
{
    prepare_memory();
    _shared_memory_size = std::max_element(params._hash_lengths.begin(), params._hash_lengths.end()) * sizeof(uint32_t);
}

GpuRangeMotifFinderInternal::~GpuRangeMotifFinderInternal()
{
    free_memory();
}

void GpuRangeMotifFinderInternal::get_occurrence_in_range(uint32_t begin, uint32_t end, std::vector<uint16_t> &out_motif_weights) const
{
    auto mots_to_copy = end - begin;
    motif_finder_gpu_internal(_weights_out,
                              _seq_hashes,
                              _params.count,
                              _params._hash_lengths,
                              _params._hash_begins,
                              mots_to_copy,
                              _gpu_params.motif_range_size
                              _shared_memory_size,
                              begin,
                              _gpu_params.threads_per_block);

    _mem->MEMCPY_TO_HOST(&(out_motif_weights[0]), _weights_out, mots_to_copy * sizeof(uint16_t));
}

void GpuRangeMotifFinderInternal::prepare_memory()
{
    _mem = create_memory_allocator(_gpu_params.unified_memory);

    _weights_out = allocate_on_device_and_init<uint16_t>(mem.get(), _gpu_params.motif_range_size, nullptr);
    _seq_hashes = allocate_on_device_and_init<uint32_t>(mem.get(), _params.hashes.size(), &(_params.hashes[0]));
    _seq_lengths = allocate_on_device_and_init<uint32_t>(mem.get(), _params.lengths.size(), &(_params.lengths[0]));
    _seq_begins = allocate_on_device_and_init<uint32_t>(mem.get(), _params.seq_begins.size(), &(_params.seq_begins[0]));
}

void GpuRangeMotifFinderInternal::free_memory()
{
    cudaFree(_weights_out);
    cudaFree(_seq_hashes);
    cudaFree(_hash_lengths);
    cudaFree(_hash_begins);
}
