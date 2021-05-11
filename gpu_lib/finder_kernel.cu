#include "finder_kernel.cuh"

#include "gpu_memory.h"

namespace {
__device__ uint16_t get_occurrence(uint32_t motif_hash,
                                   uint32_t *seq_hashes,
                                   uint32_t sequences_count,
                                   uint32_t *hash_lengths,
                                   uint32_t *hash_begins,
                                   uint32_t threads_per_block)
{
    // Shared memory size bytes: hash length * sizeof (uint32_t)
    // Max sequence length: shared memory size / 4bytes (16-25k)
    extern __shared__ uint32_t sequence[];

    // Max number of sequences: 64k
    uint16_t result = 0;

    for (uint32_t j = 0; j < sequences_count; j++) {
        uint32_t seq_offset = hash_begins[j];
        uint32_t seq_length = hash_lengths[j];

        // Load sequence to fast shared memory
        for (uint32_t i = threadIdx.x; i < seq_length; i += threads_per_block) {
            if (i < seq_length) {
                sequence[i] = seq_hashes[i + seq_offset];
            }
        }

        __syncthreads();
        // Match motif to sequence
        uint8_t found = 0;
        for (uint32_t i = 0; (i < length) && !found; i++) {
            found = (sequence[i] & motif_hash) == sequence[i];
        }
        result += found;
        __syncthreads();
    }
    return result;
}

 __device__ uint32_t motif_index_to_hash(uint32_t motif_index) {
     uint32_t result = 0;
     uint32_t mult = 1;
     for (uint32_t i = MOTIV_LEN - 1; i >= 0; i--) {
         result += (motif_index % ALPH_SIZE + 1) * mult;
         mult *= HASH_BASE;
         motif_index /= ALPH_SIZE;
     }
     return result;
 }

__global__ void motif_finder_kernel_external(
        uint16_t *weights_out,
        uint32_t *seq_hashes,
        uint32_t sequences_count,
        uint32_t *hash_lengths,
        uint32_t *hash_begins,
        uint32_t *mot_hashes,
        uint32_t mots_to_copy,
        uint32_t threads_per_block)
{
    uint32_t index = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t motif_hash = (index <= mots_to_copy) ? mot_hashes[index] : 0;
    weights_out[index] = get_occurrence(motif_hash, seq_hashes, sequences_count, hash_lengths, hash_begins, threads_per_block);
}

__global__ void motif_finder_kernel_internal(
        uint16_t *weights_out,
        uint32_t *seq_hashes,
        uint32_t sequences_count,
        uint32_t *hash_lengths,
        uint32_t *hash_begins,
        uint32_t mots_to_copy,
        uint32_t motif_idx_offset,
        uint32_t threads_per_block)
{
    uint32_t index = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t motif_hash = (index <= mots_to_copy) ? motif_index_to_hash(motif_idx_offset + index) : 0;
    weights_out[index] = get_occurrence(motif_hash, seq_hashes, sequences_count, hash_lengths, hash_begins, threads_per_block);
}
} // namespace


GpuExternalMemory::GpuExternalMemory(const GpuCudaParams &params, const SequenceHashes &sequence_hashes)
{
    allocator = create_memory_allocator(params.unified_memory);

    weights_out = allocate_on_device_and_init<uint16_t>(allocator->get(), params.motif_range_size, nullptr);
    seq_hashes = allocate_on_device_and_init<uint32_t>(allocator->get(), sequence_hashes.hashes.size(), &(sequence_hashes.hashes[0]));
    seq_lengths = allocate_on_device_and_init<uint32_t>(allocator->get(), sequence_hashes.lengths.size(), &(sequence_hashes.lengths[0]));
    seq_begins = allocate_on_device_and_init<uint32_t>(allocator->get(), sequence_hashes.seq_begins.size(), &(sequence_hashes.seq_begins[0]));
    motif_hashes = allocate_on_device_and_init<uint32_t>(allocator->get(), params.motif_range_size, nullptr);
    sequences_count = sequence_hashes.count;
    shared_memory_size = std::max_element(sequence_hashes.lengths.begin(), sequence_hashes.lengths.end()) * sizeof(uint32_t);
    weights_count = params.motif_range_size;
}

GpuExternalMemory::~GpuExternalMemory()
{
    cudaFree(weights_out);
    cudaFree(seq_hashes);
    cudaFree(seq_lengths);
    cudaFree(seq_begins);
    cudaFree(motif_hashes);
    weights_count = 0;
    sequences_count = 0;
}

void motif_finder_gpu_external(
    const std::vector<uint32_t> &motif_hashes,
    const GpuExternalMemory &mem,
    const GpuCudaParams &params,
    std::vector<uint16_t> &out_motif_weights,
    uint32_t motif_offset,
    uint32_t motifs_count,
    int device_id)
{
    cudaSetDevice(device_id);

    mem.allocator->MEMCPY_TO_DEVICE(mem.motif_hashes, &motif_hashes[motif_offset],  motifs_count * sizeof(uint32_t));
    motif_finder_kernel_external<<<params.motif_range_size/params.threads_per_block, params.threads_per_block, mem.shared_memory_size>>>(
        mem.weights_out,
        mem.seq_hashes,
        mem.sequences_count,
        mem.seq_lengths,
        mem.seq_begins,
        mem.motif_hashes,
        motifs_count,
        motif_offset,
        params.threads_per_block);
    cudaDeviceSynchronize();
    mem.allocator->MEMCPY_TO_HOST(&(out_motif_weights[motif_offset]), mem.weights_out, motifs_count * sizeof(uint16_t));
}

// INTERNAL

GpuInternalMemory::GpuInternalMemory(const GpuCudaParams &params, const SequenceHashes &sequence_hashes)
{
    allocator = create_memory_allocator(params.unified_memory);

    weights_out = allocate_on_device_and_init<uint16_t>(allocator->get(), params.motif_range_size, nullptr);
    seq_hashes = allocate_on_device_and_init<uint32_t>(allocator->get(), sequence_hashes.hashes.size(), &(sequence_hashes.hashes[0]));
    seq_lengths = allocate_on_device_and_init<uint32_t>(allocator->get(), sequence_hashes.lengths.size(), &(sequence_hashes.lengths[0]));
    seq_begins = allocate_on_device_and_init<uint32_t>(allocator->get(), sequence_hashes.seq_begins.size(), &(sequence_hashes.seq_begins[0]));
    sequences_count = sequence_hashes.count;
    shared_memory_size = std::max_element(sequence_hashes.lengths.begin(), sequence_hashes.lengths.end()) * sizeof(uint32_t);
    weights_count = params.motif_range_size;
}

GpuInternalMemory::~GpuInternalMemory()
{
    cudaFree(weights_out);
    cudaFree(seq_hashes);
    cudaFree(seq_lengths);
    cudaFree(seq_begins);
    weights_count = 0;
    sequences_count = 0;
}

void motif_finder_gpu_internal(
    const GpuInternalMemory &mem,
    const GpuCudaParams &params,
    std::vector<uint16_t> &out_motif_weights,
    uint32_t motif_idx_offset,
    uint32_t motifs_count,
    int device_id)
{
    cudaSetDevice(device_id);
    motif_finder_kernel_internal<<<params.motif_range_size/params.threads_per_block, params.threads_per_block, mem.shared_memory_size>>>(
        mem.weights_out,
        mem.seq_hashes,
        mem.sequences_count,
        mem.seq_lengths,
        mem.seq_begins,
        motifs_count,
        motif_idx_offset,
        params.threads_per_block);
    cudaDeviceSynchronize();
    mem.allocator->MEMCPY_TO_HOST(&(out_motif_weights[motif_idx_offset]), mem.weights_out, motifs_count * sizeof(uint16_t));
}
