#include "finder_kernel.cuh"

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

void motif_finder_gpu_internal(
        uint16_t *weights_out,
        uint32_t *seq_hashes,
        uint32_t sequences_count,
        uint32_t *hash_lengths,
        uint32_t *hash_begins,
        uint32_t mots_to_copy,
        uint32_t motif_range_size,
        uint32_t shared_memory_size,
        uint32_t motif_idx_offset,
        uint32_t threads_per_block)
{
    motif_finder_kernel_internal<<<motif_range_size/threads_per_block, threads_per_block, shared_memory_size>>>(
            weights_out,
            seq_hashes,
            sequences_count,
            hash_lengths,
            hash_begins,
            mots_to_copy,
            motif_idx_offset,
            threads_per_block);
    cudaDeviceSynchronize();
}

void motif_finder_gpu_external(
        uint16_t *weights_out,
        uint32_t *seq_hashes,
        uint32_t sequences_count,
        uint32_t *hash_lengths,
        uint32_t *hash_begins,
        uint32_t *motif_hashes,
        uint32_t mots_to_copy,
        uint32_t motif_range_size,
        uint32_t shared_memory_size,
        uint32_t motif_idx_offset,
        uint32_t threads_per_block)
{
    motif_finder_kernel_external<<<motif_range_size/threads_per_block, threads_per_block, shared_memory_size>>>(
            weights_out,
            seq_hashes,
            sequences_count,
            hash_lengths,
            hash_begins,
            motif_hashes,
            mots_to_copy,
            motif_idx_offset,
            threads_per_block);
    cudaDeviceSynchronize();
}

