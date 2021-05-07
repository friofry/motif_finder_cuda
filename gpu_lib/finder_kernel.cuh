#ifndef MOTIF_FINDER_FINDER_KERNEL_CUH
#define MOTIF_FINDER_FINDER_KERNEL_CUH

// these are wrapper functions for cuda kernels

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
        uint32_t threads_per_block);

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
        uint32_t threads_per_block);

#endif //MOTIF_FINDER_FINDER_KERNEL_CUH
