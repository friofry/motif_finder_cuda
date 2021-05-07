#ifndef COMMON_CONFIG_H
#define COMMON_CONFIG_H

// Motif Length
#define MOTIV_LEN 8

// IUPAC alphabet size
#define ALPH_SIZE 15

// Aligned value for a single letter in hash
#define HASH_BASE 16

// Total number of motifs (len = 8, 15^8)
#define TOTAL_MOT 2562890625u

// motif portion size per kernel run
const unsigned int MOT_PER_CHUNK = 512*64*64*4;

// CUDA unified memory
#define UNIFIED_MEMORY_SUPPORT false

// CUDA threads per block
#define THREADS_PER_BLOCK (512)

#endif // COMMON_CONFIG_H