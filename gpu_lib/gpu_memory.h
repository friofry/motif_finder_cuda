//
// Created by andrey on 06.03.2021.
//

#ifndef MOTIF_FINDER_GPU_MEMORY_H
#define MOTIF_FINDER_GPU_MEMORY_H
#include <memory>
#include <driver_types.h>

class IGpuMemory {
public:
    virtual cudaError_t MALLOC(void **devPtr, std::size_t size) = 0;
    virtual void *MEMCPY_TO_DEVICE(void *dest, const void *src, std::size_t count) = 0;
    virtual void *MEMCPY_TO_HOST(void *dest, const void *src, std::size_t count) = 0;
    virtual void *MEMSET(void *ptr, int value, std::size_t num) = 0;
    virtual ~IGpuMemory();
};

class UnifiedGpuMemory : public IGpuMemory {
public:
    cudaError_t MALLOC(void **devPtr, std::size_t size) override;
    void *MEMCPY_TO_DEVICE(void *dest, const void *src, std::size_t count) override;
    void *MEMCPY_TO_HOST(void *dest, const void *src, std::size_t count) override;
    void *MEMSET(void *ptr, int value, std::size_t num) override;
};

class NonUnifiedGpuMemory : public IGpuMemory {
public:
    cudaError_t MALLOC(void **devPtr, std::size_t size) override;
    void *MEMCPY_TO_DEVICE(void *dest, const void *src, std::size_t count) override;
    void *MEMCPY_TO_HOST(void *dest, const void *src, std::size_t count) override;
    void *MEMSET(void *ptr, int value, std::size_t num) override;
};

using GpuMemoryPtr = std::shared_ptr<IGpuMemory>;
GpuMemoryPtr create_memory_allocator(bool unified);


#endif //MOTIF_FINDER_GPU_MEMORY_H
