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

class GpuMemoryPointer {
public:
    GpuMemoryPointer();


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

// Create a new device memory, fill it from host src memory, or with zeroes
// @param mem - memory allocator
// @param count - number of elements
// @param src - if null, fill with zeroes, otherwise copy from the host memory to device memory
template<typename T>
T *allocate_on_device_and_init(IGpuMemory *mem, uint32_t count, T *src = nullptr)
{
    T *result_dev;
    uint32_t result_dev_bytes = count * sizeof (T);
    mem->MALLOC((void**)&result_dev, result_dev_bytes);
    if (src) {
        mem->MEMCPY_TO_DEVICE(result_dev, src, result_dev_bytes);
    } else {
        mem->MEMSET(result_dev, 0, result_dev_bytes);
    }
    return result_dev;
}

#endif //MOTIF_FINDER_GPU_MEMORY_H
