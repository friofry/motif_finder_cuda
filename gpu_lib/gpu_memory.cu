//
// Created by andrey on 06.03.2021.
//

#include "gpu_memory.h"
#include <cstring>

#include <cuda_runtime_api.h>
#include <cuda.h>

IGpuMemory::~IGpuMemory() {

}

// unified memory
cudaError_t UnifiedGpuMemory::MALLOC(void **devPtr, std::size_t size) {
    return cudaMallocManaged(devPtr, size);
}
void *UnifiedGpuMemory::MEMCPY_TO_DEVICE(void *dest, const void *src, std::size_t count) {
    return memcpy((wchar_t *)dest, src, count);
}

void *UnifiedGpuMemory::MEMCPY_TO_HOST(void *dest, const void *src, std::size_t count) {
    return memcpy((wchar_t *)dest, src, count);
}

void *UnifiedGpuMemory::MEMSET(void *ptr, int value, std::size_t num) {
    return memset((wchar_t *)ptr, value, num);
}

// non unified memory
cudaError_t NonUnifiedGpuMemory::MALLOC(void **devPtr, std::size_t size) {
    return cudaMalloc(devPtr, size);
}

void *NonUnifiedGpuMemory::MEMCPY_TO_DEVICE(void *dest, const void *src, std::size_t count) {
    cudaMemcpy(dest, src, count, cudaMemcpyHostToDevice);
    return dest;
}

void *NonUnifiedGpuMemory::MEMCPY_TO_HOST(void *dest, const void *src, std::size_t count) {
    cudaMemcpy(dest, src, count, cudaMemcpyDeviceToHost);
    return dest;
}

void *NonUnifiedGpuMemory::MEMSET(void *ptr, int value, std::size_t num) {
    cudaMemset(ptr, value, num);
    return ptr;
}

// factory
GpuMemoryPtr create_memory_allocator(bool unified) {
    if (unified) {
        return GpuMemoryPtr(new UnifiedGpuMemory());
    } else {
        return GpuMemoryPtr(new NonUnifiedGpuMemory());
    }
}
