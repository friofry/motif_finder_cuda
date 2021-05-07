#include "gpu_info.h"

#include <stdio.h>
#include <cuda_runtime_api.h>
#include <cuda.h>

int gpu_count()
{
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    //printf("Device Count: %d\n", nDevices);
    for (int i = 0; i < nDevices; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        if (false) {
            printf("Device Number: %d\n", i);
            printf("  Device name: %s\n", prop.name);
            printf("  Memory Clock Rate (KHz): %d\n",
                   prop.memoryClockRate);
            printf("  Memory Bus Width (bits): %d\n",
                   prop.memoryBusWidth);
            printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
                   2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
        }
    }
    return nDevices;
}