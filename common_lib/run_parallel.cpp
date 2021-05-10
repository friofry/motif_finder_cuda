#include "run_parallel.h"

#include <thread>
#include <vector>

using namespace std;

void run_parallel(uint32_t thread_count, ThreadFunction func)
{
    vector<thread> threads;
    for (uint32_t i = 0; i < thread_count; i++) {
        threads.push_back(std::thread(func, i));
    }

    // Wait for finished
    for (uint32_t i = 0; i < thread_count; i++) {
        threads[i].join();
    }
}