
#ifndef MOTIF_FINDER_RUN_PARALLEL_H
#define MOTIF_FINDER_RUN_PARALLEL_H

#include <functional>
#include <cstdint>

using ThreadFunction = std::function<void(uint32_t thread_id)>;

void run_parallel(uint32_t threads, ThreadFunction func);

#endif //MOTIF_FINDER_RUN_PARALLEL_H
