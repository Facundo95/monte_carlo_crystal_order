#include "rng.h"
#include <random>
#include <mutex>

// Thread-safe PRNG using std::mt19937_64. Seed is fixed for deterministic runs;
// change the seed value below to get different sequences.
static std::mt19937_64& global_engine() {
    static std::mt19937_64 eng(5489ULL);
    return eng;
}

static std::mutex& global_engine_mutex() {
    static std::mutex m;
    return m;
}

long RanEnt(){
    std::lock_guard<std::mutex> lock(global_engine_mutex());
    return static_cast<long>(global_engine()());
}

float Ran0a1(){
    std::lock_guard<std::mutex> lock(global_engine_mutex());
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    return dist(global_engine());
}

int RanEnt1a8(){
    std::lock_guard<std::mutex> lock(global_engine_mutex());
    std::uniform_int_distribution<int> dist(1, 8);
    return dist(global_engine());
}
