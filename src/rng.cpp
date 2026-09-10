#include "rng.h"
#include <cmath>
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

HeisenbergVector randomUnitVector() {
    constexpr double pi = 3.14159265358979323846;
    const double z = 2.0 * static_cast<double>(Ran0a1()) - 1.0;
    const double azimuth = 2.0 * pi * static_cast<double>(Ran0a1());
    const double radialComponent = std::sqrt(std::max(0.0, 1.0 - z * z));
    return {radialComponent * std::cos(azimuth),
            radialComponent * std::sin(azimuth),
            z};
}

HeisenbergVector randomGaussianVector() {
    constexpr double twoPi = 6.28318530717958647692;
    auto gaussianPair = [twoPi]() {
        const double firstUniform = std::max(static_cast<double>(Ran0a1()), 1e-12);
        const double secondUniform = static_cast<double>(Ran0a1());
        const double radius = std::sqrt(-2.0 * std::log(firstUniform));
        const double angle = twoPi * secondUniform;
        return std::array<double, 2>{radius * std::cos(angle),
                                     radius * std::sin(angle)};
    };

    const std::array<double, 2> firstPair = gaussianPair();
    const std::array<double, 2> secondPair = gaussianPair();
    return {firstPair[0], firstPair[1], secondPair[0]};
}

HeisenbergVector gaussianSpinProposal(const HeisenbergVector& moment,
                                      double sigma) {
    const HeisenbergVector gaussian = randomGaussianVector();
    HeisenbergVector proposal{
        moment[0] + sigma * gaussian[0],
        moment[1] + sigma * gaussian[1],
        moment[2] + sigma * gaussian[2]
    };
    const double norm = std::sqrt(proposal[0] * proposal[0] +
                                  proposal[1] * proposal[1] +
                                  proposal[2] * proposal[2]);
    if (norm == 0.0) {
        return moment;
    }

    for (double& component : proposal) {
        component /= norm;
    }
    return proposal;
}
