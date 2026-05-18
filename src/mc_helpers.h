#ifndef MC_HELPERS_H
#define MC_HELPERS_H

#include <cmath>

// Forward declaration
struct BoltzmannDeltaETable;

/**
 * @struct MCStepResults
 * @brief Aggregates counters and accumulators modified during a Monte Carlo step.
 * Stored as references so callers may keep their existing variables and pass
 * them bundled into a single parameter.
 */
struct MCStepResults {
    double& DeltaEAcum;
    int& changesAccepted;
    int& changesAttempted;

    MCStepResults(double& dE, int& acc, int& att)
        : DeltaEAcum(dE), changesAccepted(acc), changesAttempted(att) {}
};

/**
 * @brief Metropolis acceptance test. Returns true if the change should be accepted.
 * @param dE Energy change
 * @param table Boltzmann factor lookup table
 * @return true if change should be accepted, false otherwise
 */
bool metropolisAccept(double dE, BoltzmannDeltaETable& table);

#endif // MC_HELPERS_H
