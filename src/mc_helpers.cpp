#include "mc_helpers.h"
#include "simulation.h"
#include "rng.h"

/**
 * @brief Metropolis acceptance test. Returns true if the change should be accepted.
 */
bool metropolisAccept(double dE, BoltzmannDeltaETable& table) {
    if (dE <= 0.0) return true;
    double epsilon = Ran0a1();
    double boltz = table.lookup(dE);
    return boltz >= epsilon;
}
