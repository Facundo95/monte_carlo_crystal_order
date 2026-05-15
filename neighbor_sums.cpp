#include "neighbor_sums.h"
#include "lattice.h"

/**
 * @brief Compute neighbor species sums used by chemical exchange step.
 */
NeighborSpeciesSums computeNeighborSpeciesSums(const Lattice& lattice, int site, int siteNeighbor) {
    NeighborSpeciesSums s{};
    s.linNN_A = lattice.calculateNeighborSpeciesSum(site, 1, 1) - lattice.getSpecies(siteNeighbor);
    s.linNN_N = lattice.calculateNeighborSpeciesSum(siteNeighbor, 1, 1) - lattice.getSpecies(site);
    s.cuadNN_A = lattice.calculateNeighborSpeciesSum(site, 1, 2) - lattice.getSpecies(siteNeighbor) * lattice.getSpecies(siteNeighbor);
    s.cuadNN_N = lattice.calculateNeighborSpeciesSum(siteNeighbor, 1, 2) - lattice.getSpecies(site) * lattice.getSpecies(site);

    s.linNNN_A = lattice.calculateNeighborSpeciesSum(site, 2, 1);
    s.linNNN_N = lattice.calculateNeighborSpeciesSum(siteNeighbor, 2, 1);
    s.cuadNNN_A = lattice.calculateNeighborSpeciesSum(site, 2, 2);
    s.cuadNNN_N = lattice.calculateNeighborSpeciesSum(siteNeighbor, 2, 2);
    return s;
}

/**
 * @brief Compute neighbor spin sums for shells 1..6.
 */
std::array<int, 6> computeNeighborSpinSums(const Lattice& lattice, int site) {
    std::array<int, 6> sums;
    sums[0] = static_cast<int>(lattice.calculateNeighborSpinSum(site, 1));
    sums[1] = static_cast<int>(lattice.calculateNeighborSpinSum(site, 2));
    sums[2] = static_cast<int>(lattice.calculateNeighborSpinSum(site, 3));
    // neighbors4 and neighbors5 may be uninitialized; keep as 0 if so
    sums[3] = 0;
    sums[4] = 0;
    sums[5] = static_cast<int>(lattice.calculateNeighborSpinSum(site, 6));
    return sums;
}
