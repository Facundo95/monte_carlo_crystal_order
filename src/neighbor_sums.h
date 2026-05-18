#ifndef NEIGHBOR_SUMS_H
#define NEIGHBOR_SUMS_H

#include <array>

// Forward declaration
class Lattice;

/**
 * @struct NeighborSpeciesSums
 * @brief Holds precomputed neighbor species sums for a site and its chosen neighbor.
 */
struct NeighborSpeciesSums {
    int linNN_A;
    int linNN_N;
    int cuadNN_A;
    int cuadNN_N;
    int linNNN_A;
    int linNNN_N;
    int cuadNNN_A;
    int cuadNNN_N;
};

/**
 * @brief Compute and return neighbor species sums needed for chemical energy.
 * @param lattice The lattice
 * @param site The active site index
 * @param siteNeighbor The neighbor site index being exchanged with
 * @return Struct containing all neighbor sum components
 */
NeighborSpeciesSums computeNeighborSpeciesSums(const Lattice& lattice, int site, int siteNeighbor);

/**
 * @brief Compute and return spin neighbor sums for shells 1..6 for a given site.
 * @param lattice The lattice
 * @param site The site index
 * @return Array of 6 neighbor spin sums
 */
std::array<int, 6> computeNeighborSpinSums(const Lattice& lattice, int site);

#endif // NEIGHBOR_SUMS_H
