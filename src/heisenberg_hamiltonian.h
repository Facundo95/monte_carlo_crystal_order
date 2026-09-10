#ifndef HEISENBERG_HAMILTONIAN_H
#define HEISENBERG_HAMILTONIAN_H

#include <array>

using HeisenbergVector = std::array<double, 3>;
using HeisenbergNeighborSums = std::array<HeisenbergVector, 6>;

/**
 * @struct HeisenbergCouplings
 * @brief Heisenberg exchange couplings for up to six neighbor shells.
 */
struct HeisenbergCouplings {
    double Jm1;
    double Jm2;
    double Jm3;
    double Jm4;
    double Jm5;
    double Jm6;

    /** @brief Extract couplings in shell order. */
    std::array<double, 6> toArray() const {
        return {Jm1, Jm2, Jm3, Jm4, Jm5, Jm6};
    }
};

/** @brief Return the dot product of two three-dimensional moments. */
double dotProduct(const HeisenbergVector& first, const HeisenbergVector& second);

/**
 * @brief Calculate the energy change for replacing one moment.
 *
 * The local Hamiltonian is E = -sum(J_n Si dot Sj) - H Siz.
 * Neighbor sums contain the vector sum for each shell around the site.
 */
double calculateDeltaHeisenbergEnergy(const HeisenbergVector& currentMoment,
                                      const HeisenbergVector& proposedMoment,
                                      double externalField,
                                      const HeisenbergCouplings& couplings,
                                      const HeisenbergNeighborSums& neighborSums);

#endif // HEISENBERG_HAMILTONIAN_H
