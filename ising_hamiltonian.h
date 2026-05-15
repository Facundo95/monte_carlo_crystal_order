#ifndef ISING_HAMILTONIAN_H
#define ISING_HAMILTONIAN_H

#include <array>
#include <cmath>

/**
 * @struct IsingCouplings
 * @brief Ising Hamiltonian magnetic couplings for up to 6 neighbor shells.
 */
struct IsingCouplings {
    double Jm1; ///< 1st NN magnetic interaction energy
    double Jm2; ///< 2nd NN magnetic interaction energy
    double Jm3; ///< 3rd NN magnetic interaction energy
    double Jm4; ///< 4th NN magnetic interaction energy
    double Jm5; ///< 5th NN magnetic interaction energy
    double Jm6; ///< 6th NN magnetic interaction energy
    
    /**
     * @brief Extract couplings into array for iteration.
     * @return Array of 6 couplings: {Jm1, Jm2, Jm3, Jm4, Jm5, Jm6}
     */
    std::array<double, 6> toArray() const {
        return {Jm1, Jm2, Jm3, Jm4, Jm5, Jm6};
    }
};

/**
 * @brief Calculate the energy change for a spin flip in an Ising model with external field.
 *
 * Ising Hamiltonian: H = -H*m - ΣJm*m*mn
 * Energy change for spin flip at site:
 *   ΔE = 2*m*(H + ΣJm*sumNeighbors)
 *
 * @param spinAtSite Current spin value (±1 or 0)
 * @param externalField External magnetic field H
 * @param couplings Magnetic couplings for 6 neighbor shells
 * @param neighborSpinSums Pre-computed spin sums for each shell (array of 6)
 * @return Energy change from flipping the spin
 */
double calculateDeltaIsingEnergy(int spinAtSite,
                                 double externalField,
                                 const IsingCouplings& couplings,
                                 const std::array<int, 6>& neighborSpinSums);

#endif // ISING_HAMILTONIAN_H
