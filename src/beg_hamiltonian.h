#ifndef BEG_HAMILTONIAN_H
#define BEG_HAMILTONIAN_H

#include <cmath>

// Forward declaration
class Lattice;
struct NeighborSpeciesSums;

/**
 * @struct BEGCoefficients
 * @brief Pre-computed BEG (Binary Ternary) Hamiltonian coefficients.
 * These parameters define the interaction energy in the BEG model.
 */
struct BEGCoefficients {
    double jota1;  ///< Coefficient for w1_13 (1st NN interaction species 1-3)
    double jota2;  ///< Coefficient for w2_13 (2nd NN interaction species 1-3)
    double ka1;    ///< Interaction coefficient combining w1_12, w1_23, w1_13 (1st NN)
    double ka2;    ///< Interaction coefficient combining w2_12, w2_23, w2_13 (2nd NN)
    double ele1;   ///< Electrostatic-like term from w1_12 - w1_23 (1st NN)
    double ele2;   ///< Electrostatic-like term from w2_12 - w2_23 (2nd NN)
};

/**
 * @brief Calculate the change in chemical exchange energy (BEG Hamiltonian).
 * 
 * Computes the energy difference if two species at neighboring sites exchange places,
 * using the BEG model with first and second nearest neighbor interactions.
 *
 * @param specieA Species type at site A
 * @param specieN Species type at neighbor site N
 * @param coeff BEG Hamiltonian coefficients
 * @param sums Pre-computed neighbor species sums for both sites
 * @return Energy change from the exchange
 */
double calculateDeltaChemicalEnergy(int specieA, int specieN,
                                    const BEGCoefficients& coeff,
                                    const NeighborSpeciesSums& sums);

/**
 * @brief Calculate the total chemical (BEG) energy for the entire lattice.
 *
 * This returns the BEG chemical contribution already divided by two where
 * appropriate (to correct for double-counting of pair interactions), so the
 * returned value can be summed directly with the Ising magnetic contribution.
 */
double calculateTotalChemicalEnergy(const Lattice& lattice, const BEGCoefficients& coeff);

#endif // BEG_HAMILTONIAN_H
