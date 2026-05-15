#include "beg_hamiltonian.h"
#include "neighbor_sums.h"

/**
 * @brief Calculate the energy change for a chemical species exchange using BEG Hamiltonian.
 *
 * The BEG (Binary Ternary) model computes energy as:
 *   ΔE = (typeN - typeA) * [
 *     JOTA * (sumLinNN_A - sumLinNN_N) +
 *     KA * (sumCuadNN_A - sumCuadNN_N) * (typeN + typeA) +
 *     ELE * ((sumLinNN_A - sumLinNN_N) * (typeN + typeA) + (sumCuadNN_A - sumCuadNN_N)) +
 *     ... (same for 2nd NN)
 *   ]
 */
double calculateDeltaChemicalEnergy(int specieA, int specieN,
                                    const BEGCoefficients& coeff,
                                    const NeighborSpeciesSums& sums) {
    // Species difference and sum
    double diffTipo = specieN - specieA;
    double sumTipo  = specieN + specieA;

    // Nearest neighbor differences
    double diffLinNN  = sums.linNN_A  - sums.linNN_N;
    double diffCuadNN = sums.cuadNN_A - sums.cuadNN_N;

    // Next nearest neighbor differences
    double diffLinNNN  = sums.linNNN_A  - sums.linNNN_N;
    double diffCuadNNN = sums.cuadNNN_A - sums.cuadNNN_N;

    // BEG Hamiltonian delta computation
    double deltaEQ = diffTipo * (
        // Nearest Neighbors (NN)
        coeff.jota1 * diffLinNN +
        coeff.ka1   * diffCuadNN * sumTipo +
        coeff.ele1  * (diffLinNN * sumTipo + diffCuadNN) +

        // Next Nearest Neighbors (NNN)
        coeff.jota2 * diffLinNNN +
        coeff.ka2   * diffCuadNNN * sumTipo +
        coeff.ele2  * (diffLinNNN * sumTipo + diffCuadNNN));
    
    return deltaEQ;
}
