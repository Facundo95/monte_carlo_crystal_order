#include "beg_hamiltonian.h"
#include "neighbor_sums.h"
#include "lattice.h"

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

double calculateTotalChemicalEnergy(const Lattice& lattice, const BEGCoefficients& coeff) {
    double totalChemicalE = 0.0;
    for (int site = 0; site < lattice.totalSites(); ++site) {
        int Si = lattice.getSpecies(site);

        float sumLin1 = lattice.calculateNeighborSpeciesSum(site, 1, 1);
        float sumCuad1 = lattice.calculateNeighborSpeciesSum(site, 1, 2);
        float sumLin2 = lattice.calculateNeighborSpeciesSum(site, 2, 1);
        float sumCuad2 = lattice.calculateNeighborSpeciesSum(site, 2, 2);

        totalChemicalE += Si * (coeff.jota1 * sumLin1 + coeff.ele1 * sumCuad1) +
                          Si * Si * (coeff.ka1 * sumCuad1 + coeff.ele1 * sumLin1) +
                          Si * (coeff.jota2 * sumLin2 + coeff.ele2 * sumCuad2) +
                          Si * Si * (coeff.ka2 * sumCuad2 + coeff.ele2 * sumLin2);
    }
    // Pair interactions were summed twice across sites; apply 1/2 correction here.
    return 0.5 * totalChemicalE;
}
