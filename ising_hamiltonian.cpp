#include "ising_hamiltonian.h"

/**
 * @brief Calculate the energy change for a spin flip in Ising model with external field.
 *
 * For a spin flip at a site: ΔE = 2*m*(H + Σ Jm*neighborSpins)
 * This is twice the effective field times the spin value.
 */
double calculateDeltaIsingEnergy(int spinAtSite,
                                 double externalField,
                                 const IsingCouplings& couplings,
                                 const std::array<int, 6>& neighborSpinSums) {
    const std::array<double, 6> magneticCouplings = couplings.toArray();
    
    // Accumulate magnetic contribution from all neighbor shells
    double magneticContribution = externalField;
    for (std::size_t shell = 0; shell < magneticCouplings.size(); ++shell) {
        magneticContribution += magneticCouplings[shell] * static_cast<double>(neighborSpinSums[shell]);
    }
    
    // Energy change for spin flip: 2 * spin * (field contribution)
    double dETotal = 2.0 * static_cast<double>(spinAtSite) * magneticContribution;
    
    return dETotal;
}
