#include "ising_hamiltonian.h"
#include "lattice.h"
#include "neighbor_sums.h"

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

double calculateDeltaIsingEnergyForExchange(const Lattice& lattice,
                                            int siteA,
                                            int siteB,
                                            const IsingCouplings& couplings) {
    const std::array<double, 6> magneticCouplings = couplings.toArray();
    const double spinA = static_cast<double>(lattice.getSpin(siteA));
    const double spinB = static_cast<double>(lattice.getSpin(siteB));
    double deltaEnergy = 0.0;

    for (int shell = 1; shell <= 6; ++shell) {
        const double sumAroundA = static_cast<double>(
            lattice.calculateNeighborSpinSumExcluding(siteA, shell, siteB));
        const double sumAroundB = static_cast<double>(
            lattice.calculateNeighborSpinSumExcluding(siteB, shell, siteA));

        deltaEnergy -= magneticCouplings[shell - 1] *
                       ((spinB - spinA) * sumAroundA +
                        (spinA - spinB) * sumAroundB);
    }

    return deltaEnergy;
}

double calculateTotalIsingEnergy(const Lattice& lattice, const IsingCouplings& couplings, double externalField) {
    double totalMagneticE = 0.0;
    double totalMagnetization = 0.0;

    for (int site = 0; site < lattice.totalSites(); ++site) {
        int Mi = lattice.getSpin(site);
        totalMagnetization += static_cast<double>(Mi);

        if (Mi != 0) {
            auto sums = computeNeighborSpinSums(lattice, site);
            // map sums to individual shells
            double sum1 = static_cast<double>(sums[0]);
            double sum2 = static_cast<double>(sums[1]);
            double sum3 = static_cast<double>(sums[2]);
            double sum4 = static_cast<double>(sums[3]);
            double sum5 = static_cast<double>(sums[4]);
            double sum6 = static_cast<double>(sums[5]);

            totalMagneticE += lattice.calculateSiteMagneticEnergy(Mi,
                                                                   couplings.Jm1, couplings.Jm2, couplings.Jm3,
                                                                   couplings.Jm4, couplings.Jm5, couplings.Jm6,
                                                                   externalField,
                                                                   sum1, sum2, sum3, sum4, sum5, sum6);
        }
    }

    double fieldEnergy = - totalMagnetization * externalField;
    double interactionMagneticE = totalMagneticE - fieldEnergy;
    return 0.5 * interactionMagneticE + fieldEnergy;
}
