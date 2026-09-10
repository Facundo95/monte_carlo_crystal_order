#include "heisenberg_hamiltonian.h"

#include <cstddef>

namespace {

constexpr std::size_t zComponent = 2;

}

double dotProduct(const HeisenbergVector& first, const HeisenbergVector& second) {
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
}

double calculateDeltaHeisenbergEnergy(const HeisenbergVector& currentMoment,
                                      const HeisenbergVector& proposedMoment,
                                      double externalField,
                                      const HeisenbergCouplings& couplings,
                                      const HeisenbergNeighborSums& neighborSums) {
    const std::array<double, 6> magneticCouplings = couplings.toArray();
    double deltaEnergy = -externalField *
                         (proposedMoment[zComponent] - currentMoment[zComponent]);

    for (std::size_t shell = 0; shell < magneticCouplings.size(); ++shell) {
        deltaEnergy -= magneticCouplings[shell] *
                       (dotProduct(proposedMoment, neighborSums[shell]) -
                        dotProduct(currentMoment, neighborSums[shell]));
    }

    return deltaEnergy;
}
