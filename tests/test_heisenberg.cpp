#include "lattice.h"
#include "rng.h"
#include "simulation.h"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

constexpr double tolerance = 1e-12;

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

double norm(const HeisenbergVector& moment) {
    return std::sqrt(dotProduct(moment, moment));
}

SimulationParameters zeroCouplingParameters() {
    return SimulationParameters(
        1, 3, 1,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
        1, false, false, false);
}

void testMomentStorage() {
    Lattice lattice(1);
    const HeisenbergVector expected{0.25, -0.5, 0.75};

    lattice.setMoment(0, expected);
    const HeisenbergVector& stored = lattice.getMoment(0);
    for (std::size_t component = 0; component < expected.size(); ++component) {
        require(std::abs(stored[component] - expected[component]) < tolerance,
                "Lattice did not preserve a stored moment component");
    }
}

void testProposalNormalizationAndCentering() {
    const double z = std::sqrt(1.0 - 0.3 * 0.3 - 0.4 * 0.4);
    const HeisenbergVector selectedMoment{0.3, 0.4, z};
    constexpr int sampleCount = 20000;
    constexpr double sigma = 0.1;

    HeisenbergVector mean{0.0, 0.0, 0.0};
    for (int sample = 0; sample < sampleCount; ++sample) {
        const HeisenbergVector proposal = gaussianSpinProposal(selectedMoment, sigma);
        require(std::abs(norm(proposal) - 1.0) < tolerance,
                "Gaussian proposal is not normalized");
        for (std::size_t component = 0; component < proposal.size(); ++component) {
            mean[component] += proposal[component];
        }
    }

    for (double& component : mean) {
        component /= static_cast<double>(sampleCount);
    }
    require(std::abs(mean[0] - selectedMoment[0]) < 0.02,
            "Gaussian proposal mean is not centered on the selected x component");
    require(std::abs(mean[1] - selectedMoment[1]) < 0.02,
            "Gaussian proposal mean is not centered on the selected y component");
    require(std::abs(mean[2] - selectedMoment[2]) < 0.02,
            "Gaussian proposal mean is not centered on the selected z component");
}

void testHeisenbergSweepStoresNormalizedMoments() {
    Lattice lattice(1);
    lattice.initializeNeighbors();
    lattice.setMoment(0, {1.0, 0.0, 0.0});
    lattice.setMoment(1, {0.0, 1.0, 0.0});

    SimulationParameters params = zeroCouplingParameters();
    std::vector<double> energies;
    BoltzmannDeltaETable table(energies, 1.0);
    double sigma = 0.1;
    std::uint64_t previousAccepted = 0;
    std::uint64_t previousAttempted = 0;
    double accumulatedEnergy = 0.0;
    std::uint64_t accepted = 0;
    std::uint64_t attempted = 0;
    MCStepResults stats(accumulatedEnergy, accepted, attempted);

    MonteCarloStepHeisenberg(lattice, 0.0, params, table, sigma,
                             previousAccepted, previousAttempted, stats);

    require(attempted == 2 && accepted == 2,
            "Heisenberg sweep did not accept both zero-energy proposals");
    require(previousAttempted == 2 && previousAccepted == 2,
            "Heisenberg sweep counters are inconsistent");
    for (int site = 0; site < lattice.totalSites(); ++site) {
        require(std::abs(norm(lattice.getMoment(site)) - 1.0) < tolerance,
                "Heisenberg sweep stored a non-normalized moment");
    }
}

} // namespace

int main() {
    try {
        testMomentStorage();
        testProposalNormalizationAndCentering();
        testHeisenbergSweepStoresNormalizedMoments();
    } catch (const std::exception& error) {
        std::cerr << "Heisenberg test failed: " << error.what() << '\n';
        return 1;
    }

    std::cout << "Heisenberg tests passed\n";
    return 0;
}
