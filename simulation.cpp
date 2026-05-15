#include "simulation.h"
#include "rng.h"
#include "file_handler.h"
#include <stdexcept>
#include <cstring>
#include <cstdint>

/** @brief Performs a Monte Carlo step using chemical species exchange dynamics. 
 * @param lattice The lattice object representing the system.
 * @param params The simulation parameters.
 * @param tableBeg The pre-computed BEG site energy table.
 * @param DeltaEAcumM Accumulated energy change for magnetization.
 * @param changesAccepted Counter for accepted changes.
*/
void MonteCarloStepChemicalExchange(Lattice& lattice,
                                    const SimulationParameters& params,
                                    BoltzmannDeltaETable& table,
                                    MCStepResults& stats) {

    for (int site = 0; site < lattice.totalSites(); site++) {
        
        int SpecieAct = lattice.getSpecies(site);
        // selects a random neighbor
        int siteNeighbor = lattice.getNeighbors1(site)[RanEnt1a8()-1];
        int SpecieNeigh = lattice.getSpecies(siteNeighbor);
        
        if (SpecieAct == SpecieNeigh) continue; // Skip if same species

        stats.changesAttempted++;
        
        auto sums = computeNeighborSpeciesSums(lattice, site, siteNeighbor);

        double dETotal = lattice.calculateDeltaChemicalEnergy(SpecieAct, SpecieNeigh,
                                    params.jota1, params.jota2,
                                    params.ka1, params.ka2,
                                    params.ele1, params.ele2,
                                    sums.linNN_A, sums.linNNN_A,
                                    sums.cuadNN_A, sums.cuadNNN_A,
                                    sums.linNN_N, sums.linNNN_N,
                                    sums.cuadNN_N, sums.cuadNNN_N);
        
        // Metropolis acceptance
        if (metropolisAccept(dETotal, table)) {
            lattice.exchangeSpecies(site, siteNeighbor);
            stats.changesAccepted++;
            stats.DeltaEAcum += dETotal;

        }
    }
}

/**
 * @brief Executes one full Monte Carlo sweep only for spin with an external field(N spin flip attempts).
 * @param lattice The lattice object representing the system.
 * @param H The external magnetic field.
 * @param params The simulation parameters.
 * @param table The pre-computed spin Boltzmann table.
 * @param DeltaEAcumM Accumulated energy change for magnetization.
 * @param changesAccepted Counter for accepted changes.
 */
void MonteCarloStepSpinExtH(Lattice& lattice,
                            double H,
                            const SimulationParameters& params,
                            BoltzmannDeltaETable& table,
                            MCStepResults& stats) {

    for (int site = 0; site < lattice.totalSites(); site++) {

        int SpinAct = lattice.getSpin(site);

        if (SpinAct == 0) continue; // Skip if spin is 0 (non-magnetic species)

        stats.changesAttempted++;

        const std::array<int, 6> neighborSums = computeNeighborSpinSums(lattice, site);
        const std::array<double, 6> magneticCouplings = {
            params.Jm1, params.Jm2, params.Jm3, params.Jm4, params.Jm5, params.Jm6
        };

        double magneticContribution = H;
        for (std::size_t shell = 0; shell < magneticCouplings.size(); ++shell) {
            magneticContribution += magneticCouplings[shell] * static_cast<double>(neighborSums[shell]);
        }

        double dETotal = 2.0 * static_cast<double>(SpinAct) * magneticContribution;
        
        // Metropolis acceptance
        if (metropolisAccept(dETotal, table)) {
            lattice.flipSpin(site);
            stats.changesAccepted++;
            stats.DeltaEAcum += dETotal;
        }
    }
}


/**
 * @brief Main simulation loop, iterating over Temperature and Field.
 * @param params The simulation parameters.
 * @param nombrefile The base name for output files.
*/
void SimulationLoop(const SimulationParameters& params, 
                    const char* file_in, 
                    const char* file_out) {
    
    bool verbose = false;

    Lattice lattice(params.lattice_side);
    // set element symbols from params (if any)
    lattice.setAtomNames(params.atom1, params.atom2, params.atom3);
    // 1. Initialization
    try {
        lattice.loadInitialConfiguration(file_in);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return;
    }

    // 2. Setup Output
    std::ofstream parout;
    if (!OpenLROParametersFile(file_out, parout, lattice)) {
        return; 
    }
    
    // 3. Monte Carlo Loop Setup
    lattice.initializeNeighbors();

    if (params.simulation_method != 0 && params.simulation_method != 1 && params.simulation_method != 2) {
        std::cerr << "ERROR: Método de simulación desconocido: " << params.simulation_method << std::endl;
        return;
    }

    
    std::vector<double> listaCampos = createSweepList(params.H_start, params.H_end, params.step_H, params.flag_loop);
    std::vector<double> listaTemperaturas = createSweepList(params.T_start, params.T_end, params.step_T, false);
    int output_count = 0; // Counter for final file naming

    for (double T : listaTemperaturas) {

        std::vector<double> dEs = {};
        auto table = BoltzmannDeltaETable(dEs, T);
        
        for (double H: listaCampos) {
            double currentTotalEnergy = lattice.calculateTotalEnergy(params, H);

            if (verbose) {
                std::cout << "----------------------------------------" << std::endl;
                std::cout << "Trabajando a T = " << T << " y H = " << H << std::endl;
            }

            int spinChangesAccepted = 0;
            int chemicalChangesAccepted = 0;
            int spinChangesAttempted = 0;
            int chemicalChangesAttempted = 0;
            double DeltaEAcumM = 0.0;
            double DeltaEAcumC = 0.0;

            // Bundle counters into MCStepResults for cleaner function calls
            MCStepResults chemStats(DeltaEAcumC, chemicalChangesAccepted, chemicalChangesAttempted);
            MCStepResults spinStats(DeltaEAcumM, spinChangesAccepted, spinChangesAttempted);

            for (int contador = 1; contador <= params.num_steps; contador++) {
                
                if (params.simulation_method == 0) {
                    MonteCarloStepChemicalExchange(lattice, params, table, chemStats);
                } else if (params.simulation_method == 1) {
                    MonteCarloStepSpinExtH(lattice, H, params, table, spinStats);
                } else if (params.simulation_method == 2) {
                    MonteCarloStepChemicalExchange(lattice, params, table, chemStats);
                    MonteCarloStepSpinExtH(lattice, H, params, table, spinStats);
                }

                double energyAtStep = currentTotalEnergy + DeltaEAcumM + DeltaEAcumC;

                    // 3b. Measurement and Output
                if (contador > (params.num_steps - params.steps_to_output)) {
                    lattice.calculateAndWriteLRO(parout, contador, T, H, energyAtStep);
                }
            }

            // 3c. Final Configuration Save
            if (params.flag_save_config) {
                bool ok = lattice.saveFinalConfiguration(file_out, H, T, output_count);
                if (!ok) std::cerr << "WARNING: could not save final configuration for output_count=" << output_count << std::endl;
            }

            if (verbose) {
                std::cout << "Intercambios de espines aceptados / intentado: " << spinChangesAccepted << "/" << spinChangesAttempted << " en " << params.num_steps << " pasos." << std::endl;
                std::cout << "Intercambios químicos aceptados / intentado: " << chemicalChangesAccepted << "/" << chemicalChangesAttempted << " en " << params.num_steps << " pasos." << std::endl;
            }

            output_count++;

        }
    }

    parout.close();

}