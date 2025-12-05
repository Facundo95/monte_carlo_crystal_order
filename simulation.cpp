#include "simulation.h"
#include "rng.h"
#include "file_handler.h"
#include <stdexcept>
#include <cstring>
#include <cstdint>

/**
 * @brief Creates the vector of values for the sweep list.
 */
std::vector<float> createSweepList(float start, float end, float step, bool loop) {
    std::vector<float> list;
    
    if (step <= 0) {
        throw std::invalid_argument("Step must be positive.");
    }
    
    if (start == end) {
        list.push_back(start);
        return list;
    }

    if (!loop) {
        if (start < end) {
            for (float val = start; val <= end; val += step) {
                list.push_back(val);
            }
        } else {
            for (float val = start; val >= end; val -= step) {
                list.push_back(val);
            }
        }
        return list;
    }
    // Looping case
    if (start < end) {
        for (float val = start; val <= end; val += step) {
            list.push_back(val);
        }
        for (float val = end - step; val >= start; val -= step) {
            list.push_back(val);
        }
    } else {
        for (float val = start; val >= end; val -= step) {
            list.push_back(val);
        }
        for (float val = end + step; val <= start; val += step) {
            list.push_back(val);
        }
    }
    return list;
}


/**
 * @brief Creates the vector of temperatures values for the T-sweep loop.
 */
std::vector<float> createTSweepList(const SimulationParameters& params) {
    std::vector<float> list;
    
    if (params.step_T <= 0) {
        throw std::invalid_argument("Step_T must be positive.");
    }
    
    if (params.T_start == params.T_end) {
        list.push_back(params.T_start);
        return list;
    }

    if (params.T_start < params.T_end) {
        for (float t = params.T_start; t <= params.T_end; t += params.step_T) {
            list.push_back(t);
        }
        return list;
    }

    for (float t = params.T_start; t >= params.T_end; t -= params.step_T) {
        list.push_back(t);
    }

    return list;
}

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
                                    float& DeltaEAcumM,
                                    int& changesAccepted,
                                    int& changesAttempted) {

    float jota1 = 0.25 * params.w1_13;
    float jota2 = 0.25 * params.w2_13;
    float ka1 = 0.25 * (2 * params.w1_12 + 2 * params.w1_23 - params.w1_13);
    float ka2 = 0.25 * (2 * params.w2_12 + 2 * params.w2_23 - params.w2_13);
    float ele1 = 0.25 * (params.w1_12 - params.w1_23);
    float ele2 = 0.25 * (params.w2_12 - params.w2_23);

    for (int site = 0; site < lattice.totalSites(); site++) {
        
        int SpecieAct = lattice.getSpecies(site);
        // selects a random neighbor
        int siteNeighbor = lattice.getNeighbors1(site)[RanEnt1a8()-1];
        int SpecieNeigh = lattice.getSpecies(siteNeighbor);
        
        if (SpecieAct == SpecieNeigh) continue; // Skip if same species

        changesAttempted++;
        
        int SumLinNN_A = lattice.calculateNeighborSpeciesSum(site, 1, 1) - lattice.getSpecies(siteNeighbor);
        int SumLinNN_N = lattice.calculateNeighborSpeciesSum(siteNeighbor, 1, 1) - lattice.getSpecies(site);
        int SumCuadNN_A = lattice.calculateNeighborSpeciesSum(site, 1, 2) - lattice.getSpecies(siteNeighbor)*lattice.getSpecies(siteNeighbor);
        int SumCuadNN_N = lattice.calculateNeighborSpeciesSum(siteNeighbor, 1, 2) - lattice.getSpecies(site)*lattice.getSpecies(site);

        int SumLinNNN_A = lattice.calculateNeighborSpeciesSum(site, 2, 1);
        int SumLinNNN_N = lattice.calculateNeighborSpeciesSum(siteNeighbor, 2, 1);
        int SumCuadNNN_A = lattice.calculateNeighborSpeciesSum(site, 2, 2);
        int SumCuadNNN_N = lattice.calculateNeighborSpeciesSum(siteNeighbor, 2, 2);

        double dETotal = lattice.calculateDeltaChemicalEnergy(SpecieAct, SpecieNeigh,
                                                            jota1, jota2,
                                                            ka1, ka2,
                                                            ele1, ele2,
                                                            SumLinNN_A, SumLinNNN_A,
                                                            SumCuadNN_A, SumCuadNNN_A,
                                                            SumLinNN_N, SumLinNNN_N,
                                                            SumCuadNN_N, SumCuadNNN_N);
        
        // Metropolis Algorithm
        if (dETotal > 0) {
            // Accept: Probabilistically
            double epsilon = Ran0a1();
            
            float BoltzmannChem = table.lookup(dETotal);
            
            if (BoltzmannChem == 0.0f) {
                std::cout << "Warning: dE not found in Boltzmann table: dE = " << dETotal << std::endl;
                std::cout << "Species on site: " << SpecieAct << ", Species at neighbor " << siteNeighbor << ": " << SpecieNeigh << std::endl;
                std::cout << "SumLinNN_A: " << SumLinNN_A << ", SumLinNNN_A: " << SumLinNNN_A
                          << ", SumCuadNN_A: " << SumCuadNN_A << ", SumCuadNNN_A: " << SumCuadNNN_A << std::endl;
                std::cout << "SumLinNN_N: " << SumLinNN_N << ", SumLinNNN_N: " << SumLinNNN_N
                          << ", SumCuadNN_N: " << SumCuadNN_N << ", SumCuadNNN_N: " << SumCuadNNN_N << std::endl;
                std::cout << "----------------------------------------" << std::endl;
            }

            if (BoltzmannChem >= epsilon) {
                lattice.exchangeSpecies(site, siteNeighbor);
                changesAccepted++;
                DeltaEAcumM += dETotal;
            }
            continue;
        }

        // Accept: Exchange species
        lattice.exchangeSpecies(site, siteNeighbor);
        changesAccepted++;
        DeltaEAcumM += dETotal;
    
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
                            float H,
                            const SimulationParameters& params, 
                            BoltzmannDeltaETable& table,
                            float& DeltaEAcumM,
                            int& changesAccepted,
                            int& changesAttempted) {
    
    for (int site = 0; site < lattice.totalSites(); site++) {
        
        int SpinAct = lattice.getSpin(site);
        
        if (SpinAct == 0) continue; // Skip if spin is 0 (non-magnetic species)

        changesAttempted++;

        int Sum3N = lattice.calculateNeighborSpinSum(site, 3);
        int Sum6N = lattice.calculateNeighborSpinSum(site, 6);

        float dETotal = 2 * SpinAct * (params.Jm3 * Sum3N + params.Jm6 * Sum6N + H);
        
        // Metropolis Algorithm
        if (dETotal > 0) {
        
            // Accept: Probabilistically
            double epsilon = Ran0a1();
            float Boltzmann = table.lookup(dETotal);

            if (Boltzmann == 0.0f) {
                std::cout << "Warning: dE not found in Boltzmann table: dE = " << dETotal << std::endl;
            }

            if (Boltzmann >= epsilon) {
        
                lattice.flipSpin(site);
                changesAccepted++;
                DeltaEAcumM += dETotal;
        
            }
        
            continue;
        }
        
        // Accept: Flip spin
        lattice.flipSpin(site);
        changesAccepted++;
        DeltaEAcumM += dETotal;
    
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
    // 1. Initialization
    try {
        lattice.loadInitialConfiguration(file_in);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return;
    }

    // 2. Setup Output
    std::ofstream parout;
    if (!OpenLROParametersFile(file_out, parout)) {
        return; 
    }
    
    // 3. Monte Carlo Loop Setup
    lattice.initializeNeighbors();

    if (params.simulation_method != 0 && params.simulation_method != 1) {
        std::cerr << "ERROR: Método de simulación desconocido: " << params.simulation_method << std::endl;
        return;
    }

    
    std::vector<float> listaCampos = createSweepList(params.H_start, params.H_end, params.step_H, params.flag_loop);
    std::vector<float> listaTemperaturas = createSweepList(params.T_start, params.T_end, params.step_T, false);
    int output_count = 0; // Counter for final file naming

    for (float T : listaTemperaturas) {

        std::vector<double> dEs = {};
        auto table = BoltzmannDeltaETable(dEs, T);
        
        for (float H: listaCampos) {
            
            if (verbose) {
                std::cout << "----------------------------------------" << std::endl;
                std::cout << "Trabajando a T = " << T << " y H = " << H << std::endl;
            }

            int changesAccepted = 0;
            int changesAttempted = 0;
            float DeltaEAcumM = 0.0f;

            for (int contador = 1; contador <= params.num_steps; contador++) {
                
                if (params.simulation_method == 0) {
                    MonteCarloStepChemicalExchange(lattice, params, table, DeltaEAcumM, changesAccepted, changesAttempted);
                } else if (params.simulation_method == 1) {
                    MonteCarloStepSpinExtH(lattice, H, params, table, DeltaEAcumM, changesAccepted, changesAttempted);
                }
                    // 3b. Measurement and Output
                if (contador > (params.num_steps - params.steps_to_output)) {
                    lattice.calculateAndWriteLRO(parout, contador, T, H, DeltaEAcumM);
                }
            }

            // 3c. Final Configuration Save
            if (params.flag_save_config) {
                bool ok = lattice.saveFinalConfiguration(file_out, H, T, output_count);
                if (!ok) std::cerr << "WARNING: could not save final configuration for output_count=" << output_count << std::endl;
            }

            if (verbose) {
                std::cout << "Intercambios aceptados / intentado: " << changesAccepted << "/" << params.num_steps << " en " << params.num_steps << " pasos." << std::endl;
            }

            output_count++;

        }
    }

    parout.close();

}