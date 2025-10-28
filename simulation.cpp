#include "simulation.h"
#include "rng.h"
#include "file_handler.h"
#include <stdexcept>

using namespace std;

// --- Main Simulation Flow Functions ---

/**
 * @brief Creates the vector of magnetic field values for the H-sweep loop.
 */
std::vector<float> createHSweepList(const SimulationParameters& params) {
    std::vector<float> list;

    // A. 0 to HUpper
    for (float h = 0.0; h <= params.H_upper; h += params.step_H) {
        list.push_back(h);
    }
    // B. HUpper to HLower
    for (float h = params.H_upper; h >= params.H_lower; h -= params.step_H) {
        list.push_back(h);
    }
    // C. HLower to 0 (or HUpper, original code was ambiguous, using HLower to HUpper)
    // Adjusting to match the user's explicit original loop bounds:
    for (float h = params.H_lower + params.step_H; h <= params.H_upper; h += params.step_H) {
        list.push_back(h);
    }

    return list;
}

void MonteCarloStepChemicalExchange(Lattice& lattice, 
                                    const SimulationParameters& params, 
                                    const FastBoltzmannTable& table,
                                    float& DeltaEAcumM) {
    lattice.initializeNeighbors();
    for (int site = 0; site < LATTICE_TOTAL_SITES; site++) {
        int SpecieAct = lattice.getSpecies(site);
        int SpecieNeigh = lattice.getNeighbors1(site)[RanEnt1a8()];
        if (SpecieAct != SpecieNeigh) {
        // placeholder for chemical exchange logic
        // Currently not implemented in the original code snippet
        }

    }
}

/**
 * @brief Executes one full Monte Carlo sweep only for spin with an external field(N spin flip attempts).
 */
void MonteCarloStepSpinExtH(Lattice& lattice, 
                            float H, 
                            const SimulationParameters& params, 
                            const FastBoltzmannTable& table,
                            float& DeltaEAcumM) {
    
    lattice.initializeNeighbors();
    for (int site = 0; site < LATTICE_TOTAL_SITES; site++) {
        int SpinAct = lattice.getSpin(site);
        
        // Calculate neighbor sums for 3rd and 6th NN (as Jm1=Jm2=0)
        int Sum3N = lattice.calculateNeighborSpinSum(site, 3);
        int Sum6N = lattice.calculateNeighborSpinSum(site, 6);
        
        // Magnetic Energy Difference Calculation
        // Delta EM1: Spin-Spin Interactions
        float deltaEM1 = 2 * params.Jm3 * SpinAct * Sum3N 
                        + 2 * params.Jm6 * SpinAct * Sum6N;
        
        // Delta EM2: Field Interaction
        float deltaEM2 = 2 * H * SpinAct;
        float deltaEM = deltaEM1 + deltaEM2;

        // Metropolis Algorithm
        if (SpinAct != 0) {
            if (deltaEM <= 0) {
                // Accept: Flip spin
                lattice.flipSpin(site);
                DeltaEAcumM += deltaEM;
            } else {
                // Accept: Probabilistically
                double epsilonM = Ran0a1();
                float BoltzmannM = table.lookup(SpinAct, Sum3N, Sum6N);
                if (BoltzmannM >= epsilonM) {
                    lattice.flipSpin(site);
                    DeltaEAcumM += deltaEM;
                }
            }
        }
    }
}


/**
 * @brief Main simulation loop, iterating over Temperature and Field.
*/
void SimulationLoop(const SimulationParameters& params, const char* nombrefile) {
    
    Lattice lattice; 
    string init_file = string(nombrefile) + ".txt";
    
    // 1. Initialization
    try {
        lattice.loadInitialConfiguration(init_file);
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return;
    }

    // 2. Setup Output
    ofstream parout;
    if (!OpenLROParametersFile(nombrefile, parout)) {
        return; 
    }
    
    // 3. Monte Carlo Loop Setup
    vector<float> listaCampos = createHSweepList(params);
    int H_count = 0; // Counter for final file naming

    for (float TEMPERA = params.T_lower; TEMPERA <= params.T_upper; TEMPERA += params.step_T) {
        for (float Hache : listaCampos) {
            cout << "Trabajando a T = " << TEMPERA << " y H = " << Hache << endl;

            float DeltaEAcumM = 0;
            auto table = FastBoltzmannTable(params.Jm3, params.Jm6, Hache, TEMPERA);

            for (int contador = 1; contador <= params.num_steps; contador++) {
                
                // 3a. Single-Spin Update (Metropolis)
                MonteCarloStepSpinExtH(lattice, Hache, params, table, DeltaEAcumM);

                // 3b. Measurement and Output (Occurs only in the last 200 steps)
                if (contador > (params.num_steps - 200)) {
                    lattice.calculateAndWriteLRO(parout, contador, TEMPERA, Hache, DeltaEAcumM);
                }
            }
            
            // 3c. Final Configuration Save
            if (params.flag_save_config) {
                lattice.saveFinalConfiguration(nombrefile, Hache, H_count);
            }
            H_count++;
        }
    }

    parout.close();
}