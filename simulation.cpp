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
    
    if (params.step_H <= 0) {
        throw std::invalid_argument("Step_H must be positive.");
    }
    
    if (params.H_upper == params.H_lower) {
        list.push_back(params.H_upper);
        return list;
    }
    
    // A. 0 to HUpper
    for (float h = 0.0; h <= params.H_upper; h += params.step_H) {
        list.push_back(h);
    }
    
    // B. HUpper to HLower
    for (float h = params.H_upper; h >= params.H_lower; h -= params.step_H) {
        list.push_back(h);
    }
    
    // C. HLower to 0
    for (float h = params.H_lower + params.step_H; h <= params.H_upper; h += params.step_H) {
        list.push_back(h);
    }

    return list;
}


/**
 * @brief Creates the vector of magnetic field values for the H-sweep loop.
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

void MonteCarloStepChemicalExchange(Lattice& lattice,
                                    float H,
                                    const SimulationParameters& params, 
                                    const SiteEnergyTableBEG& tableBeg,
                                    const FastBoltzmannTableSpin& tableSpin,
                                    float& DeltaEAcumM) {

    float jota1 = 0.25*params.w1_12;
    float jota2 = 0.25*params.w2_12;
    float ka1 = 0.25*(2*params.w1_13 + 2*params.w1_23 - params.w1_12);
    float ka2 = 0.25*(2*params.w2_13 + 2*params.w2_23 - params.w2_12);
    float ele1 = 0.25*(params.w1_13 - params.w1_23);
    float ele2 = 0.25*(params.w2_13 - params.w2_23);

    for (int site = 0; site < lattice.totalSites(); site++) {
        
        int SpecieAct = lattice.getSpecies(site);
        // selects a random neighbor
        int siteNeighbor = lattice.getNeighbors1(site)[RanEnt1a8()-1];
        int SpecieNeigh = lattice.getSpecies(siteNeighbor);
        
        if (SpecieAct == SpecieNeigh) continue; // Skip if same species
        
        // Calculate neighbor sums for 3rd and 6th NN (as Jm1=Jm2=0)
        int SumSpin3N_Act = lattice.calculateNeighborSpinSum(site, 3);
        int SumSpin6N_Act = lattice.calculateNeighborSpinSum(site, 6);
        int SumSpin3N_Neigh = lattice.calculateNeighborSpinSum(siteNeighbor, 3);
        int SumSpin6N_Neigh = lattice.calculateNeighborSpinSum(siteNeighbor, 6);
        
        // Magnetic Energy Difference Calculation
        // Delta EM1: Spin-Spin Interactions
        int SpinAct = lattice.getSpin(site);
        int SpinNeigh = lattice.getSpin(siteNeighbor);
                    
        float dEM = lattice.calculateSiteMagneticEnergy(SpinNeigh, params.Jm3, params.Jm6, H, SumSpin3N_Act, SumSpin6N_Act) +
                    lattice.calculateSiteMagneticEnergy(SpinAct, params.Jm3, params.Jm6, H, SumSpin3N_Neigh, SumSpin6N_Neigh) -
                    lattice.calculateSiteMagneticEnergy(SpinAct, params.Jm3, params.Jm6, H, SumSpin3N_Act, SumSpin6N_Act) -
                    lattice.calculateSiteMagneticEnergy(SpinNeigh, params.Jm3, params.Jm6, H, SumSpin3N_Neigh, SumSpin6N_Neigh);
        
        int SumSpecie1N_Act_linear = lattice.calculateNeighborSpeciesSum(site, 1, 1) - lattice.getSpecies(siteNeighbor);
        int SumSpecie1N_Neigh_linear = lattice.calculateNeighborSpeciesSum(siteNeighbor, 1, 1) - lattice.getSpecies(site);
        int SumSpecie1N_Act_quadratic = lattice.calculateNeighborSpeciesSum(site, 1, 2) - lattice.getSpecies(siteNeighbor)*lattice.getSpecies(siteNeighbor);
        int SumSpecie1N_Neigh_quadratic = lattice.calculateNeighborSpeciesSum(siteNeighbor, 1, 2) - lattice.getSpecies(site)*lattice.getSpecies(site);

        int SumSpecie2N_Act_linear = lattice.calculateNeighborSpeciesSum(site, 2, 1);
        int SumSpecie2N_Neigh_linear = lattice.calculateNeighborSpeciesSum(siteNeighbor, 2, 1);
        int SumSpecie2N_Act_quadratic = lattice.calculateNeighborSpeciesSum(site, 2, 2);
        int SumSpecie2N_Neigh_quadratic = lattice.calculateNeighborSpeciesSum(siteNeighbor, 2, 2);

        float dEChem_i = lattice.calculateSiteChemicalEnergy(SpecieAct,
                                                    jota1, jota2,
                                                    ka1, ka2,
                                                    ele1, ele2,
                                                    SumSpecie1N_Act_linear, SumSpecie1N_Act_quadratic,
                                                    SumSpecie2N_Act_linear, SumSpecie2N_Act_quadratic)+
                        lattice.calculateSiteChemicalEnergy(SpecieNeigh,
                                                    jota1, jota2,
                                                    ka1, ka2,
                                                    ele1, ele2,
                                                    SumSpecie1N_Neigh_linear, SumSpecie1N_Neigh_quadratic,
                                                    SumSpecie2N_Neigh_linear, SumSpecie2N_Neigh_quadratic);
        
        float dEChem_f = lattice.calculateSiteChemicalEnergy(SpecieNeigh,
                                                    jota1, jota2,
                                                    ka1, ka2,
                                                    ele1, ele2,
                                                    SumSpecie1N_Act_linear, SumSpecie1N_Act_quadratic,
                                                    SumSpecie2N_Act_linear, SumSpecie2N_Act_quadratic)+
                        lattice.calculateSiteChemicalEnergy(SpecieAct,
                                                    jota1, jota2,
                                                    ka1, ka2,
                                                    ele1, ele2,
                                                    SumSpecie1N_Neigh_linear, SumSpecie1N_Neigh_quadratic,
                                                    SumSpecie2N_Neigh_linear, SumSpecie2N_Neigh_quadratic);

        float dEChem = dEChem_f - dEChem_i;
        float dETotal = dEM + dEChem;
        
        // Metropolis Algorithm
        if (dETotal > 0) {
            // Accept: Probabilistically
            double epsilonM = Ran0a1();
            
            float BoltzmannM = tableSpin.lookup(SpinNeigh, SumSpin3N_Act, SumSpin6N_Act) * 
                                tableSpin.lookup(SpinAct, SumSpin3N_Neigh, SumSpin6N_Neigh) /
                                (tableSpin.lookup(SpinAct, SumSpin3N_Act, SumSpin6N_Act) * 
                                tableSpin.lookup(SpinNeigh, SumSpin3N_Neigh, SumSpin6N_Neigh));
            
            float BoltzmannChem = tableBeg.lookup(SpecieAct, SumSpecie1N_Act_linear, SumSpecie1N_Act_quadratic, SumSpecie2N_Act_linear, SumSpecie2N_Act_quadratic)*
                                    tableBeg.lookup(SpecieNeigh, SumSpecie1N_Neigh_linear, SumSpecie1N_Neigh_quadratic, SumSpecie2N_Neigh_linear, SumSpecie2N_Neigh_quadratic)/
                                    (tableBeg.lookup(SpecieNeigh, SumSpecie1N_Act_linear, SumSpecie1N_Act_quadratic, SumSpecie2N_Act_linear, SumSpecie2N_Act_quadratic)*
                                    tableBeg.lookup(SpecieAct, SumSpecie1N_Neigh_linear, SumSpecie1N_Neigh_quadratic, SumSpecie2N_Neigh_linear, SumSpecie2N_Neigh_quadratic));
            
            float Boltzmann = std::min(1.0f, BoltzmannM * BoltzmannChem);
            
            if (Boltzmann >= epsilonM) {
                lattice.exchangeSpecies(site, siteNeighbor);
                DeltaEAcumM += dETotal;
            }
            continue;
        }

        // Accept: Exchange species
        lattice.exchangeSpecies(site, siteNeighbor);
        DeltaEAcumM += dETotal;
    
    }
}

/**
 * @brief Executes one full Monte Carlo sweep only for spin with an external field(N spin flip attempts).
 */
void MonteCarloStepSpinExtH(Lattice& lattice, 
                            float H,
                            const SimulationParameters& params, 
                            const FastBoltzmannTableSpin& table,
                            float& DeltaEAcumM) {
    
    for (int site = 0; site < lattice.totalSites(); site++) {
        
        int SpinAct = lattice.getSpin(site);
        
        if (SpinAct == 0) continue; // Skip if spin is 0 (non-magnetic species)

        int Sum3N = lattice.calculateNeighborSpinSum(site, 3);
        int Sum6N = lattice.calculateNeighborSpinSum(site, 6);

        float deltaEM = lattice.calculateSiteMagneticEnergy(-SpinAct, params.Jm3, params.Jm6, H, Sum3N, Sum6N)\
                    - lattice.calculateSiteMagneticEnergy(SpinAct, params.Jm3, params.Jm6, H, Sum3N, Sum6N);
        
        // Metropolis Algorithm
        if (deltaEM > 0) {
        
            // Accept: Probabilistically
            double epsilonM = Ran0a1();
            float BoltzmannMf = table.lookup(-SpinAct, Sum3N, Sum6N);
            float BoltzmannMi = table.lookup(SpinAct, Sum3N, Sum6N);
            
            float Boltzmann = std::min(1.0f, BoltzmannMf / BoltzmannMi);

            if (Boltzmann >= epsilonM) {
        
                lattice.flipSpin(site);
                DeltaEAcumM += deltaEM;
        
            }
        
            continue;
        }
        
        // Accept: Flip spin
        lattice.flipSpin(site);
        DeltaEAcumM += deltaEM;
    
    }
}


/**
 * @brief Main simulation loop, iterating over Temperature and Field.
*/
void SimulationLoop(const SimulationParameters& params, const char* nombrefile) {
    
    Lattice lattice(params.lattice_side); 
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
    lattice.initializeNeighbors();
    vector<float> listaCampos = createHSweepList(params);
    vector<float> listaTemperaturas = createTSweepList(params);
    int output_count = 0; // Counter for final file naming

    for (float T : listaTemperaturas) {

        auto tableBeg = SiteEnergyTableBEG(params.w1_12, params.w2_12,
                                            params.w1_13, params.w2_13,
                                            params.w1_23, params.w2_23,
                                            T);

        for (float H : listaCampos) {
            cout << "Trabajando a T = " << T << " y H = " << H << endl;

            float DeltaEAcumM = 0.0f;
            
            auto tableSpin = FastBoltzmannTableSpin(params.Jm3, params.Jm6, H, T);

            for (int contador = 1; contador <= params.num_steps; contador++) {
                
                // 3a. Single-Spin Update (Metropolis)
                if (params.simulation_method == 1)
                    MonteCarloStepChemicalExchange(lattice, H, params, tableBeg, tableSpin, DeltaEAcumM);
                else if(params.simulation_method == 0) {   
                    MonteCarloStepSpinExtH(lattice, H, params, tableSpin, DeltaEAcumM);
                }

                // 3b. Measurement and Output (Occurs only in the last 200 steps)
                if (contador > (params.num_steps - params.steps_to_output)) {
                    lattice.calculateAndWriteLRO(parout, contador, T, H, DeltaEAcumM);
                }
            }
            
            // 3c. Final Configuration Save
            if (params.flag_save_config) {
                lattice.saveFinalConfiguration(nombrefile, H, T, output_count);
            }
            output_count++;
        }
    }

    parout.close();
}