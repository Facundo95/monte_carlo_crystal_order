#include "simulation.h"
#include "rng.h"
#include "file_handler.h"
#include <stdexcept>

// --- Utility Functions for Lattice ---

// Wraps X/Y coordinates using periodic boundary conditions
int Lattice::wrapCoord(int coord, int max_size) const {
    if (coord < 0) return coord + max_size;
    if (coord >= max_size) return coord - max_size;
    return coord;
}

// Wraps Z coordinate (size 2*side) using periodic boundary conditions
int Lattice::wrapCoordZ(int coord) const {
    if (coord < 0) return coord + LATTICE_DEPTH;
    if (coord >= LATTICE_DEPTH) return coord - LATTICE_DEPTH;
    return coord;
}

/**
 * @brief Loads the initial configuration from a .txt file.
 */
void Lattice::loadInitialConfiguration(const std::string& filename) {
    std::ifstream redin(filename, std::ios::in | std::ios::binary);
    if (!redin.is_open()) {
        throw std::runtime_error("No se pudo abrir el archivo de entrada inicial: " + filename);
    }

    float aux;
    redin.seekg(0, std::ios::beg);
    for (int i = 0; i < LATTICE_SIDE; i++) {
        for (int j = 0; j < LATTICE_SIDE; j++) {
            for (int k = 0; k < LATTICE_DEPTH; k++) {
                if (!(redin >> aux)) {
                    throw std::runtime_error("Error de lectura o archivo incompleto en la red inicial.");
                }
                // Polynomial transformation from the original code
                red[i][j][k] = (-1. / 12) * aux * aux * aux + (1. / 3) * aux * aux + (7. / 12) * aux - (5. / 6);
                magn[i][j][k] = (-1. / 12) * aux * aux * aux - (1. / 3) * aux * aux + (7. / 12) * aux + (5. / 6);
            }
        }
    }
    redin.close();
}

/**
 * @brief Calculates the sum of spins for a specific neighbor shell (3rd or 6th NN).
 * @param shell_type Must be 3 (3rd NN) or 6 (6th NN).
 * (Original Jm1 and Jm2 (1st and 2nd NN) were 0, so only 3 and 6 are implemented here).
 */
float Lattice::calculateNeighborSpinSum(int XLoc, int YLoc, int ZLoc, int shell_type) const {
    float sum = 0.0;
    
    if (shell_type == 3) {
        // 3rd Neighbors (12 sites)
        // (X+-1, Y+-1, Z), (X, Y+-1, Z+-1), (X+-1, Y, Z+-1) -- equivalent to FCC face diagonals
        
        // Neighbors in XY plane
        sum += magn[wrapCoord(XLoc + 1, LATTICE_SIDE)][wrapCoord(YLoc + 1, LATTICE_SIDE)][ZLoc];
        sum += magn[wrapCoord(XLoc + 1, LATTICE_SIDE)][wrapCoord(YLoc - 1, LATTICE_SIDE)][ZLoc];
        sum += magn[wrapCoord(XLoc - 1, LATTICE_SIDE)][wrapCoord(YLoc + 1, LATTICE_SIDE)][ZLoc];
        sum += magn[wrapCoord(XLoc - 1, LATTICE_SIDE)][wrapCoord(YLoc - 1, LATTICE_SIDE)][ZLoc];
        
        // Neighbors in YZ plane
        sum += magn[XLoc][wrapCoord(YLoc + 1, LATTICE_SIDE)][wrapCoordZ(ZLoc + 2)]; // Z+2, Z-2 are 3rd NN
        sum += magn[XLoc][wrapCoord(YLoc - 1, LATTICE_SIDE)][wrapCoordZ(ZLoc + 2)];
        sum += magn[XLoc][wrapCoord(YLoc + 1, LATTICE_SIDE)][wrapCoordZ(ZLoc - 2)];
        sum += magn[XLoc][wrapCoord(YLoc - 1, LATTICE_SIDE)][wrapCoordZ(ZLoc - 2)];

        // Neighbors in XZ plane
        sum += magn[wrapCoord(XLoc + 1, LATTICE_SIDE)][YLoc][wrapCoordZ(ZLoc + 2)];
        sum += magn[wrapCoord(XLoc - 1, LATTICE_SIDE)][YLoc][wrapCoordZ(ZLoc + 2)];
        sum += magn[wrapCoord(XLoc + 1, LATTICE_SIDE)][YLoc][wrapCoordZ(ZLoc - 2)];
        sum += magn[wrapCoord(XLoc - 1, LATTICE_SIDE)][YLoc][wrapCoordZ(ZLoc - 2)];
        
    } else if (shell_type == 6) {
        // 6th Neighbors (6 sites) - (X+-2, Y, Z), (X, Y+-2, Z), (X, Y, Z+-4)
        sum += magn[wrapCoord(XLoc + 2, LATTICE_SIDE)][YLoc][ZLoc];
        sum += magn[wrapCoord(XLoc - 2, LATTICE_SIDE)][YLoc][ZLoc];
        sum += magn[XLoc][wrapCoord(YLoc + 2, LATTICE_SIDE)][ZLoc];
        sum += magn[XLoc][wrapCoord(YLoc - 2, LATTICE_SIDE)][ZLoc];
        sum += magn[XLoc][YLoc][wrapCoordZ(ZLoc + 4)];
        sum += magn[XLoc][YLoc][wrapCoordZ(ZLoc - 4)];
    }
    // Note: The original implementation in PasoMC defined SumaLin3NMagnActual 
    // using neighbors that seem to correspond to the 3rd NN shell.
    // The implementation here follows the structure suggested by the original indices (X+-1, Y+-1, Z, etc. and Z+-2)
    
    return sum;
}

/**
 * @brief Calculates and writes the LRO parameters and Magnetization to the output file.
 */
void Lattice::calculateAndWriteLRO(std::ofstream& parout, int step_count, float T, float H, float DeltaEAcumM) const {
    
    int CuI=0, CuII=0, CuIII=0, CuIV=0;
    int MnUpI=0, MnUpII=0, MnUpIII=0, MnUpIV=0, MnDownI=0, MnDownII=0, MnDownIII=0, MnDownIV=0; 
    int AlI=0, AlII=0, AlIII=0, AlIV=0;

    long Magnetizacion = 0;
    
    // Traversal and counting for LRO parameters
    for (int i=0; i<LATTICE_SIDE;i++){
        for (int j=0; j<LATTICE_SIDE; j++){
            for (int k=0; k<LATTICE_DEPTH; k++){
                
                Magnetizacion += magn[i][j][k];

                // Determine sublattice (I, II, III, or IV)
                bool Z_is_even = (k % 2 == 0);
                bool XY_sum_is_even = ((i + j + (k/2)) % 2 == 0);

                int* Cu_ptr = nullptr;
                int* MnUp_ptr = nullptr;
                int* MnDown_ptr = nullptr;
                int* Al_ptr = nullptr;
                
                // Simplified conditional logic using pointers
                if (Z_is_even) {
                    if (XY_sum_is_even) { // Sublattice I (Z even, X+Y+Z/2 even)
                        Cu_ptr = &CuI; MnUp_ptr = &MnUpI; MnDown_ptr = &MnDownI; Al_ptr = &AlI;
                    } else { // Sublattice II (Z even, X+Y+Z/2 odd)
                        Cu_ptr = &CuII; MnUp_ptr = &MnUpII; MnDown_ptr = &MnDownII; Al_ptr = &AlII;
                    }
                } else {
                    if (XY_sum_is_even) { // Sublattice III (Z odd, X+Y+Z/2 even)
                        Cu_ptr = &CuIII; MnUp_ptr = &MnUpIII; MnDown_ptr = &MnDownIII; Al_ptr = &AlIII;
                    } else { // Sublattice IV (Z odd, X+Y+Z/2 odd)
                        Cu_ptr = &CuIV; MnUp_ptr = &MnUpIV; MnDown_ptr = &MnDownIV; Al_ptr = &AlIV;
                    }
                }
                
                // Increment counters based on species (red[i][j][k] = 1, 0, or -1)
                if (red[i][j][k] == 1) { // Copper
                    (*Cu_ptr)++;
                } else if (red[i][j][k] == 0) { // Manganese
                    if (magn[i][j][k] == 1) {
                        (*MnUp_ptr)++;
                    } else {
                        (*MnDown_ptr)++;
                    }
                } else if (red[i][j][k] == -1) { // Aluminum
                    (*Al_ptr)++;
                }
            }
        }
    }
    
    // Calculate LRO Parameters (X, Y, Z) and normalized Magnetization
    const float ENE = (float)LATTICE_TOTAL_SITES;
    
    // X parameters (I+II vs III+IV)
    float Xcu = (float)((CuI+CuII-CuIII-CuIV)/(ENE));
    float Xmnup = (float)((MnUpI+MnUpII-MnUpIII-MnUpIV)/(ENE));
    float Xmndown = (float)((MnDownI+MnDownII-MnDownIII-MnDownIV)/(ENE));
    float Xal = (float)((AlI+AlII-AlIII-AlIV)/(ENE));
    
    // Y parameters (I vs II)
    float Ycu = (float)((CuI-CuII)*2/(ENE));
    float Ymnup = (float)((MnUpI-MnUpII)*2/(ENE));
    float Ymndown = (float)((MnDownI-MnDownII)*2/(ENE));
    float Yal = (float)((AlI-AlII)*2/(ENE));
    
    // Z parameters (III vs IV)
    float Zcu = (float)((CuIII-CuIV)*2/(ENE));
    float Zmnup = (float)((MnUpIII-MnUpIV)*2/(ENE));
    float Zmndown = (float)((MnDownIII-MnDownIV)*2/(ENE));
    float Zal = (float)((AlIII-AlIV)*2/(ENE));
    
    // Write results to file
    parout << step_count << "\t" << H << "\t" << T << "\t"
           << Xcu << "\t" << Xmnup << "\t" << Xmndown << "\t" << Xal << "\t"
           << Ycu << "\t" << Ymnup << "\t" << Ymndown << "\t" << Yal << "\t"
           << Zcu << "\t" << Zmnup << "\t" << Zmndown << "\t" << Zal << "\t"
           << Magnetizacion/ENE << "\t"
           << DeltaEAcumM << "\t" << std::endl;
}

/**
 * @brief Saves the final structural configuration to a text file.
 */
void Lattice::saveFinalConfiguration(const char* nombrefile, float Hache, int count) {
    std::ofstream redout;
    if (!OpenFinalRedFile(nombrefile, Hache, count, redout)) {
        return; // Error already reported by helper function
    }
    
    redout.seekp(0, std::ios::beg);
    float aux2, aux3;
    for (int i = 0; i < LATTICE_SIDE; i++) {
        for (int j = 0; j < LATTICE_SIDE; j++) {
            for (int k = 0; k < LATTICE_DEPTH; k++) {
                aux2 = red[i][j][k];
                aux3 = magn[i][j][k];
                // Transformation used in original code for output
                redout << (0.5 * aux2 * aux2 + 1.5 * aux2 - 0.5 * aux3 * aux3 + 1.5 * aux3) << std::endl;
            }
        }
    }
    redout.close();
}


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


/**
 * @brief Executes one full Monte Carlo sweep (N spin flip attempts).
 */
void MonteCarloStep(Lattice& lattice, float T, float H, const SimulationParameters& params, float& DeltaEAcumM) {
    
    for (int XLoc = 0; XLoc < LATTICE_SIDE; XLoc++) {
        for (int YLoc = 0; YLoc < LATTICE_SIDE; YLoc++) {
            for (int ZLoc = 0; ZLoc < LATTICE_DEPTH; ZLoc++) {

                float SpinAct = lattice.getSpin(XLoc, YLoc, ZLoc);
                
                // Calculate neighbor sums for 3rd and 6th NN (as Jm1=Jm2=0)
                float Sum3N = lattice.calculateNeighborSpinSum(XLoc, YLoc, ZLoc, 3);
                float Sum6N = lattice.calculateNeighborSpinSum(XLoc, YLoc, ZLoc, 6);
                
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
                        lattice.flipSpin(XLoc, YLoc, ZLoc);
                        DeltaEAcumM += deltaEM;
                    } else {
                        // Accept: Probabilistically
                        double epsilonM = Ran0a1();
                        float BoltzmannM = exp(-deltaEM / T);
                        if (BoltzmannM >= epsilonM) {
                            lattice.flipSpin(XLoc, YLoc, ZLoc);
                            DeltaEAcumM += deltaEM;
                        }
                    }
                }
            }
        }
    }
}


/**
 * @brief Main simulation loop, iterating over Temperature and Field. (Replaces PasoMC)
 */
void SimulationLoop(const SimulationParameters& params, const char* nombrefile) {
    
    Lattice lattice; 
    std::string init_file = std::string(nombrefile) + ".txt";
    
    // 1. Initialization
    try {
        lattice.loadInitialConfiguration(init_file);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return;
    }

    // 2. Setup Output
    std::ofstream parout;
    if (!OpenLROParametersFile(nombrefile, parout)) {
        return; 
    }
    
    // 3. Monte Carlo Loop Setup
    std::vector<float> listaCampos = createHSweepList(params);
    int H_count = 0; // Counter for final file naming

    std::cout << "J3=" << params.Jm3 << "; J6=" << params.Jm6 << std::endl;

    for (float TEMPERA = params.T_lower; TEMPERA <= params.T_upper; TEMPERA += params.step_T) {
        std::cout << std::endl << "Trabajando a T = " << TEMPERA << std::endl;

        for (float Hache : listaCampos) {
            std::cout << std::endl << "Trabajando a H = " << Hache << std::endl;

            float DeltaEAcumM = 0;

            for (int contador = 1; contador <= params.num_steps; contador++) {
                
                // 3a. Single-Spin Update (Metropolis)
                MonteCarloStep(lattice, TEMPERA, Hache, params, DeltaEAcumM);

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
