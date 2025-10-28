#include "lattice.h"
#include "simulation.h"
#include "rng.h"
#include "file_handler.h"
#include <stdexcept>

// --- Utility Functions for Lattice ---

/**
 * @brief Loads the initial configuration from a .txt file.
 */
void Lattice::loadInitialConfiguration(const std::string& filename) {
    std::ifstream redin(filename, std::ios::in | std::ios::binary);
    if (!redin.is_open()) {
        throw std::runtime_error("No se pudo abrir el archivo de entrada inicial: " + filename);
    }

    int aux;
    int count=0;
    redin.seekg(0, std::ios::beg);
    while (redin >> aux) {
        if (count >= LATTICE_TOTAL_SITES) {
            throw std::runtime_error("El archivo de entrada tiene más datos de los esperados.");
        }
        // Polynomial transformation from the original code
        float temp_magn = (-1./12)*aux*aux*aux-(1./3)*aux*aux+(7./12)*aux+(5./6);
        float temp_red = (-1./12)*aux*aux*aux+(1./3)*aux*aux+(7./12)*aux-(5./6);
        magn_flat[count] = static_cast<int>(temp_magn);
        red_flat[count] = static_cast<int>(temp_red);
        count++;
    }
    redin.close();
}

void Lattice::initializeNeighbors() {
    for (int site = 0; site < LATTICE_TOTAL_SITES; ++site) {
        int x, y, z;
        idxToXYZ(site, LATTICE_SIDE, LATTICE_DEPTH, x, y, z);

        std::array<int, 8> n1;
        std::array<int, 6> n2;
        std::array<int, 12> n3;
        std::array<int, 6> n6;

        // 1st Neighbors (8 sites)
        int p = 0;
        n1[p++] = idx3D(x, y, wrap(z + 1, LATTICE_DEPTH));
        n1[p++] = idx3D(x, y, wrap(z - 1, LATTICE_DEPTH));
        if (z % 2 == 0) {
            n1[p++] = idx3D(wrap(x + 1, LATTICE_SIDE), y, wrap(z + 1, LATTICE_DEPTH));
            n1[p++] = idx3D(wrap(x + 1, LATTICE_SIDE), y, wrap(z - 1, LATTICE_DEPTH));
            n1[p++] = idx3D(x, wrap(y + 1, LATTICE_SIDE), wrap(z + 1, LATTICE_DEPTH));
            n1[p++] = idx3D(x, wrap(y + 1, LATTICE_SIDE), wrap(z - 1, LATTICE_DEPTH));
            n1[p++] = idx3D(wrap(x + 1, LATTICE_SIDE), wrap(y + 1, LATTICE_SIDE), wrap(z + 1, LATTICE_DEPTH));
            n1[p++] = idx3D(wrap(x + 1, LATTICE_SIDE), wrap(y + 1, LATTICE_SIDE), wrap(z - 1, LATTICE_DEPTH));
        } else {
            n1[p++] = idx3D(wrap(x - 1, LATTICE_SIDE), y, wrap(z + 1, LATTICE_DEPTH));
            n1[p++] = idx3D(wrap(x - 1, LATTICE_SIDE), y, wrap(z - 1, LATTICE_DEPTH));
            n1[p++] = idx3D(x, wrap(y - 1, LATTICE_SIDE), wrap(z + 1, LATTICE_DEPTH));
            n1[p++] = idx3D(x, wrap(y - 1, LATTICE_SIDE), wrap(z - 1, LATTICE_DEPTH));
            n1[p++] = idx3D(wrap(x - 1, LATTICE_SIDE), wrap(y - 1, LATTICE_SIDE), wrap(z + 1, LATTICE_DEPTH));
            n1[p++] = idx3D(wrap(x - 1, LATTICE_SIDE), wrap(y - 1, LATTICE_SIDE), wrap(z - 1, LATTICE_DEPTH));
        }
        neighbors1[site] = n1;

        // 2nd Neighbors (6 sites)
        p = 0;
        n2[p++] = idx3D(wrap(x + 1, LATTICE_SIDE), wrap(y, LATTICE_SIDE), z);
        n2[p++] = idx3D(wrap(x - 1, LATTICE_SIDE), wrap(y, LATTICE_SIDE), z);
        n2[p++] = idx3D(wrap(x, LATTICE_SIDE), wrap(y + 1, LATTICE_SIDE), z);
        n2[p++] = idx3D(wrap(x, LATTICE_SIDE), wrap(y - 1, LATTICE_SIDE), z);
        n2[p++] = idx3D(wrap(x, LATTICE_SIDE), wrap(y, LATTICE_SIDE), wrap(z + 2, LATTICE_DEPTH));
        n2[p++] = idx3D(wrap(x, LATTICE_SIDE), wrap(y, LATTICE_SIDE), wrap(z - 2, LATTICE_DEPTH));
        neighbors2[site] = n2;

        // 3rd Neighbors (12 sites)
        int p = 0;
        n3[p++] = idx3D(wrap(x + 1, LATTICE_SIDE), wrap(y + 1, LATTICE_SIDE), z);
        n3[p++] = idx3D(wrap(x + 1, LATTICE_SIDE), wrap(y - 1, LATTICE_SIDE), z);
        n3[p++] = idx3D(wrap(x - 1, LATTICE_SIDE), wrap(y + 1, LATTICE_SIDE), z);
        n3[p++] = idx3D(wrap(x - 1, LATTICE_SIDE), wrap(y - 1, LATTICE_SIDE), z);
        n3[p++] = idx3D(x, wrap(y + 1, LATTICE_SIDE), wrap(z + 2, LATTICE_DEPTH));
        n3[p++] = idx3D(x, wrap(y - 1, LATTICE_SIDE), wrap(z + 2, LATTICE_DEPTH));
        n3[p++] = idx3D(x, wrap(y + 1, LATTICE_SIDE), wrap(z - 2, LATTICE_DEPTH));
        n3[p++] = idx3D(x, wrap(y - 1, LATTICE_SIDE), wrap(z - 2, LATTICE_DEPTH));
        n3[p++] = idx3D(wrap(x + 1, LATTICE_SIDE), y, wrap(z + 2, LATTICE_DEPTH));
        n3[p++] = idx3D(wrap(x - 1, LATTICE_SIDE), y, wrap(z + 2, LATTICE_DEPTH));
        n3[p++] = idx3D(wrap(x + 1, LATTICE_SIDE), y, wrap(z - 2, LATTICE_DEPTH));
        n3[p++] = idx3D(wrap(x - 1, LATTICE_SIDE), y, wrap(z - 2, LATTICE_DEPTH));
        neighbors3[site] = n3;

        // 6th Neighbors (6 sites)
        p = 0;
        n6[p++] = idx3D(wrap(x + 2, LATTICE_SIDE), y, z);
        n6[p++] = idx3D(wrap(x - 2, LATTICE_SIDE), y, z);
        n6[p++] = idx3D(x, wrap(y + 2, LATTICE_SIDE), z);
        n6[p++] = idx3D(x, wrap(y - 2, LATTICE_SIDE), z);
        n6[p++] = idx3D(x, y, wrap(z + 4, LATTICE_DEPTH));
        n6[p++] = idx3D(x, y, wrap(z - 4, LATTICE_DEPTH));
        neighbors6[site] = n6;
    }
}


/**
 * @brief Calculates the sum of spins for a specific neighbor shell (3rd or 6th NN).
 * @param shell_type Must be 3 (3rd NN) or 6 (6th NN).
 * (Original Jm1 and Jm2 (1st and 2nd NN) were 0, so only 3 and 6 are implemented here).
 */
float Lattice::calculateNeighborSpinSum(int site, int shell_type) const {
    if (shell_type != 1 && shell_type != 2 && shell_type != 3 && shell_type != 6) {
        throw std::invalid_argument("shell_type options: 1,2,3 or 6.");
    }
    float sum = 0;
    if (shell_type == 1) {
        const auto &n1 = neighbors1[site];
        for (int i = 0; i < 8; ++i) sum += magn_flat[n1[i]];
    } else if (shell_type == 2) {
        const auto &n2 = neighbors2[site];
        for (int i = 0; i < 6; ++i) sum += magn_flat[n2[i]];
    } else if (shell_type == 3) {
        const auto &n3 = neighbors3[site];
        for (int i = 0; i < 12; ++i) sum += magn_flat[n3[i]];
    } else if (shell_type == 6) {
        const auto &n6 = neighbors6[site];
        for (int i = 0; i < 6; ++i) sum += magn_flat[n6[i]];
    }   
    return sum;
}

/**
 * @brief Calculates the sum of species (red_flat) for a specific neighbor shell.
 * @param shell_type Must be 1, 2, 3, or 6.
 * @param order Must be 1 (linear) or 2 (quadratic).
 */
float Lattice::calculateNeighborSpeciesSum(int site, int shell_type, int order) const {
    if (shell_type != 1 && shell_type != 2 && shell_type != 3 && shell_type != 6) {
        throw std::invalid_argument("shell_type options: 1,2,3 or 6.");
    }
    if (order != 1 && order != 2) {
        throw std::invalid_argument("order options: 1 (linear) or 2 (quadratic).");
    }
    float sum = 0;
    if (shell_type == 1) {
        const auto &n1 = neighbors1[site];
        for (int i = 0; i < 8; ++i) {
            float val = static_cast<float>(red_flat[n1[i]]);
            sum += (order == 1) ? val : val * val;
        }
    } else if (shell_type == 2) {
        const auto &n2 = neighbors2[site];
        for (int i = 0; i < 6; ++i) {
            float val = static_cast<float>(red_flat[n2[i]]);
            sum += (order == 1) ? val : val * val;
        }
    } else if (shell_type == 3) {
        const auto &n3 = neighbors3[site];
        for (int i = 0; i < 12; ++i) {
            float val = static_cast<float>(red_flat[n3[i]]);
            sum += (order == 1) ? val : val * val;
        }
    } else if (shell_type == 6) {
        const auto &n6 = neighbors6[site];
        for (int i = 0; i < 6; ++i) {
            float val = static_cast<float>(red_flat[n6[i]]);
            sum += (order == 1) ? val : val * val;
        }
    }   
    return sum;
}

float Lattice::calculateSiteChemicalEnergy(int type, 
                                        float JOTA1, float JOTA2, 
                                        float KA1, float KA2, 
                                        float ELE1, float ELE2,
                                        float sumLin1, float sumCuad1,
                                        float sumLin2, float sumCuad2) const {
        float typeSqr = type * type;

        // NN contribution
        float C1_NN = JOTA1 * type + ELE1 * typeSqr;
        float C2_NN = KA1 * typeSqr + ELE1 * type;
        float energyNN = C1_NN * sumLin1 + C2_NN * sumCuad1;

        // NNN contribution
        float C1_NNN = JOTA2 * type + ELE2 * typeSqr;
        float C2_NNN = KA2 * typeSqr + ELE2 * type;
        float energyNNN = C1_NNN * sumLin2 + C2_NNN * sumCuad2;

        return energyNN + energyNNN;
}

float Lattice::calculateSiteMagneticEnergy(int Spin, float j3, float j6, float H, float sum3, float sum6) const {
    float dEM = - Spin * (j3 * sum3 + j6 * sum6) - H * Spin;
    return dEM;
}

/**
 * @brief Calculates and writes the LRO parameters and Magnetization to the output file.
 */
void Lattice::calculateAndWriteLRO(std::ofstream& parout, int step_count, float T, float H, float DeltaEAcumM) const {
    
    int CuI=0, CuII=0, CuIII=0, CuIV=0;
    int MnUpI=0, MnUpII=0, MnUpIII=0, MnUpIV=0, MnDownI=0, MnDownII=0, MnDownIII=0, MnDownIV=0; 
    int AlI=0, AlII=0, AlIII=0, AlIV=0;

    float Magnetizacion = 0;
    
    // Traversal and counting for LRO parameters
    for (int site = 0; site<LATTICE_TOTAL_SITES;site++){
        int i, j, k;
        idxToXYZ(site, LATTICE_SIDE, LATTICE_DEPTH, i, j, k);
        Magnetizacion += magn_flat[site];

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
        
        // Increment counters based on species (red[site] = 1, 0, or -1)
        if (red_flat[site] == 1) { // Copper
            (*Cu_ptr)++;
        } else if (red_flat[site] == 0) { // Manganese
            if (magn_flat[site] == 1) {
                (*MnUp_ptr)++;
            } else {
                (*MnDown_ptr)++;
            }
        } else if (red_flat[site] == -1) { // Aluminum
            (*Al_ptr)++;
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
void Lattice::saveFinalConfiguration(const char* nombrefile, float Hache, float TEMPERA, int count) {
    std::ofstream redout;
    if (!OpenFinalRedFile(nombrefile, Hache, count, redout)) {
        return; // Error already reported by helper function
    }
    
    redout.seekp(0, std::ios::beg);
    float aux2, aux3;
    for (int site = 0; site < LATTICE_SIDE; site++) {
        aux2 = red_flat[site];
        aux3 = magn_flat[site];
        // Transformation used in original code for output
        redout << (0.5 * aux2 * aux2 + 1.5 * aux2 - 0.5 * aux3 * aux3 + 1.5 * aux3) << std::endl;        
    }
    redout.close();
}