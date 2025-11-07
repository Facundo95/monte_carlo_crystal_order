#include "lattice.h"
#include "simulation.h"
#include "rng.h"
#include "file_handler.h"
#include <stdexcept>

/** @brief Loads the initial configuration from a specified file.
 * The file should contain the species and spin states for each lattice site.
 * @param filename The path to the input configuration file.
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
        if (count >= m_total_sites) {
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

/** @brief Initializes the neighbor lists for each lattice site. */
void Lattice::initializeNeighbors() {
    for (int site = 0; site < m_total_sites; ++site) {
        int x, y, z;
        idxToXYZ(site, m_side, x, y, z);

        auto& n1 = neighbors1[site];
        auto& n2 = neighbors2[site];
        auto& n3 = neighbors3[site];
        auto& n6 = neighbors6[site];

        // 1st Neighbors (8 sites)
        int p = 0;
        n1[p++] = idx3D(x, y, wrap(z + 1, m_depth));
        n1[p++] = idx3D(x, y, wrap(z - 1, m_depth));
        if (z % 2 == 0) {
            n1[p++] = idx3D(wrap(x + 1, m_side), y, wrap(z + 1, m_depth));
            n1[p++] = idx3D(wrap(x + 1, m_side), y, wrap(z - 1, m_depth));
            n1[p++] = idx3D(x, wrap(y + 1, m_side), wrap(z + 1, m_depth));
            n1[p++] = idx3D(x, wrap(y + 1, m_side), wrap(z - 1, m_depth));
            n1[p++] = idx3D(wrap(x + 1, m_side), wrap(y + 1, m_side), wrap(z + 1, m_depth));
            n1[p++] = idx3D(wrap(x + 1, m_side), wrap(y + 1, m_side), wrap(z - 1, m_depth));
        } else {
            n1[p++] = idx3D(wrap(x - 1, m_side), y, wrap(z + 1, m_depth));
            n1[p++] = idx3D(wrap(x - 1, m_side), y, wrap(z - 1, m_depth));
            n1[p++] = idx3D(x, wrap(y - 1, m_side), wrap(z + 1, m_depth));
            n1[p++] = idx3D(x, wrap(y - 1, m_side), wrap(z - 1, m_depth));
            n1[p++] = idx3D(wrap(x - 1, m_side), wrap(y - 1, m_side), wrap(z + 1, m_depth));
            n1[p++] = idx3D(wrap(x - 1, m_side), wrap(y - 1, m_side), wrap(z - 1, m_depth));
        }

        // 2nd Neighbors (6 sites)
        p = 0;
        n2[p++] = idx3D(wrap(x + 1, m_side), wrap(y, m_side), z);
        n2[p++] = idx3D(wrap(x - 1, m_side), wrap(y, m_side), z);
        n2[p++] = idx3D(x, wrap(y + 1, m_side), z);
        n2[p++] = idx3D(x, wrap(y - 1, m_side), z);
        n2[p++] = idx3D(x, wrap(y, m_side), wrap(z + 2, m_depth));
        n2[p++] = idx3D(x, wrap(y, m_side), wrap(z - 2, m_depth));

        // 3rd Neighbors (12 sites)
        p = 0;
        n3[p++] = idx3D(wrap(x + 1, m_side), wrap(y + 1, m_side), z);
        n3[p++] = idx3D(wrap(x + 1, m_side), wrap(y - 1, m_side), z);
        n3[p++] = idx3D(wrap(x - 1, m_side), wrap(y + 1, m_side), z);
        n3[p++] = idx3D(wrap(x - 1, m_side), wrap(y - 1, m_side), z);
        n3[p++] = idx3D(x, wrap(y + 1, m_side), wrap(z + 2, m_depth));
        n3[p++] = idx3D(x, wrap(y - 1, m_side), wrap(z + 2, m_depth));
        n3[p++] = idx3D(x, wrap(y + 1, m_side), wrap(z - 2, m_depth));
        n3[p++] = idx3D(x, wrap(y - 1, m_side), wrap(z - 2, m_depth));
        n3[p++] = idx3D(wrap(x + 1, m_side), y, wrap(z + 2, m_depth));
        n3[p++] = idx3D(wrap(x - 1, m_side), y, wrap(z + 2, m_depth));
        n3[p++] = idx3D(wrap(x + 1, m_side), y, wrap(z - 2, m_depth));
        n3[p++] = idx3D(wrap(x - 1, m_side), y, wrap(z - 2, m_depth));

        // 6th Neighbors (6 sites)
        p = 0;
        n6[p++] = idx3D(wrap(x + 2, m_side), y, z);
        n6[p++] = idx3D(wrap(x - 2, m_side), y, z);
        n6[p++] = idx3D(x, wrap(y + 2, m_side), z);
        n6[p++] = idx3D(x, wrap(y - 2, m_side), z);
        n6[p++] = idx3D(x, y, wrap(z + 4, m_depth));
        n6[p++] = idx3D(x, y, wrap(z - 4, m_depth));
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

/** @brief Calculates the chemical energy contribution for a given site. 
 * @param type Species type at the site.
 * @param JOTA1, JOTA2 Interaction parameters for 1st and 2nd NN.
 * @param KA1, KA2 Interaction parameters for 1st and 2nd NN.
 * @param ELE1, ELE2 Interaction parameters for 1st and 2nd NN.
 * @param sumLin1, sumCuad1 Linear and quadratic sums over 1st NN species.
 * @param sumLin2, sumCuad2 Linear and quadratic sums over 2nd NN species.
 */
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

/** @brief Calculates the magnetic energy contribution for a given site. 
 * @param Spin Spin state at the site.
 * @param j3, j6 Interaction parameters for 3rd and 6th NN.
 * @param H External magnetic field.
 * @param sum3, sum6 Sums over neighbor spins for 3rd and 6th NN.
 */
float Lattice::calculateSiteMagneticEnergy(int Spin, 
                                        float j3, float j6, 
                                        float H, 
                                        float sum3, float sum6) const {
    return - Spin * (j3 * sum3 + j6 * sum6 + H);
}

/**
 * @brief Calculates and writes the LRO parameters and Magnetization to the output file.
 * @param parout Output file stream.
 * @param step_count Current simulation step count.
 * @param T Current temperature.
 * @param H Current magnetic field.
 * @param DeltaEAcumM Accumulated energy change for magnetization.
 */
void Lattice::calculateAndWriteLRO(std::ofstream& parout, 
                                int step_count, float T, 
                                float H, float DeltaEAcumM) const {
    
    int CuI=0, CuII=0, CuIII=0, CuIV=0;
    int MnUpI=0, MnUpII=0, MnUpIII=0, MnUpIV=0, MnDownI=0, MnDownII=0, MnDownIII=0, MnDownIV=0; 
    int AlI=0, AlII=0, AlIII=0, AlIV=0;

    float Magnetizacion = 0;
    
    // Traversal and counting for LRO parameters
    for (int site = 0; site<m_total_sites;site++){
        int i, j, k;
        idxToXYZ(site, m_side, i, j, k);
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
    const float ENE = (float)m_total_sites;
    
    // X parameters (I+II vs III+IV)
    float Xcu = (CuI+CuII-CuIII-CuIV) / ENE;
    float Xmnup = (MnUpI+MnUpII-MnUpIII-MnUpIV) / ENE;
    float Xmndown = (MnDownI+MnDownII-MnDownIII-MnDownIV) / ENE;
    float Xal = (AlI+AlII-AlIII-AlIV) / ENE;
    
    // Y parameters (I vs II)
    float Ycu = 2 * (CuI-CuII) / ENE;
    float Ymnup = 2 * (MnUpI-MnUpII) / ENE;
    float Ymndown = 2 * (MnDownI-MnDownII) / ENE;
    float Yal = 2 * (AlI-AlII) / ENE;
    
    // Z parameters (III vs IV)
    float Zcu = 2 * (CuIII-CuIV) / ENE;
    float Zmnup = 2 * (MnUpIII-MnUpIV) / ENE;
    float Zmndown = 2 * (MnDownIII-MnDownIV) / ENE;
    float Zal = 2 * (AlIII-AlIV) / ENE;
    
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
 * @param nombrefile Name of the output file.
 * @param Hache Current magnetic field.
 * @param TEMPERA Current temperature.
 * @param count Current simulation step count.
 */
void Lattice::saveFinalConfiguration(const char* nombrefile, 
                                    float Hache, float TEMPERA, int count) {
    std::ofstream redout;
    if (!OpenFinalRedFile(nombrefile, Hache, TEMPERA, count, redout)) {
        return; // Error already reported by helper function
    }
    
    redout.seekp(0, std::ios::beg);
    float aux2, aux3;
    for (int site = 0; site < m_total_sites; site++) {
        aux2 = red_flat[site];
        aux3 = magn_flat[site];
        // Transformation used in original code for output
        redout << (0.5 * aux2 * aux2 + 1.5 * aux2 - 0.5 * aux3 * aux3 + 1.5 * aux3) << std::endl;        
    }
    redout.close();
}