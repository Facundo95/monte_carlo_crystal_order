#include "lattice.h"
#include "simulation.h"
#include "rng.h"
#include "file_handler.h"
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <cmath>

/** @brief Loads the initial configuration from a specified file.
 * The file should contain the species and spin states for each lattice site.
 * @param filename The path to the input configuration file.
 */
void Lattice::loadInitialConfiguration(const std::string& filename) {
    // Ensure the filename has a .txt extension (case-insensitive)
    auto pos = filename.find_last_of('.');
    std::string ext;
    if (pos == std::string::npos) ext = "";
    else ext = filename.substr(pos);
    // lower-case ext for case-insensitive compare
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    if (ext == ".txt") {
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
        if (count != m_total_sites) {
            throw std::runtime_error("El archivo de entrada no tiene el numero esperado de sitios.");
        }
        redin.close();
    } else if (ext == ".xyz") {
        std::ifstream redin(filename);
        if (!redin.is_open()) {
            throw std::runtime_error("No se pudo abrir el archivo de entrada inicial: " + filename);
        }

        std::string line;
        // 1) read number of atoms (first line)
        if (!std::getline(redin, line)) {
            throw std::runtime_error("Archivo .xyz vacio: " + filename);
        }
        // optionally parse and check the number
        int natoms = 0;
        try {
            natoms = std::stoi(line);
        } catch (...) {
            // ignore parse errors, we'll validate by counting lines
            natoms = -1;
        }

        // 2) read comment/header line
        if (!std::getline(redin, line)) {
            throw std::runtime_error("Archivo .xyz mal formado (sin linea de comentario): " + filename);
        }

        int filled = 0;
        while (std::getline(redin, line)) {
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::string atomToken;
            double x,y,z;
            int spin;
            if (!(ss >> atomToken >> x >> y >> z >> spin)) {
                throw std::runtime_error("Linea .xyz mal formada (se esperaba: atom x y z spin): " + line);
            }

            // atomToken must be an element symbol (no atomic number token allowed)
            std::string tok = atomToken;
            std::transform(tok.begin(), tok.end(), tok.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
            int specie = 0;
            if (tok == "co") specie = 1;
            else if (tok == "ni") specie = 0;
            else if (tok == "al") specie = -1;
            else {
                throw std::runtime_error("Tipo de atomo no reconocido en .xyz (solo símbolos permitidos): " + atomToken);
            }

            // Map coordinates back to integer lattice indices
            int ix = static_cast<int>(std::round(x));
            int iy = static_cast<int>(std::round(y));
            int kz = static_cast<int>(std::round(z / 0.5));

            if (ix < 0 || ix >= m_side || iy < 0 || iy >= m_side || kz < 0 || kz >= m_depth) {
                throw std::runtime_error("Coordenadas fuera de rango en .xyz: " + line);
            }

            int idx = idx3D(ix, iy, kz);
            red_flat[idx] = specie;
            magn_flat[idx] = spin;
            ++filled;
        }

        if (natoms != -1 && natoms != filled) {
            throw std::runtime_error("El numero de atomos en cabecera no coincide con las lineas leidas en: " + filename);
        }
        if (filled != m_total_sites) {
            throw std::runtime_error("El archivo .xyz no contiene el numero esperado de sitios.");
        }

        redin.close();
    } else {
        throw std::runtime_error("Formato de archivo no soportado para el archivo de entrada inicial: " + filename);
    }
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
double Lattice::calculateDeltaChemicalEnergy(int type_A, int type_N, 
                                            float JOTA1, float JOTA2, 
                                            float KA1, float KA2, 
                                            float ELE1, float ELE2,
                                            int sumLinNN_A, int sumLinNNN_A,
                                            int sumCuadNN_A, int sumCuadNNN_A,
                                            int sumLinNN_N, int sumLinNNN_N,
                                            int sumCuadNN_N, int sumCuadNNN_N) const {
    // 1. Pre-calculate state differences
    double diffTipo = type_N - type_A;
    double sumTipo  = type_N + type_A;

    // 2. Pre-calculate Sum differences (Nearest Neighbors)
    double diffLinNN  = sumLinNN_A  - sumLinNN_N;
    double diffCuadNN = sumCuadNN_A - sumCuadNN_N;

    // 3. Pre-calculate Sum differences (Next Nearest Neighbors)
    double diffLinNNN  = sumLinNNN_A  - sumLinNNN_N;
    double diffCuadNNN = sumCuadNNN_A - sumCuadNNN_N;

    // 4. Calculate Delta directly
    // Notice that 'diffTipo' factors out of the entire equation
    float deltaEQ = diffTipo * (
        // Nearest Neighbors (NN)
        JOTA1 * diffLinNN +
        KA1   * diffCuadNN * sumTipo +
        ELE1  * (diffLinNN * sumTipo + diffCuadNN) +

        // Next Nearest Neighbors (NNN)
        JOTA2 * diffLinNNN +
        KA2   * diffCuadNNN * sumTipo +
        ELE2  * (diffLinNNN * sumTipo + diffCuadNNN));
    
    return deltaEQ;
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
    
    int AI=0, AII=0, AIII=0, AIV=0;
    int BUpI=0, BUpII=0, BUpIII=0, BUpIV=0, BDownI=0, BDownII=0, BDownIII=0, BDownIV=0; 
    int CI=0, CII=0, CIII=0, CIV=0;

    float Magnetizacion = 0;
    
    // Traversal and counting for LRO parameters
    for (int site = 0; site<m_total_sites;site++){
        int i, j, k;
        idxToXYZ(site, m_side, i, j, k);
        Magnetizacion += magn_flat[site];

        // Determine sublattice (I, II, III, or IV)
        bool Z_is_even = (k % 2 == 0);
        bool XY_sum_is_even = ((i + j + (k/2)) % 2 == 0);

        int* A_ptr = nullptr;
        int* BUp_ptr = nullptr;
        int* BDown_ptr = nullptr;
        int* C_ptr = nullptr;
        
        // Simplified conditional logic using pointers
        if (Z_is_even) {
            if (XY_sum_is_even) { // Sublattice I (Z even, X+Y+Z/2 even)
                A_ptr = &AI; BUp_ptr = &BUpI; BDown_ptr = &BDownI; C_ptr = &CI;
            } else { // Sublattice II (Z even, X+Y+Z/2 odd)
                A_ptr = &AII; BUp_ptr = &BUpII; BDown_ptr = &BDownII; C_ptr = &CII;
            }
        } else {
            if (XY_sum_is_even) { // Sublattice III (Z odd, X+Y+Z/2 even)
                A_ptr = &AIII; BUp_ptr = &BUpIII; BDown_ptr = &BDownIII; C_ptr = &CIII;
            } else { // Sublattice IV (Z odd, X+Y+Z/2 odd)
                A_ptr = &AIV; BUp_ptr = &BUpIV; BDown_ptr = &BDownIV; C_ptr = &CIV;
            }
        }
        
        // Increment counters based on species (red[site] = 1, 0, or -1)
        if (red_flat[site] == 1) { // Copper
            (*A_ptr)++;
        } else if (red_flat[site] == 0) { // Manganese
            if (magn_flat[site] == 1) {
                (*BUp_ptr)++;
            } else {
                (*BDown_ptr)++;
            }
        } else if (red_flat[site] == -1) { // Aluminum
            (*C_ptr)++;
        }
    }
    
    // Calculate LRO Parameters (X, Y, Z) and normalized Magnetization
    const float ENE = (float)m_total_sites;
    
    // X parameters (I+II vs III+IV)
    float X_A = (AI+AII-AIII-AIV) / ENE;
    float X_Bup = (BUpI+BUpII-BUpIII-BUpIV) / ENE;
    float X_Bdown = (BDownI+BDownII-BDownIII-BDownIV) / ENE;
    float X_C = (CI+CII-CIII-CIV) / ENE;
    
    // Y parameters (I vs II)
    float Y_A = 2 * (AI-AII) / ENE;
    float Y_Bup = 2 * (BUpI-BUpII) / ENE;
    float Y_Bdown = 2 * (BDownI-BDownII) / ENE;
    float Y_C = 2 * (CI-CII) / ENE;
    
    // Z parameters (III vs IV)
    float Z_A = 2 * (AIII-AIV) / ENE;
    float Z_Bup = 2 * (BUpIII-BUpIV) / ENE;
    float Z_Bdown = 2 * (BDownIII-BDownIV) / ENE;
    float Z_C = 2 * (CIII-CIV) / ENE;
    
    // Write results to file
    parout << step_count << "\t" << H << "\t" << T << "\t"
           << X_A << "\t" << X_Bup << "\t" << X_Bdown << "\t" << X_C << "\t"
           << Y_A << "\t" << Y_Bup << "\t" << Y_Bdown << "\t" << Y_C << "\t"
           << Z_A << "\t" << Z_Bup << "\t" << Z_Bdown << "\t" << Z_C << "\t"
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
bool Lattice::saveFinalConfiguration(const char* nombrefile, 
                                    float Hache, float TEMPERA, int count) {
    std::ofstream redout;
    if (!OpenFinalRedFile(nombrefile, Hache, TEMPERA, count, redout)) {
        return false;
    }

    std::string filename(nombrefile);
    auto pos = filename.find_last_of('.');
    std::string ext;
    if (pos == std::string::npos) ext = "";
    else ext = filename.substr(pos);
    // lower-case ext for case-insensitive compare
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });

    if (ext == ".xyz") {
        // Write extended .xyz header
        redout << m_total_sites << "\n";
        redout << "Lattice=\"" << m_side << " 0.0 0.0 0.0 " << m_side << " 0.0 0.0 0.0" << m_depth << "\"\n";
        redout << "Properties=\"species:S:1:pos:R:3:spin:I:1\"" << "\n";

        // For each site, write: Element x y z species spin
        for (int site = 0; site < m_total_sites; ++site) {
            int i, j, k;
            idxToXYZ(site, m_side, i, j, k);
            // Simple geometric mapping: unit spacing in x,y and 0.5 in z to reflect depth
            float x;
            float y;
            float z = static_cast<float>(k) * 0.5f;

            if ((k % 2) != 0) {
                float x = static_cast<float>(i) + 0.5f;
                float y = static_cast<float>(j) + 0.5f;
            } else {
                float x = static_cast<float>(i);
                float y = static_cast<float>(j);
            }

            int specie = red_flat[site];
            int spin = magn_flat[site];
            const char* elem = (specie == 1) ? "Co" : (specie == 0 ? "Ni" : "Al");

            // Write: ElementSymbol x y z specie spin
            redout << elem << " " << x << " " << y << " " << z << " " << spin << "\n";
        }

        redout.close();
        return true;
    } else {    
        float aux2, aux3;
        for (int site = 0; site < m_total_sites; site++) {
            aux2 = red_flat[site];
            aux3 = magn_flat[site];
            // Transformation used in original code for output
            redout << (0.5 * aux2 * aux2 + 1.5 * aux2 - 0.5 * aux3 * aux3 + 1.5 * aux3) << std::endl;        
        }
        redout.close();
        return true;
    }
    return false;
}