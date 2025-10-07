#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>

// Define system dimensions
static const int LATTICE_SIDE = 32;
static const int LATTICE_DEPTH = 64; // 2 * LATTICE_SIDE
static const int LATTICE_TOTAL_SITES = LATTICE_SIDE * LATTICE_SIDE * LATTICE_DEPTH;


/**
 * * @brief Look up table for exp(-ΔE/T) values to optimize Metropolis acceptance checks.
 */
struct BoltzmannTable {
    std::unordered_map<int, double> expTable;
    double T;
};

// Build the lookup table for exp(-ΔE/T)
BoltzmannTable buildBoltzmannTable(double T, double J3, double J6, double H);

// Retrieve precomputed exp(-ΔE/T)
double getBoltzmannFactor(double deltaE, const BoltzmannTable& table);

/**
 * @brief Groups all constant simulation parameters for cleaner function calls.
 */
struct SimulationParameters {
    int num_steps;
    float Jm3; // 3rd NN interaction energy
    float Jm6; // 6th NN interaction energy
    float T_upper;
    float T_lower;
    float step_T;
    float H_upper;
    float H_lower;
    float step_H;
    bool flag_save_config;

    SimulationParameters(int steps, float j3, float j6, float t_up, float t_low, float dt,
                         float h_up, float h_low, float dh, bool flag_red)
        : num_steps(steps), Jm3(j3), Jm6(j6), T_upper(t_up), T_lower(t_low), step_T(dt),
          H_upper(h_up), H_lower(h_low), step_H(dh), flag_save_config(flag_red) {}
};

inline std::ostream& operator<<(std::ostream& os, const SimulationParameters& p) {
    os << "  NUM_STEPS: " << p.num_steps << '\n'
       << "  J_M3: " << p.Jm3 << '\n'
       << "  J_M6: " << p.Jm6 << '\n'
       << "  T_UPPER: " << p.T_upper << '\n'
       << "  T_LOWER: " << p.T_lower << '\n'
       << "  STEP_T: " << p.step_T << '\n'
       << "  H_UPPER: " << p.H_upper << '\n'
       << "  H_LOWER: " << p.H_lower << '\n'
       << "  STEP_H: " << p.step_H << '\n'
       << "  FLAG_SAVE_CONFIG: " << (p.flag_save_config ? "true" : "false");
    return os;
}

/**
 * @brief Encapsulates the lattice data and complex coordinate/neighbor logic.
 */
class Lattice {
private:
    float red[LATTICE_SIDE][LATTICE_SIDE][LATTICE_DEPTH];  // Structural/Species state
    float magn[LATTICE_SIDE][LATTICE_SIDE][LATTICE_DEPTH]; // Magnetic state (Spin)

    // Helper functions for boundary conditions (private utility)
    int wrapCoord(int coord, int max_size) const;
    int wrapCoordZ(int coord) const; // Specific wrapper for 2*side Z coordinate

public:
    // Initialization
    void loadInitialConfiguration(const std::string& filename);
    
    // Core simulation utilities
    float getSpin(int X, int Y, int Z) const { return magn[X][Y][Z]; }
    void flipSpin(int X, int Y, int Z) { magn[X][Y][Z] = -magn[X][Y][Z]; }
    
    // Calculates the required spin sum (used in delta E calculation)
    float calculateNeighborSpinSum(int X, int Y, int Z, int shell_type) const; 
    
    // Output utility
    void saveFinalConfiguration(const char* nombrefile, float Hache, int count);
    
    // Measurement utility (implemented in simulation.cpp)
    void calculateAndWriteLRO(std::ofstream& parout, int step_count, float T, float H, float DeltaEAcumM) const;
};

// Main simulation flow functions
std::vector<float> createHSweepList(const SimulationParameters& params);
void MonteCarloStep(Lattice& lattice, 
                    float H, 
                    const SimulationParameters& params, 
                    const BoltzmannTable& table,
                    float& DeltaEAcumM);
void SimulationLoop(const SimulationParameters& params, const char* nombrefile);

#endif // SIMULATION_H
