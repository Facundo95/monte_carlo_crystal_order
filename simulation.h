#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <array>
#include <cstdint> // For int8_t
#include "lattice.h"

using namespace std;

/**
 * @struct BoltzmannTable
 * @brief A structure to pre-calculate and store Boltzmann factors using direct indexing.
 *
 * This table uses a flat std::vector where each possible state change (defined
 * by the spin being flipped and the state of its neighbors) maps to a unique index.
 * This is highly efficient for systems with a discrete and finite number of
 * possible energy changes.
 */
struct FastBoltzmannTable {
    // The flat vector storing the pre-computed Boltzmann factors.
    vector<float> table;

    // Dimensions for the 3D to 1D index mapping
    const int N_si = 2;   // s_i can be {-1, 1}
    const int N_S3 = 13;  // Sum of 12 neighbors can be {-12, -10, ..., 12}
    const int N_S6 = 7;   // Sum of 6 neighbors can be {-6, -4, ..., 6}

    /**
     * @brief Constructs and populates the BoltzmannTable for the specific BCC model.
     * @param J3 Coupling constant for 3rd neighbors.
     * @param J6 Coupling constant for 6th neighbors.
     * @param h External magnetic field.
     * @param temperature The temperature T (or kT) of the system.
     */
    FastBoltzmannTable(double J3, double J6, double h, double temperature) {
        // Resize the table to hold all possible combinations
        const int table_size = N_si * N_S3 * N_S6; // 2 * 13 * 7 = 182
        table.resize(table_size);

        for (int index = 0; index < table_size; ++index) {
            // Perform the inverse mapping from the flat index back to the physical state (s_i, S3, S6)
            
            // 1. Deconstruct the index
            const int term_S3_S6 = N_S3 * N_S6;
            const int idx_si = index / term_S3_S6;
            const int temp_index = index - idx_si * term_S3_S6;
            const int idx_S3 = temp_index / N_S6;
            const int idx_S6 = temp_index - idx_S3 * N_S6;

            // 2. Convert zero-based indices back to physical spin values
            const int s_i_val = (idx_si * 2) - 1;      // Maps {0, 1} to {-1, 1}
            const int S3_val = (idx_S3 * 2) - 12;     // Maps {0, 1,..., 12} to {-12, -10,...}
            const int S6_val = (idx_S6 * 2) - 6;      // Maps {0, 1,..., 6} to {-6, -4,...}

            // 3. Calculate the energy change for this specific combination
            double dE = -2.0 * s_i_val * (J3 * S3_val + J6 * S6_val + h);

            // 4. Store the pre-computed Boltzmann factor at the current index
            table[index] = exp(-dE / temperature);
        }
    }

    /**
     * @brief Calculates the flat array index from the physical parameters.
     * This function maps a 3D state (s_i, S3, S6) to a 1D index.
     * @param s_i The spin being flipped (-1 or +1).
     * @param S3 The sum of the 12 3rd-nearest neighbors.
     * @param S6 The sum of the 6 6th-nearest neighbors.
     * @return The corresponding index in the flat lookup table.
     */
    inline int get_index(int s_i, int S3, int S6) const {
        // Map physical values to zero-based indices
        const int idx_si = (s_i + 1) / 2;    // Maps {-1, 1} to {0, 1}
        const int idx_S3 = (S3 + 12) / 2;   // Maps {-12, -10,...} to {0, 1,..., 12}
        const int idx_S6 = (S6 + 6) / 2;    // Maps {-6, -4,...} to {0, 1,..., 6}
        
        // Calculate the flat index using row-major order logic
        return idx_si * (N_S3 * N_S6) + idx_S3 * N_S6 + idx_S6;
    }

    /**
     * @brief Looks up the Boltzmann factor for a given state change.
     * @param s_i The spin being flipped (-1 or +1).
     * @param S3 The sum of the 12 3rd-nearest neighbors.
     * @param S6 The sum of the 6 6th-nearest neighbors.
     * @return The pre-computed Boltzmann factor.
     */
    double lookup(int s_i, int S3, int S6) const {
        return table[get_index(s_i, S3, S6)];
    }
};

/**
* @struct SimulationParameters
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

inline ostream& operator<<(ostream& os, const SimulationParameters& p) {
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

// Main simulation flow functions
std::vector<float> createHSweepList(const SimulationParameters& params);
void MonteCarloStepSpinExtH(Lattice& lattice, 
                            float H, 
                            const SimulationParameters& params, 
                            const FastBoltzmannTable& table,
                            float& DeltaEAcumM);
void SimulationLoop(const SimulationParameters& params, const char* nombrefile);

#endif // SIMULATION_H
