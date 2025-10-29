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
struct FastBoltzmannTableSpin {
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
    FastBoltzmannTableSpin(double J3, double J6, double h, double temperature) {
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
 * @struct SiteEnergyTableBEG
 * @brief Pre-calculates the partial Boltzmann factor e^(-E/T) for a single
 * site in a BEG model. Designed for efficient Kawasaki dynamics.
 *
 * This table maps a 5D state (sigma_i, S_Lin_NN, S_Cuad_NN, S_Lin_NNN, S_Cuad_NNN)
 * to a 1D index. The sums *must* exclude the site being swapped with.
 */
struct SiteEnergyTableBEG {
    
    // --- Store Hamiltonian constants ---
    double J1, K1, E1, J2, K2, E2;
    double T;

    // --- Table storing the pre-computed factors e^(-E/T) ---
    vector<double> table;

    // --- State space dimensions for indexing ---
    static constexpr int N_sigma = 3; // sigma_i in {-1, 0, 1}

    // NN sum ranges (over 7 neighbors)
    static constexpr int N_S_Lin_NN = 15;  // {-7, ..., 7}
    static constexpr int N_S_Cuad_NN = 8;   // {0, ..., 7}

    // NNN sum ranges (over 6 neighbors)
    static constexpr int N_S_Lin_NNN = 13; // {-6, ..., 6}
    static constexpr int N_S_Cuad_NNN = 7;  // {0, ..., 6}

    // --- Offsets to map physical values to indices {0, ...} ---
    static constexpr int OFFSET_sigma = 1;     // Maps {-1, 0, 1} to {0, 1, 2}
    static constexpr int OFFSET_S_Lin_NN = 7;  // Maps {-7, ..., 7} to {0, ..., 14}
    static constexpr int OFFSET_S_Cuad_NN = 0; // Maps {0, ..., 7} to {0, ..., 7}
    static constexpr int OFFSET_S_Lin_NNN = 6; // Maps {-6, ..., 6} to {0, ..., 12}
    static constexpr int OFFSET_S_Cuad_NNN = 0; // Maps {0, ..., 6} to {0, ..., 6}

    // --- Pre-calculate multipliers for flat index calculation ---
    // M_SCNN = N_S_Lin_NNN * N_S_Cuad_NNN
    static constexpr size_t M_SCNN = 13 * 7; // 91
    // M_SLNN = N_S_Cuad_NN * M_SCNN
    static constexpr size_t M_SLNN = 8 * M_SCNN; // 728
    // M_SIGMA = N_S_Lin_NN * M_SLNN
    static constexpr size_t M_SIGMA = 15 * M_SLNN; // 10920

private:
    /**
     * @brief Internal helper to calculate site energy.
     */
    double calculateSiteEnergy(double tipo,
                               double sumaLinNN, double sumaCuadNN,
                               double sumaLinNNN, double sumaCuadNNN) {
        
        double tipoSqr = tipo * tipo;

        // NN contribution
        double C1_NN = J1 * tipo + E1 * tipoSqr;
        double C2_NN = K1 * tipoSqr + E1 * tipo;
        double energyNN = C1_NN * sumaLinNN + C2_NN * sumaCuadNN;

        // NNN contribution
        double C1_NNN = J2 * tipo + E2 * tipoSqr;
        double C2_NNN = K2 * tipoSqr + E2 * tipo;
        double energyNNN = C1_NNN * sumaLinNNN + C2_NNN * sumaCuadNNN;

        return energyNN + energyNNN;
    }

public:
    /**
     * @brief Constructs and populates the BEG site energy table.
     */
    SiteEnergyTableBEG(double jota1, double ka1, double ele1,
                       double jota2, double ka2, double ele2,
                       double temperature) 
    {
        // Store parameters
        J1 = jota1; K1 = ka1; E1 = ele1;
        J2 = jota2; K2 = ka2; E2 = ele2;
        T = temperature;
        if (T == 0) T = 1e-9; // Avoid division by zero

        // Resize the table to hold all possible state combinations
        const size_t table_size = N_sigma * N_S_Lin_NN * N_S_Cuad_NN * N_S_Lin_NNN * N_S_Cuad_NNN;
        table.resize(table_size);

        for (size_t index = 0; index < table_size; ++index) {
            
            // --- 1. Deconstruct the flat index back to 5D state indices ---
            const int idx_sigma = index / M_SIGMA;
            size_t temp_idx = index % M_SIGMA;

            const int idx_slnn = temp_idx / M_SLNN;
            temp_idx = temp_idx % M_SLNN;

            const int idx_scnn = temp_idx / M_SCNN;
            temp_idx = temp_idx % M_SCNN;

            const int idx_slnnn = temp_idx / N_S_Cuad_NNN; // or % M_SLNNN
            const int idx_scnnn = temp_idx % N_S_Cuad_NNN;

            // --- 2. Convert zero-based indices back to physical spin/sum values ---
            const int sigma_i = idx_sigma - OFFSET_sigma;
            const int S_Lin_NN = idx_slnn - OFFSET_S_Lin_NN;
            const int S_Cuad_NN = idx_scnn - OFFSET_S_Cuad_NN;
            const int S_Lin_NNN = idx_slnnn - OFFSET_S_Lin_NNN;
            const int S_Cuad_NNN = idx_scnnn - OFFSET_S_Cuad_NNN;

            // --- 3. Calculate energy for this state ---
            double E = calculateSiteEnergy(sigma_i, S_Lin_NN, S_Cuad_NN, S_Lin_NNN, S_Cuad_NNN);

            // --- 4. Store the partial Boltzmann factor ---
            table[index] = exp(-E / T);
        }
    }

    /**
     * @brief Calculates the flat array index from the 5D physical parameters.
     * @return The corresponding index in the flat lookup table.
     */
    inline size_t get_index(int sigma_i, int S_Lin_NN, int S_Cuad_NN, 
                            int S_Lin_NNN, int S_Cuad_NNN) const 
    {
        // 1. Map physical values to zero-based indices
        const size_t idx_sigma = sigma_i + OFFSET_sigma;
        const size_t idx_slnn = S_Lin_NN + OFFSET_S_Lin_NN;
        const size_t idx_scnn = S_Cuad_NN + OFFSET_S_Cuad_NN;
        const size_t idx_slnnn = S_Lin_NNN + OFFSET_S_Lin_NNN;
        const size_t idx_scnnn = S_Cuad_NNN + OFFSET_S_Cuad_NNN;

        // 2. Calculate the flat index (row-major order)
        // (Pre-calculated multipliers are faster than 4 `double` multiplications)
        return (idx_sigma * M_SIGMA) +
               (idx_slnn  * M_SLNN) +
               (idx_scnn  * M_SCNN) +
               (idx_slnnn * N_S_Cuad_NNN) +
               idx_scnnn;
    }

    /**
     * @brief Looks up the partial Boltzmann factor e^(-E/T) for a given state.
     * @return The pre-computed factor.
     */
    inline double lookup(int sigma_i, int S_Lin_NN, int S_Cuad_NN, 
                         int S_Lin_NNN, int S_Cuad_NNN) const 
    {
        return table[get_index(sigma_i, S_Lin_NN, S_Cuad_NN, S_Lin_NNN, S_Cuad_NNN)];
    }
};

/**
* @struct SimulationParameters
* @brief Groups all constant simulation parameters for cleaner function calls.
*/
struct SimulationParameters {
    int num_steps;
    int simulation_method;
    float w1_12; // Interaction energy between species 1 and 2 at 1st NN
    float w2_12; // Interaction energy between species 1 and 2 at 2nd NN
    float w1_13; // Interaction energy between species 1 and 3 at 1st NN
    float w2_13; // Interaction energy between species 1 and 3 at 2nd NN
    float w1_23; // Interaction energy between species 2 and 3 at 1st NN
    float w2_23; // Interaction energy between species 2 and 3 at 2nd NN
    float Jm3; // 3rd NN interaction energy
    float Jm6; // 6th NN interaction energy
    float T_start;
    float T_end;
    float step_T;
    float H_upper;
    float H_lower;
    float step_H;
    bool flag_save_config;

    SimulationParameters(int steps, int sim, float w1_12, float w2_12, float w1_13, 
                         float w2_13, float w1_23, float w2_23, float j3, float j6, 
                         float t_s, float t_e, float dt, float h_up, float h_low, 
                         float dh, bool flag_red)
        : num_steps(steps), simulation_method(sim), w1_12(w1_12), w2_12(w2_12), w1_13(w1_13), 
        w2_13(w2_13), w1_23(w1_23), w2_23(w2_23), Jm3(j3), Jm6(j6), T_start(t_s), T_end(t_e), 
        step_T(dt),H_upper(h_up), H_lower(h_low), step_H(dh), flag_save_config(flag_red) {}
};

inline ostream& operator<<(ostream& os, const SimulationParameters& p) {
    os << "  NUM_STEPS: " << p.num_steps << '\n'
       << "  SIMULATION_METHOD: " << p.simulation_method << '\n'
       << "  W1_12: " << p.w1_12 << '\n'
       << "  W2_12: " << p.w2_12 << '\n'
       << "  W1_13: " << p.w1_13 << '\n'
       << "  W2_13: " << p.w2_13 << '\n'
       << "  W1_23: " << p.w1_23 << '\n'
       << "  W2_23: " << p.w2_23 << '\n'
       << "  J_M3: " << p.Jm3 << '\n'
       << "  J_M6: " << p.Jm6 << '\n'
       << "  T_START: " << p.T_start << '\n'
       << "  T_END: " << p.T_end << '\n'
       << "  STEP_T: " << p.step_T << '\n'
       << "  H_UPPER: " << p.H_upper << '\n'
       << "  H_LOWER: " << p.H_lower << '\n'
       << "  STEP_H: " << p.step_H << '\n'
       << "  FLAG_SAVE_CONFIG: " << (p.flag_save_config ? "true" : "false");
    return os;
}

// Main simulation flow functions
std::vector<float> createHSweepList(const SimulationParameters& params);

std::vector<float> createTSweepList(const SimulationParameters& params);

void MonteCarloStepChemicalExchange(Lattice& lattice,
                                    float H,    
                                    const SimulationParameters& params, 
                                    const SiteEnergyTableBEG& tableBeg,
                                    const FastBoltzmannTableSpin& tableSpin,
                                    float& DeltaEAcumM);

void MonteCarloStepSpinExtH(Lattice& lattice, 
                            float H,
                            const SimulationParameters& params, 
                            const FastBoltzmannTableSpin& table,
                            float& DeltaEAcumM);

void SimulationLoop(const SimulationParameters& params, const char* nombrefile);

#endif // SIMULATION_H
