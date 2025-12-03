#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>
#include <cstring>
#include <cstdint> // For int8_t

#include "lattice.h"

#include <unordered_map>
#include <vector>
#include <cmath>
#include <cstring>

struct BoltzmannDeltaETable {
    std::unordered_map<std::uint64_t, double> table;

    // --- codificar un double como entero de 64 bits (bitwise exact) ---
    static inline std::uint64_t encode(double x) {
        std::uint64_t bits;
        static_assert(sizeof(double) == sizeof(std::uint64_t),
                      "double must be 64-bit");
        std::memcpy(&bits, &x, sizeof(double));
        return bits;
    }

    BoltzmannDeltaETable(const std::vector<double>& all_dEs, double T) {
        if (T == 0.0) T = 1e-12;

        for (double dE : all_dEs) {
            std::uint64_t key = encode(dE);
            table[key] = std::exp(-dE / T);
        }
    }

    inline double lookup(double dE) const {
        std::uint64_t key = encode(dE);
        auto it = table.find(key);
        if (it != table.end()) return it->second;

        // --- si no se encontró: robust fallback ---
        return 0.0;  // o 1.0 dependiendo del criterio que prefieras
    }
};

struct FastBoltzmannTableSpin {

    BoltzmannDeltaETable table;

    double J3_, J6_, h_, T_;

    FastBoltzmannTableSpin(double J3, double J6, double h, double T)
        : table(compute_dEs(J3, J6, h), T), J3_(J3), J6_(J6), h_(h), T_(T)
    {}

    static std::vector<double> compute_dEs(double J3, double J6, double h) {

        static constexpr int N_si = 2;
        static constexpr int N_S3 = 26;
        static constexpr int N_S6 = 14;

        std::vector<double> dEs;
        dEs.reserve(N_si * N_S3 * N_S6);

        int s_vals[2] = {-1, 1};

        int S3_vals[N_S3];
        for (int i = 0; i < N_S3; i++)
            S3_vals[i] = -12 + i;

        int S6_vals[N_S6];
        for (int i = 0; i < N_S6; i++)
            S6_vals[i] = -6 + i;

        for (int si : s_vals)
            for (int S3 : S3_vals)
                for (int S6 : S6_vals)
                    dEs.push_back( 2.0 * si * (J3*S3 + J6*S6 + h) );

        return dEs;
    }

    inline double lookup_from_dE(double dE) const {
        return table.lookup(dE);
    }
};

struct FastBoltzmannTableBEG {

    BoltzmannDeltaETable table;

    float w1_12, w1_13, w1_23;
    float w2_12, w2_13, w2_23;
    float T;

    FastBoltzmannTableBEG(float w1_12_, float w1_13_, float w1_23_,
                         float w2_12_, float w2_13_, float w2_23_,
                         float T_)
        : table(compute_dEs(w1_12_, w1_13_, w1_23_,
                            w2_12_, w2_13_, w2_23_), T_), 
          w1_12(w1_12_), w1_13(w1_13_), w1_23(w1_23_),
          w2_12(w2_12_), w2_13(w2_13_), w2_23(w2_23_), T(T_)
    {}

    static std::vector<double> compute_dEs(float w1_12, float w1_13, float w1_23,
                                           float w2_12, float w2_13, float w2_23) {
        
        float jota1 = 0.25 * w1_13;
        float jota2 = 0.25 * w2_13;
        float ka1 = 0.25 * (2 * w1_12 + 2 * w1_23 - w1_13);
        float ka2 = 0.25 * (2 * w2_12 + 2 * w2_23 - w2_13);
        float ele1 = 0.25 * (w1_12 - w1_23);
        float ele2 = 0.25 * (w2_12 - w2_23);

        std::vector<double> dEs;
        
        int diffLinNN[33] = {};
        int diffCuadNN[16] = {};
        int diffLinNNN[25] = {};
        int diffCuadNNN[13] = {};
        int diffTipo[4] = {-2, -1, 1, 2};
        int sumTipo[3] = {-1, 0, 1};
        for (int i = 0; i < 33; ++i) diffLinNN[i] = i - 16;
        for (int i = 0; i < 16; ++i) diffCuadNN[i] = i;
        for (int i = 0; i < 25; ++i) diffLinNNN[i] = i - 12;
        for (int i = 0; i < 13; ++i) diffCuadNNN[i] = i;

        for (int diff_tipo : diffTipo) {
            for (int sum_tipo : sumTipo) {
                for (int dLinNN : diffLinNN) {
                    for (int dCuadNN : diffCuadNN) {
                        for (int dLinNNN : diffLinNNN) {
                            for (int dCuadNNN : diffCuadNNN) {
                                double dE = diff_tipo * (
                                    // Nearest Neighbors (NN)
                                    jota1 * dLinNN +
                                    ka1   * dCuadNN * sum_tipo +
                                    ele1  * (dLinNN * sum_tipo + dCuadNN) + 
                                    // Next Nearest Neighbors (NNN)
                                    jota2 * dLinNNN +
                                    ka2   * dCuadNNN * sum_tipo +
                                    ele2  * (dLinNNN * sum_tipo + dCuadNNN)
                                );
                                
                                dEs.push_back(dE);
                            }
                        }
                    }
                }
            }
        }
        return dEs;
    }

    inline double look_from_dE(double dE) const {
        return table.lookup(dE);
    }
};

/**
* @struct SimulationParameters
* @brief Groups all constant simulation parameters for cleaner function calls.
*/
struct SimulationParameters {
    int num_steps;
    int simulation_method;
    int lattice_side;
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
    float H_start;
    float H_end;
    float step_H;
    int steps_to_output;
    bool flag_save_config;
    bool flag_loop;

    SimulationParameters(int steps, int sim, int side, float w1_12, float w2_12, float w1_13, 
                         float w2_13, float w1_23, float w2_23, float j3, float j6, 
                         float t_s, float t_e, float dt, float h_start, float h_end, 
                         float dh, int step_out, bool flag_red, bool loop)
        : num_steps(steps), simulation_method(sim), lattice_side(side), 
        w1_12(w1_12), w2_12(w2_12), w1_13(w1_13), w2_13(w2_13), w1_23(w1_23), w2_23(w2_23), 
        Jm3(j3), Jm6(j6), T_start(t_s), T_end(t_e), step_T(dt),H_start(h_start), H_end(h_end), 
        step_H(dh), steps_to_output(step_out), flag_save_config(flag_red), flag_loop(loop) {}
};

inline std::ostream& operator<<(std::ostream& os, const SimulationParameters& p) {
    os << "  NUM_STEPS: " << p.num_steps << '\n'
       << "  SIMULATION_METHOD: " << p.simulation_method << '\n'
       << "  LATTICE_SIDE: " << p.lattice_side << '\n'
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
       << "  H_UPPER: " << p.H_start << '\n'
       << "  H_LOWER: " << p.H_end << '\n'
       << "  STEP_H: " << p.step_H << '\n'
       << "  LOOP_FLAG: " << (p.flag_loop ? "true" : "false") << '\n'
       << "  STEPS_TO_OUTPUT: " << p.steps_to_output << '\n'
       << "  FLAG_SAVE_CONFIG: " << (p.flag_save_config ? "true" : "false");
    return os;
}

/** @brief Creates a list of magnetic field values for the sweep. */
std::vector<float> createSweepList(float start, float end, float step, bool loop);

/** @brief Performs a Monte Carlo step using chemical species exchange dynamics. 
 * @param lattice The lattice object representing the system.
 * @param params The simulation parameters.
 * @param tableBeg The pre-computed BEG site energy table.
 * @param DeltaEAcumM Accumulated energy change for magnetization.
 * @param changesAccepted Counter for accepted changes.
*/
void MonteCarloStepChemicalExchange(Lattice& lattice,
                                    const SimulationParameters& params,
                                    const FastBoltzmannTableBEG& table, 
                                    float& DeltaEAcumM,
                                    int& changesAccepted,
                                    int& changesAttempted);

/** @brief Performs a Monte Carlo step using spin flip dynamics.
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
                            const FastBoltzmannTableSpin& table,
                            float& DeltaEAcumM,
                            int& changesAccepted,
                            int& changesAttempted);

/** @brief Main simulation loop handling temperature and magnetic field sweeps.
 * @param params The simulation parameters.
 * @param nombrefile The base name for output files.
*/
void SimulationLoop(const SimulationParameters& params, 
                    const char* file_in, 
                    const char* file_out);

#endif // SIMULATION_H
