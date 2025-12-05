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
    // Use quantized integer keys to avoid exact-bit mismatches caused by
    // floating point rounding differences between the precompute and runtime.
    std::unordered_map<long long, double> table;
    double eps; // quantization resolution used to generate keys
    double T;   // temperature used to compute Boltzmann factors (for fallback)

    static inline long long key_of(double x, double eps) {
        return static_cast<long long>(std::llround(x / eps));
    }

    BoltzmannDeltaETable(const std::vector<double>& all_dEs, double T_, double eps_ = 1e-7)
        : eps(eps_), T(T_)
    {
        if (T == 0.0) T = 1e-12;
        for (double dE : all_dEs) {
            long long key = key_of(dE, eps);
            table[key] = std::exp(-dE / T);
        }
    }

    inline double lookup(double dE) {
        long long key = key_of(dE, eps);
        auto it = table.find(key);
        if (it != table.end()) return it->second;

        // Fallback: compute Boltzmann factor on the fly and insert it
        double Tloc = T;
        if (Tloc == 0.0) Tloc = 1e-12;
        double val = std::exp(-dE / Tloc);
        table[key] = val;
        return val;
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
                                    BoltzmannDeltaETable& table,
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
                            BoltzmannDeltaETable& table,
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
