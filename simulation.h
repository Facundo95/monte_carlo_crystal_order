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

    BoltzmannDeltaETable(const std::vector<double>& all_dEs, double T_, double eps_ = 1e-6)
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
    double w1_12; // Interaction energy between species 1 and 2 at 1st NN
    double w2_12; // Interaction energy between species 1 and 2 at 2nd NN
    double w1_13; // Interaction energy between species 1 and 3 at 1st NN
    double w2_13; // Interaction energy between species 1 and 3 at 2nd NN
    double w1_23; // Interaction energy between species 2 and 3 at 1st NN
    double w2_23; // Interaction energy between species 2 and 3 at 2nd NN
    double Jm1; // 1st NN magnetic interaction energy
    double Jm2; // 2nd NN magnetic interaction energy
    double Jm3; // 3rd NN interaction energy
    double Jm4; // 4th NN magnetic interaction energy
    double Jm5; // 5th NN magnetic interaction energy
    double Jm6; // 6th NN interaction energy
    double T_start;
    double T_end;
    double step_T;
    double H_start;
    double H_end;
    double step_H;
    int steps_to_output;
    bool flag_save_config;
    bool flag_loop;
    // Element symbols defined in the input file (e.g., "Cu", "Ni", "Al").
    std::string atom1;
    std::string atom2;
    std::string atom3;
    
    // Pre-computed BEG energy coefficients (computed once at initialization)
    double jota1; ///< Coefficient for w1_13 interaction
    double jota2; ///< Coefficient for w2_13 interaction
    double ka1;   ///< Interaction coefficient combining w1_12, w1_23, w1_13
    double ka2;   ///< Interaction coefficient combining w2_12, w2_23, w2_13
    double ele1;  ///< Electrostatic-like coefficient from w1_12 - w1_23
    double ele2;  ///< Electrostatic-like coefficient from w2_12 - w2_23

    SimulationParameters(int steps, int sim, int side, double w1_12, double w2_12, double w1_13, 
                         double w2_13, double w1_23, double w2_23, double j1, double j2, double j3, double j4, double j5, double j6,
                         double t_s, double t_e, double dt, double h_start, double h_end, 
                         double dh, int step_out, bool flag_red, bool loop)
        : num_steps(steps), simulation_method(sim), lattice_side(side), 
        w1_12(w1_12), w2_12(w2_12), w1_13(w1_13), w2_13(w2_13), w1_23(w1_23), w2_23(w2_23), 
        Jm1(j1), Jm2(j2), Jm3(j3), Jm4(j4), Jm5(j5), Jm6(j6), T_start(t_s), T_end(t_e), step_T(dt),H_start(h_start), H_end(h_end), 
        step_H(dh), steps_to_output(step_out), flag_save_config(flag_red), flag_loop(loop),
        jota1(0.25 * w1_13),
        jota2(0.25 * w2_13),
        ka1(0.25 * (2 * w1_12 + 2 * w1_23 - w1_13)),
        ka2(0.25 * (2 * w2_12 + 2 * w2_23 - w2_13)),
        ele1(0.25 * (w1_12 - w1_23)),
        ele2(0.25 * (w2_12 - w2_23)) {}
};

inline std::ostream& operator<<(std::ostream& os, const SimulationParameters& p) {
    os << "  NUM_STEPS: " << p.num_steps << '\n'
       << "  SIMULATION_METHOD: " << p.simulation_method << ((p.simulation_method == 0) ? " (Chemical Exchange)" : ((p.simulation_method == 1) ? " (Spin Flip)" : " Chemical + Spin Flip")) << '\n'
       << "  LATTICE_SIDE: " << p.lattice_side << '\n'
       << "  W1_12: " << p.w1_12 << " kB" << '\n'
       << "  W2_12: " << p.w2_12 << " kB" << '\n'
       << "  W1_13: " << p.w1_13 << " kB" << '\n'
       << "  W2_13: " << p.w2_13 << " kB" << '\n'
       << "  W1_23: " << p.w1_23 << " kB" << '\n'
       << "  W2_23: " << p.w2_23 << " kB" << '\n'
       << "  J_M1: " << p.Jm1 << " kB" << '\n'
       << "  J_M2: " << p.Jm2 << " kB" << '\n'
       << "  J_M3: " << p.Jm3 << " kB" << '\n'
       << "  J_M4: " << p.Jm4 << " kB" << '\n'
       << "  J_M5: " << p.Jm5 << " kB" << '\n'
       << "  J_M6: " << p.Jm6 << " kB" << '\n'
       << "  T_START: " << p.T_start << " K" << '\n'
       << "  T_END: " << p.T_end << " K" << '\n'
       << "  STEP_T: " << p.step_T << " K" << '\n'
       << "  H_UPPER: " << p.H_start << " muB" << '\n'
       << "  H_LOWER: " << p.H_end << " muB" << '\n'
       << "  STEP_H: " << p.step_H << " muB" << '\n'
       << "  LOOP_FLAG: " << (p.flag_loop ? "true" : "false") << '\n'
       << "  STEPS_TO_OUTPUT: " << p.steps_to_output << '\n'
       << "  FLAG_SAVE_CONFIG: " << (p.flag_save_config ? "true" : "false");
    return os;
}

/** 
 * @brief Generic template to create a sweep list for any parameter.
 * @tparam T The numeric type (double, int, etc.)
 * @param start Starting value
 * @param end Ending value
 * @param step Step size (must be positive)
 * @param loop If true, sweep goes from start to end and back to start. If false, just start to end.
 * @return Vector containing sweep values
 * @throws std::invalid_argument if step <= 0
 */
template<typename T>
std::vector<T> createSweepList(T start, T end, T step, bool loop = false);

/** @brief Performs a Monte Carlo step using chemical species exchange dynamics. 
 * @param lattice The lattice object representing the system.
 * @param params The simulation parameters.
 * @param tableBeg The pre-computed BEG site energy table.
 * @param DeltaEAcumM Accumulated energy change for magnetization.
 * @param changesAccepted Counter for accepted changes.
*/
// (See MCStepResults below) -- prototype using stats struct provided further down.

/**
 * @struct MCStepResults
 * @brief Aggregates counters and accumulators modified during a Monte Carlo step.
 * Stored as references so callers may keep their existing variables and pass
 * them bundled into a single parameter.
 */
struct MCStepResults {
    double& DeltaEAcum;
    int& changesAccepted;
    int& changesAttempted;

    MCStepResults(double& dE, int& acc, int& att)
        : DeltaEAcum(dE), changesAccepted(acc), changesAttempted(att) {}
};

// Updated prototype: pass results as a single struct reference
void MonteCarloStepChemicalExchange(Lattice& lattice,
                                    const SimulationParameters& params,
                                    BoltzmannDeltaETable& table,
                                    MCStepResults& stats);

/** @brief Performs a Monte Carlo step using spin flip dynamics.
 * @param lattice The lattice object representing the system.
 * @param H The external magnetic field.
 * @param params The simulation parameters.
 * @param table The pre-computed spin Boltzmann table.
 * @param DeltaEAcumM Accumulated energy change for magnetization.
 * @param changesAccepted Counter for accepted changes.
*/
// Updated prototype: pass results as a single struct reference
void MonteCarloStepSpinExtH(Lattice& lattice,
                            double H,
                            const SimulationParameters& params,
                            BoltzmannDeltaETable& table,
                            MCStepResults& stats);

/**
 * @struct NeighborSpeciesSums
 * @brief Holds precomputed neighbor species sums for a site and its chosen neighbor.
 */
struct NeighborSpeciesSums {
    int linNN_A;
    int linNN_N;
    int cuadNN_A;
    int cuadNN_N;
    int linNNN_A;
    int linNNN_N;
    int cuadNNN_A;
    int cuadNNN_N;
};

/**
 * @brief Compute and return neighbor species sums needed for chemical energy.
 * @param lattice The lattice
 * @param site The active site index
 * @param siteNeighbor The neighbor site index being exchanged with
 */
NeighborSpeciesSums computeNeighborSpeciesSums(const Lattice& lattice, int site, int siteNeighbor);

/**
 * @brief Compute and return spin neighbor sums for shells 1..6 for a given site.
 */
std::array<int,6> computeNeighborSpinSums(const Lattice& lattice, int site);

/** @brief Main simulation loop handling temperature and magnetic field sweeps.
 * @param params The simulation parameters.
 * @param nombrefile The base name for output files.
*/
void SimulationLoop(const SimulationParameters& params, 
                    const char* file_in, 
                    const char* file_out);

#endif // SIMULATION_H
