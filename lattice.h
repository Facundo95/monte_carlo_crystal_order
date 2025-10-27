#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <array>
#include <string>

using namespace std;

// Define system dimensions
static constexpr int LATTICE_SIDE = 32;
static constexpr int LATTICE_DEPTH = 64; // 2 * LATTICE_SIDE
static constexpr int LATTICE_TOTAL_SITES = LATTICE_SIDE * LATTICE_SIDE * LATTICE_DEPTH;

/**
 * @brief Encapsulates the lattice data and complex coordinate/neighbor logic.
 */
class Lattice {
private:
    int magn_flat[LATTICE_TOTAL_SITES]; // Flattened magnetic state for memory efficiency
    int red_flat[LATTICE_TOTAL_SITES];  // Flattened species state for memory efficiency

    vector<array<int, 8>>  neighbors1;
    vector<array<int, 6>>  neighbors2;
    vector<array<int, 12>> neighbors3;
    vector<array<int, 6>>  neighbors6;

    static inline int idx3D(int x, int y, int z, int side = LATTICE_SIDE, int depth = LATTICE_DEPTH) {
        // index plano: (x * side + y) * depth + z
        return (x * side + y) * depth + z;
    }

    static inline void idxToXYZ(int index, int side, int depth, int &x, int &y, int &z) {
        const int planeSize = side * depth;

        x = index / planeSize;
        int tmp = index - x * planeSize;
        y = tmp / depth;
        z = tmp - y * depth;
    }
    
    static inline int wrap(int v, int M) {
        if (v < 0) return v + M;
        if (v >= M) return v - M;
        return v;
    }

public:
    Lattice() {
        const int N = LATTICE_SIDE * LATTICE_SIDE * LATTICE_DEPTH;
        neighbors1.resize(N);
        neighbors2.resize(N);
        neighbors3.resize(N);
        neighbors6.resize(N);
    }
    // Initialization
    void loadInitialConfiguration(const string& filename);
    void initializeNeighbors(); // Pre-compute neighbor indices for efficiency

    // Core simulation utilities
    int getSpin(int site) const { return magn_flat[site]; }
    void flipSpin(int site) { magn_flat[site] = -magn_flat[site]; }
    
    // Calculates the required spin sum (used in delta E calculation)
    float calculateNeighborSpinSum(int site, int shell_type) const; 
    
    // Output utility
    void saveFinalConfiguration(const char* nombrefile, float Hache, int count);
    
    // Measurement utility (implemented in simulation.cpp)
    void calculateAndWriteLRO(ofstream& parout, 
                            int step_count, 
                            float T, 
                            float H, 
                            float DeltaEAcumM) const;
};

#endif // LATTICE_H