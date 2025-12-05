#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <array>
#include <string>
#include <fstream>

/***
 * @class Lattice
 * @brief Represents a 3D BCC lattice for Monte Carlo simulations.
 *
 * This class manages the lattice structure, including spin and species
 * configurations, neighbor relationships, and energy calculations.
 */
class Lattice {
public:
    explicit Lattice(int LATTICE_SIDE)
        : m_side(LATTICE_SIDE),
          m_depth(2 * LATTICE_SIDE),
          m_total_sites(2 * LATTICE_SIDE * LATTICE_SIDE * LATTICE_SIDE),
          magn_flat(m_total_sites),
          red_flat(m_total_sites),
          neighbors1(m_total_sites),
          neighbors2(m_total_sites),
          neighbors3(m_total_sites),
          neighbors6(m_total_sites)
    {}

    /** @brief Loads the initial configuration from a specified file.
     * The file should contain the species and spin states for each lattice site.
     * @param filename The path to the input configuration file.
     */
    void loadInitialConfiguration(const std::string& filename);
    
    /** @brief Initializes the neighbor lists for each lattice site. */
    void initializeNeighbors();

    /** @brief Returns the total number of sites in the lattice. */
    int totalSites() const { return m_total_sites; }

    /** @brief Accessor methods for neighbors. */
    const std::array<int, 8>& getNeighbors1(int site) const { return neighbors1[site]; }
    
    /** @brief Accessors for spin. */
    int getSpin(int site) const { return magn_flat[site]; }

    /** @brief Accessors for species. */
    int getSpecies(int site) const { return red_flat[site]; }   
    
    /** @brief Mutator for spin. */
    void flipSpin(int site) { magn_flat[site] = -magn_flat[site]; }

    /** @brief Mutator for species exchange. */
    void exchangeSpecies(int site1, int site2) {
        std::swap(red_flat[site1], red_flat[site2]);
        std::swap(magn_flat[site1], magn_flat[site2]);
    }

    /** @brief Calculates the chemical energy contribution for a given site. */
    double calculateDeltaChemicalEnergy(int type_A, int type_N, 
                                    float JOTA1, float JOTA2, 
                                    float KA1, float KA2, 
                                    float ELE1, float ELE2,
                                    int sumLinNN_A, int sumLinNNN_A,
                                    int sumCuadNN_A, int sumCuadNNN_A,
                                    int sumLinNN_N, int sumLinNNN_N,
                                    int sumCuadNN_N, int sumCuadNNN_N) const;
    
    /** @brief Calculates the magnetic energy contribution for a given site. */
    float calculateSiteMagneticEnergy(int Spin, 
                                      float j3, float j6, 
                                      float H, 
                                      float sum3, float sum6) const;
    
    /** @brief Calculates the sum of neighbor spins for a given shell type. */
    float calculateNeighborSpinSum(int site, int shell_type) const;
    
    /** @brief Calculates the sum of neighbor species for a given shell type and order. */
    float calculateNeighborSpeciesSum(int site, int shell_type, int order) const;
    
    /** @brief Saves the final configuration to a file and returns the opened ofstream.
     * The output is written in an extended .xyz-like format: first line is the
     * number of atoms, second line is a comment containing H/T/count, and each
     * following line contains: Element x y z species spin
     * The returned std::ofstream is moved to the caller and remains open.
     */
    bool saveFinalConfiguration(const char* nombrefile, 
                                float Hache, float TEMPERA, 
                                int count);
    
    /** @brief Calculates and writes LRO parameters to the output file. */
    void calculateAndWriteLRO(std::ofstream& parout, 
                              int step_count, 
                              float T, 
                              float H, 
                              float DeltaEAcumM) const;

private:
    int m_side;
    int m_depth;
    int m_total_sites;

    std::vector<int> magn_flat;
    std::vector<int> red_flat;
    std::vector<std::array<int, 8>>  neighbors1;
    std::vector<std::array<int, 6>>  neighbors2;
    std::vector<std::array<int, 12>> neighbors3;
    std::vector<std::array<int, 6>>  neighbors6;
    
    /** @brief Helper to convert 3D indices to flat index */
    inline int idx3D(int x, int y, int z) const {
        return (x * m_side + y) * m_depth + z;
    }

    /** @brief Helper to convert flat index to 3D indices */
    static inline void idxToXYZ(int index, int side, int &x, int &y, int &z) {
        int depth = 2 * side;
        const int planeSize = side * depth;
        x = index / planeSize;
        int tmp = index - x * planeSize;
        y = tmp / depth;
        z = tmp - y * depth;
    }

    /** @brief Helper to wrap coordinates with periodic boundary conditions */
    static inline int wrap(int v, int M) {
        if (v < 0) return v + M;
        if (v >= M) return v - M;
        return v;
    }
};

#endif // LATTICE_H