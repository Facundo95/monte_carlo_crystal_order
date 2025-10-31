#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <array>
#include <string>
#include <fstream>

/**
 * @brief Encapsulates the lattice data and complex coordinate/neighbor logic.
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

    // Initialization
    void loadInitialConfiguration(const std::string& filename);
    void initializeNeighbors();

    // Core simulation utilities
    int totalSites() const { return m_total_sites; }
    const std::array<int, 8>& getNeighbors1(int site) const { return neighbors1[site]; }
    int getSpin(int site) const { return magn_flat[site]; }
    int getSpecies(int site) const { return red_flat[site]; }   
    void flipSpin(int site) { magn_flat[site] = -magn_flat[site]; }

    void exchangeSpecies(int site1, int site2) {
        std::swap(red_flat[site1], red_flat[site2]);
        std::swap(magn_flat[site1], magn_flat[site2]);
    }

    float calculateSiteChemicalEnergy(int type, 
                                      float jota1, float jota2, 
                                      float ka1, float ka2, 
                                      float ele1, float ele2,
                                      float sumLin1, float sumCuad1,
                                      float sumLin2, float sumCuad2) const;

    float calculateSiteMagneticEnergy(int Spin, 
                                      float j3, float j6, 
                                      float H, 
                                      float sum3, float sum6) const;

    float calculateNeighborSpinSum(int site, int shell_type) const; 
    float calculateNeighborSpeciesSum(int site, int shell_type, int order) const;
    
    void saveFinalConfiguration(const char* nombrefile, float Hache, float TEMPERA, int count);
    
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

    inline int idx3D(int x, int y, int z) const {
        return (x * m_side + y) * m_depth + z;
    }

    static inline void idxToXYZ(int index, int side, int &x, int &y, int &z) {
        int depth = 2 * side;
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
};

#endif // LATTICE_H