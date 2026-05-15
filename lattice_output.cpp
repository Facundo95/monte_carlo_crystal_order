#include "lattice_output.h"

#include <iostream>

namespace lattice_output {

namespace {
// no-op helpers removed; extras support was removed.
}

void writeLROParameters(std::ofstream& parout,
                        int step_count,
                        double T,
                        double H,
                        const Lattice::LROParameters& lro,
                        double normalizedMagnetization,
                        double energyValue,
                        bool printToConsole) {
    parout << step_count << "\t" << H << "\t" << T << "\t"
           << lro.X_A << "\t" << lro.X_Bup << "\t" << lro.X_Bdown << "\t" << lro.X_C << "\t"
           << lro.Y_A << "\t" << lro.Y_Bup << "\t" << lro.Y_Bdown << "\t" << lro.Y_C << "\t"
           << lro.Z_A << "\t" << lro.Z_Bup << "\t" << lro.Z_Bdown << "\t" << lro.Z_C << "\t"
           << normalizedMagnetization << "\t"
           << energyValue;
    parout << "\t" << std::endl;

    if (printToConsole) {
        std::cout << "Step=" << step_count << " H=" << H << " T=" << T << " | "
                  << "X_A=" << lro.X_A << " X_Bup=" << lro.X_Bup << " X_Bdown=" << lro.X_Bdown << " X_C=" << lro.X_C
                  << " | mag=" << normalizedMagnetization << " | E=" << energyValue;
        std::cout << std::endl;
    }
}

void writeReducedOutput(std::ofstream& parout,
                        int step_count,
                        double T,
                        double H,
                        double normalizedMagnetization,
                        double energyValue,
                        bool printToConsole) {
    parout << step_count << "\t" << H << "\t" << T << "\t"
           << normalizedMagnetization << "\t" << energyValue;
    parout << "\t" << std::endl;

    if (printToConsole) {
        std::cout << "Step=" << step_count << " H=" << H << " T=" << T << " | mag=" << normalizedMagnetization << " | E=" << energyValue;
        std::cout << std::endl;
    }
}

} // namespace lattice_output