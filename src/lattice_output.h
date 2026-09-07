#ifndef LATTICE_OUTPUT_H
#define LATTICE_OUTPUT_H

#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "lattice.h"

namespace lattice_output {

// Output functions do not accept extra labeled parameters anymore.
void writeProgressHeader();

void writeProgressRow(const Lattice& lattice,
                      int sweep,
                      double T,
                      double H,
                      double acceptancePercentage,
                      double energyValue);

void writeLROParameters(std::ofstream& parout,
                        int step_count,
                        double T,
                        double H,
                        const Lattice::LROParameters& lro,
                        double normalizedMagnetization,
                        double energyValue,
                        bool printToConsole = false);

void writeReducedOutput(std::ofstream& parout,
                        int step_count,
                        double T,
                        double H,
                        double normalizedMagnetization,
                        double energyValue,
                        bool printToConsole = false);

}

#endif // LATTICE_OUTPUT_H