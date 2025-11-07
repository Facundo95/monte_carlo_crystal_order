#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "simulation.h"

/**
 * @brief Reads simulation parameters and a list of input filenames from a text file.
 * * @param input_filename The path to the input configuration file (e.g., "input.txt").
 * @param params_out Reference to store the constructed SimulationParameters struct.
 * @param files_out Reference to store the list of base filenames for simulation runs.
 * @return true if the file was successfully read and all critical parameters were found.
 * @return false otherwise.
 */

bool readInputFile(const std::string& input_filename, 
                    SimulationParameters& params_out, 
                    std::vector<std::string>& files_out);

#endif // INPUT_PARSER_H