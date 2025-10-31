#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "simulation.h"

bool readInputFile(const std::string& input_filename, 
                    SimulationParameters& params_out, 
                    std::vector<std::string>& files_out);

#endif // INPUT_PARSER_H