#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "simulation.h"

using namespace std;

bool readInputFile(const string& input_filename, 
                    SimulationParameters& params_out, 
                    vector<string>& files_out);

#endif // INPUT_PARSER_H