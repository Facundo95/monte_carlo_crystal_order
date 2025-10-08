#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "simulation.h"

using namespace std;

/**
 * @brief Reads simulation parameters and a list of input filenames from a text file.
 * * @param input_filename The path to the input configuration file (e.g., "input.txt").
 * @param params_out Reference to store the constructed SimulationParameters struct.
 * @param files_out Reference to store the list of base filenames for simulation runs.
 * @return true if the file was successfully read and all critical parameters were found.
 * @return false otherwise.
 */
bool readInputFile(const string& input_filename, SimulationParameters& params_out, vector<string>& files_out) {
    ifstream input_file(input_filename);
    if (!input_file.is_open()) {
        cerr << "ERROR: No se pudo abrir el archivo de entrada: " << input_filename << endl;
        return false;
    }

    // Temporary storage for required parameters before constructing SimulationParameters
    int num_steps = 0;
    float Jm3 = 0.0, Jm6 = 0.0;
    float T_upper = 0.0, T_lower = 0.0, step_T = 0.0;
    float H_upper = 0.0, H_lower = 0.0, step_H = 0.0;
    bool flag_save_config = false;
    int found_count = 0; // Tracks critical parameters found

    string line;
    while (getline(input_file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        string key;
        ss >> key;

        if (key == "NUM_STEPS" && ss >> num_steps) found_count++;
        else if (key == "J_M3" && ss >> Jm3) found_count++;
        else if (key == "J_M6" && ss >> Jm6) found_count++;
        else if (key == "T_UPPER" && ss >> T_upper) found_count++;
        else if (key == "T_LOWER" && ss >> T_lower) found_count++;
        else if (key == "STEP_T" && ss >> step_T) found_count++;
        else if (key == "H_UPPER" && ss >> H_upper) found_count++;
        else if (key == "H_LOWER" && ss >> H_lower) found_count++;
        else if (key == "STEP_H" && ss >> step_H) found_count++;
        else if (key == "FLAG_SAVE_CONFIG") {
            string val;
            if (ss >> val) {
                flag_save_config = (val == "true" || val == "TRUE" || val == "1");
                found_count++;
            }
        }
        else if (key == "FILE_ENTRY") {
            string filename;
            if (ss >> filename) {
                files_out.push_back(filename);
            }
        }
    }

    input_file.close();

    // Check if all 10 critical parameters were found
    if (found_count < 10) {
        cerr << "ERROR: Parámetros críticos incompletos en el archivo de entrada. Encontrados: " << found_count << "/10." << endl;
        return false;
    }
    if (files_out.empty()) {
        cerr << "ERROR: No se encontraron archivos de simulación (FILE_ENTRY) en el archivo de entrada." << endl;
        return false;
    }

    // If successful, initialize the SimulationParameters struct
    params_out = SimulationParameters(
        num_steps, Jm3, Jm6, 
        T_upper, T_lower, step_T, 
        H_upper, H_lower, step_H, 
        flag_save_config
    );

    return true;
}