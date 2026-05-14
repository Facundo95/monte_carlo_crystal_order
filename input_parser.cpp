#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include "simulation.h"

/**
 * @brief Helper function to extract the next value from a stringstream,
 * skipping "=" or ":" delimiters if present.
 * 
 * @param ss The stringstream to read from.
 * @param value Reference to store the extracted value.
 * @return true if a value was successfully extracted, false otherwise.
 */
static bool extractValue(std::stringstream& ss, std::string& value) {
    std::string token;
    // Skip any "=" or ":" delimiters
    while (ss >> token) {
        if (token == "=" || token == ":") {
            continue; // Skip delimiter, try next token
        }
        value = token;
        return true;
    }
    return false;
}

/**
 * @brief Helper function to extract an integer value, skipping delimiters.
 * 
 * @param ss The stringstream to read from.
 * @param value Reference to store the extracted integer.
 * @return true if an integer was successfully extracted, false otherwise.
 */
static bool extractIntValue(std::stringstream& ss, int& value) {
    std::string token;
    while (ss >> token) {
        if (token == "=" || token == ":") {
            continue; // Skip delimiter, try next token
        }
        try {
            value = std::stoi(token);
            return true;
        } catch (const std::invalid_argument&) {
            return false;
        }
    }
    return false;
}

/**
 * @brief Helper function to extract a float value, skipping delimiters.
 * 
 * @param ss The stringstream to read from.
 * @param value Reference to store the extracted float.
 * @return true if a float was successfully extracted, false otherwise.
 */
static bool extractFloatValue(std::stringstream& ss, float& value) {
    std::string token;
    while (ss >> token) {
        if (token == "=" || token == ":") {
            continue; // Skip delimiter, try next token
        }
        try {
            value = std::stof(token);
            return true;
        } catch (const std::invalid_argument&) {
            return false;
        }
    }
    return false;
}

/**
 * @brief Helper function to extract a boolean value, skipping delimiters.
 * 
 * @param ss The stringstream to read from.
 * @param value Reference to store the extracted boolean.
 * @return true if a boolean was successfully extracted, false otherwise.
 */
static bool extractBoolValue(std::stringstream& ss, bool& value) {
    std::string token;
    while (ss >> token) {
        if (token == "=" || token == ":") {
            continue; // Skip delimiter, try next token
        }
        value = (token == "true" || token == "TRUE" || token == "1");
        return true;
    }
    return false;
}

/**
 * @brief Validates if a string is a valid atomic symbol (1-2 alphabetic characters).
 * Valid examples: "H", "He", "Cu", "Ni", "Al", etc.
 * Invalid examples: "Niquel", "nikel", "123", etc.
 * 
 * @param symbol The string to validate.
 * @return true if the symbol is valid, false otherwise.
 */
static bool isValidAtomicSymbol(const std::string& symbol) {
    // Must be 1 or 2 characters long
    if (symbol.length() < 1 || symbol.length() > 2) {
        return false;
    }
    
    // All characters must be alphabetic
    for (char c : symbol) {
        if (!std::isalpha(c)) {
            return false;
        }
    }
    
    return true;
}

/**
 * @brief Reads simulation parameters and a list of input filenames from a text file.
 * * @param input_filename The path to the input configuration file (e.g., "input.txt").
 * @param params_out Reference to store the constructed SimulationParameters struct.
 * @param file_in Reference to filename for simulation run.
 * @param file_out Reference to filename for output.
 * @return true if the file was successfully read and all critical parameters were found.
 * @return false otherwise.
 */
bool readInputFile(const std::string& input_filename, SimulationParameters& params_out, 
                    std::string& file_in, std::string& file_out) {

    std::ifstream input_file(input_filename);
    if (!input_file.is_open()) {
        std::cerr << "ERROR: No se pudo abrir el archivo de entrada: " << input_filename << std::endl;
        return false;
    }

    // Temporary storage for required parameters before constructing SimulationParameters
    int num_steps = 100;  // Default value
    int simulation_method = 1;  // Default: Ising Hamiltonian
    int lattice_side = 32;  // Default small lattice
    float w1_12 = 0.0, w2_12 = 0.0, w1_13 = 0.0, w2_13 = 0.0, w1_23 = 0.0, w2_23 = 0.0;
    float Jm3 = 0.0, Jm6 = 0.0;
    float T_start = 0.0, T_end = 0.0, step_T = 0.0;
    float H_start = 0.0, H_end = 0.0, step_H = 0.0;
    bool flag_save_config = false;
    bool flag_loop = false;
    int steps_to_output = 100;  // Default value
    int found_count = 0; // Tracks MANDATORY parameters found
    std::string atom_1 = "";
    std::string atom_2 = "";
    std::string atom_3 = "";
    bool has_t_end = false;  // Track if T_END was explicitly set
    bool has_h_end = false;  // Track if H_END was explicitly set

    std::string line;
    int line_number = 0;
    while (getline(input_file, line)) {
        line_number++;
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string key;
        ss >> key;

        if (key == "NUM_STEPS") {
            extractIntValue(ss, num_steps);
        }
        else if (key == "SIMULATION_METHOD") {
            extractIntValue(ss, simulation_method);
        }
        else if (key == "LATTICE_SIDE") {
            extractIntValue(ss, lattice_side);
        }
        else if (key == "W1_12") {
            extractFloatValue(ss, w1_12);
        }
        else if (key == "W1_13") {
            extractFloatValue(ss, w1_13);
        }
        else if (key == "W1_23") {
            extractFloatValue(ss, w1_23);
        }
        else if (key == "W2_12") {
            extractFloatValue(ss, w2_12);
        }
        else if (key == "W2_13") {
            extractFloatValue(ss, w2_13);
        }
        else if (key == "W2_23") {
            extractFloatValue(ss, w2_23);
        }
        else if (key == "J_M3") {
            if (extractFloatValue(ss, Jm3)) {
                found_count++;  // MANDATORY
            } else {
                std::cerr << "WARNING: J_M3 has invalid value at line " << line_number << std::endl;
            }
        }
        else if (key == "J_M6") {
            extractFloatValue(ss, Jm6);
        }
        else if (key == "T_START") {
            if (extractFloatValue(ss, T_start)) {
                found_count++;  // MANDATORY
            } else {
                std::cerr << "WARNING: T_START has invalid value at line " << line_number << std::endl;
            }
        }
        else if (key == "T_END") {
            if (extractFloatValue(ss, T_end)) {
                has_t_end = true;
            } else {
                std::cerr << "WARNING: T_END has invalid value at line " << line_number << std::endl;
            }
        }
        else if (key == "STEP_T") {
            extractFloatValue(ss, step_T);
        }
        else if (key == "H_START") {
            if (extractFloatValue(ss, H_start)) {
                found_count++;  // MANDATORY
            } else {
                std::cerr << "WARNING: H_START has invalid value at line " << line_number << std::endl;
            }
        }
        else if (key == "H_END") {
            if (extractFloatValue(ss, H_end)) {
                has_h_end = true;
            } else {
                std::cerr << "WARNING: H_END has invalid value at line " << line_number << std::endl;
            }
        }
        else if (key == "STEP_H") {
            extractFloatValue(ss, step_H);
        }
        else if (key == "STEPS_TO_OUTPUT") {
            extractIntValue(ss, steps_to_output);
        }
        else if (key == "FLAG_SAVE_CONFIG") {
            extractBoolValue(ss, flag_save_config);
        }
        else if (key == "LOOP") {
            extractBoolValue(ss, flag_loop);
        }
        else if (key == "FILE_ENTRY") {
            std::string filename;
            if (!extractValue(ss, filename)) {
                std::cerr << "WARNING: FILE_ENTRY has no value at line " << line_number << std::endl;
            } else {
                file_in = filename;
                found_count++;  // MANDATORY
            }
        }
        else if (key == "FILE_OUTPUT") {
            std::string filename;
            if (!extractValue(ss, filename)) {
                std::cerr << "WARNING: FILE_OUTPUT has no value at line " << line_number << std::endl;
            } else {
                file_out = filename;
            }
        }
        else if (key == "ATOM_1") {
            if (!extractValue(ss, atom_1)) {
                std::cerr << "WARNING: ATOM_1 has no value at line " << line_number << std::endl;
            } else if (!isValidAtomicSymbol(atom_1)) {
                std::cerr << "WARNING: ATOM_1 '" << atom_1 << "' is not a valid atomic symbol (1-2 alphabetic characters) at line " << line_number << std::endl;
                atom_1 = ""; // Discard invalid value
            } else {
                found_count++;  // MANDATORY
            }
        }
        else if (key == "ATOM_2") {
            if (!extractValue(ss, atom_2)) {
                std::cerr << "WARNING: ATOM_2 has no value at line " << line_number << std::endl;
            } else if (!isValidAtomicSymbol(atom_2)) {
                std::cerr << "WARNING: ATOM_2 '" << atom_2 << "' is not a valid atomic symbol (1-2 alphabetic characters) at line " << line_number << std::endl;
                atom_2 = ""; // Discard invalid value
            } else {
                found_count++;  // MANDATORY
            }
        }
        else if (key == "ATOM_3") {
            if (!extractValue(ss, atom_3)) {
                std::cerr << "WARNING: ATOM_3 has no value at line " << line_number << std::endl;
            } else if (!isValidAtomicSymbol(atom_3)) {
                std::cerr << "WARNING: ATOM_3 '" << atom_3 << "' is not a valid atomic symbol (1-2 alphabetic characters) at line " << line_number << std::endl;
                atom_3 = ""; // Discard invalid value
            } else {
                found_count++;  // MANDATORY
            }
        }
    }

    input_file.close();

    if (file_in.empty()) {
        std::cerr << "ERROR: No se encontraron archivos de simulación (FILE_ENTRY) en el archivo de entrada." << std::endl;
        return false;
    }

    // Check for mandatory parameters (ATOM_1, ATOM_2, ATOM_3, J_M3, T_START, H_START, FILE_ENTRY = 7 total)
    if (found_count < 7) {
        std::cerr << "ERROR: Missing mandatory parameters. Expected 7, found " << found_count << std::endl;
        std::cerr << "Mandatory parameters: ATOM_1, ATOM_2, ATOM_3, J_M3, T_START, H_START, FILE_ENTRY" << std::endl;
        return false;
    }

    // If T_END was not explicitly set, use T_START
    if (!has_t_end) {
        T_end = T_start;
    }

    // If H_END was not explicitly set, use H_START
    if (!has_h_end) {
        H_end = H_start;
    }

    // If successful, initialize the SimulationParameters struct
    params_out = SimulationParameters(
        num_steps, simulation_method, lattice_side,
        w1_12, w2_12, w1_13, w2_13, w1_23, w2_23, 
        Jm3, Jm6, 
        T_start, T_end, step_T, 
        H_start, H_end, step_H, 
        steps_to_output, flag_save_config, flag_loop
    );

    // Attach element symbols (if provided) to the params object
    params_out.atom1 = atom_1;
    params_out.atom2 = atom_2;
    params_out.atom3 = atom_3;

    return true;
}