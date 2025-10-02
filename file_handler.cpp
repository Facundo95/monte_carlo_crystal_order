#include <iostream>
#include <fstream>
#include <string>

/**
 * @brief Constructs the output filename and attempts to open the ofstream.
 * * The file is opened in append mode (ios::app) so that subsequent runs 
 * with the same name continue writing to the existing file.
 * * @param nombrefile The base name of the simulation file (e.g., "cu-al-mn_...").
 * @param output_stream The std::ofstream object to be initialized and opened.
 * @return bool True if the file was successfully opened, false otherwise.
 */
bool OpenLROParametersFile(const char* nombrefile, std::ofstream& output_stream) {
    // 1. Construct the full filename: e.g., "cu-al-mn_..._out.txt"
    std::string fileOUT = std::string(nombrefile) + "_out.txt";
    
    // 2. Attempt to open the file in output and append mode
    output_stream.open(fileOUT, std::ios::out | std::ios::app);
    
    // 3. Check for errors
    if (!output_stream.is_open()) {
        std::cerr << "ERROR: No se pudo abrir el archivo de salida para parámetros LRO: " 
                  << fileOUT << std::endl;
        return false;
    }
    
    return true;
}

/**
 * @brief Constructs the final configuration filename and attempts to open the ofstream.
 * * @param nombrefile The base name of the simulation file.
 * @param Hache The magnetic field value to include in the filename.
 * @param count A counter or index to include in the filename.
 * @param output_stream The std::ofstream object to be initialized and opened.
 * @return bool True if the file was successfully opened, false otherwise.
 */
bool OpenFinalRedFile(const char* nombrefile, float Hache, int count, std::ofstream& output_stream) {
    // Construct the full filename: e.g., "cu-al-mn_..._200.0_0.txt"
    std::string filefinal_h = std::string(nombrefile) + "_" + std::to_string(Hache) + "_" + std::to_string(count) + ".txt";
    
    // Attempt to open the file in output (default overwrite mode)
    output_stream.open(filefinal_h, std::ios::out);
    
    // Check for errors
    if (!output_stream.is_open()) {
        std::cerr << "ERROR: No se pudo abrir el archivo de salida de la red final: "
                  << filefinal_h << std::endl;
        return false;
    }
    return true;
}
