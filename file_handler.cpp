#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

/**
 * @brief Constructs the output filename and attempts to open the ofstream.
 * * The file is opened in append mode (ios::app) so that subsequent runs 
 * with the same name continue writing to the existing file.
 * * @param nombrefile The base name of the simulation file (e.g., "cu-al-mn_...").
 * @param output_stream The std::ofstream object to be initialized and opened.
 * @return bool True if the file was successfully opened, false otherwise.
 */
bool OpenLROParametersFile(const char* nombrefile, std::ofstream& output_stream) {
    // 1. Construct the full filename: e.g., "LRO_cu-al-mn_....txt"
    std::string filename(nombrefile);
    auto pos = filename.find_last_of('.');
    std::string fileOUT = filename.substr(0, pos) + ".out";
    std::ifstream in(fileOUT);
    if (in.good()) {
        in.close();
        fileOUT = fileOUT + "_new.out";
        std::cerr << "WARNING: El archivo de salida para parámetros LRO ya existe, cambiando el nombre a: " 
                  << fileOUT << std::endl;
    } else {
        in.close();
    }

    // 2. Attempt to open the file in output and append mode
    output_stream.open(fileOUT, std::ios::out | std::ios::app);

    // 3. Check for errors
    if (!output_stream.is_open()) {
        std::cerr << "ERROR: No se pudo abrir el archivo de salida para parámetros LRO: " 
                  << fileOUT << std::endl;
        return false;
    }
    
    // 4. If opened successfully, write header

    output_stream << "# Step\tH\tT\t"
                  << "X_A\tX_BUp\tX_BDown\tX_C\t"
                  << "Y_A\tY_BUp\tY_BDown\tY_C\t"
                  << "Z_A\tZ_BUp\tZ_BDown\tZ_C\t"
                  << "Magnetization\tEtotal\n";

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
bool OpenFinalRedFile(const char* nombrefile, float Hache, float TEMPERA, 
                    int count, std::ofstream& output_stream) {
    // Construct the full filename: e.g., "DUMP_cu-al-mn_..._200.0_0.txt"
    // Format Hache with a single significant decimal (one digit after the decimal point)
    std::ostringstream hs;
    hs << std::fixed << std::setprecision(1) << Hache;
    std::string Hstr = hs.str();

    std::ostringstream ts;
    ts << std::fixed << std::setprecision(1) << TEMPERA;
    std::string TEMPERAstr = ts.str();

    std::string filefinal_h = "dump_" + Hstr + "H_" + TEMPERAstr + "K_" + std::to_string(count) + std::string(nombrefile);
    
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
