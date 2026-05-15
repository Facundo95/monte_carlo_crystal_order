#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "lattice.h"

/**
 * @brief Constructs the output filename and attempts to open the ofstream.
 * * The file is opened in append mode (ios::app) so that subsequent runs 
 * with the same name continue writing to the existing file.
 * * @param nombrefile The base name of the simulation file (e.g., "cu-al-mn_...").
 * @param output_stream The std::ofstream object to be initialized and opened.
 * @param lattice The Lattice object to access atom names.
 * @return bool True if the file was successfully opened, false otherwise.
 */
bool OpenLROParametersFile(const char* nombrefile, std::ofstream& output_stream, const Lattice& lattice) {
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
    
    // 4. If opened successfully, write header using atom names from lattice
    std::string X_A = "x_" + lattice.getAtom1();
    std::string X_BUp = "x_" + lattice.getAtom2() + "_up";
    std::string X_BDown = "x_" + lattice.getAtom2() + "_down";
    std::string X_C = "x_" + lattice.getAtom3();
    std::string Y_A = "y_" + lattice.getAtom1();
    std::string Y_BUp = "y_" + lattice.getAtom2() + "_up";
    std::string Y_BDown = "y_" + lattice.getAtom2() + "_down";
    std::string Y_C = "y_" + lattice.getAtom3();
    std::string Z_A = "z_" + lattice.getAtom1();
    std::string Z_BUp = "z_" + lattice.getAtom2() + "_up";
    std::string Z_BDown = "z_" + lattice.getAtom2() + "_down";
    std::string Z_C = "z_" + lattice.getAtom3();

    output_stream << "# step\th\ttemperature\t"
                  << X_A << "\t" << X_BUp << "\t" << X_BDown << "\t" << X_C << "\t"
                  << Y_A << "\t" << Y_BUp << "\t" << Y_BDown << "\t" << Y_C << "\t"
                  << Z_A << "\t" << Z_BUp << "\t" << Z_BDown << "\t" << Z_C << "\t"
                  << "magnetization\tetotal\n";

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
bool OpenFinalRedFile(const char* nombrefile, double Hache, double TEMPERA, 
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
