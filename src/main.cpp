#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream> // For file input
#include <sstream> // For string parsing
#include <chrono>
#include <ctime>
#include <algorithm>

#include "simulation.h"
#include "input_parser.h"
#include "time_calculator.h"
#include "logger.h"


int main(int argc, char* argv[]) {
    std::string inputFile;
    std::string logFilename = "log.output";

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-in" && i + 1 < argc) {
            inputFile = argv[++i];
        } else if (a == "-log" && i + 1 < argc) {
            logFilename = argv[++i];
        }
    }

    if (inputFile.empty()) {
        std::cerr << "# Modo de uso: ./mc_simulation -in <input_file_name> [-log <log_file_name>]" << std::endl;
        return 1;
    }

    // start logging (rotates existing file)
    if (!slog::start(logFilename)) {
        std::cerr << "Warning: could not open " << logFilename << " for writing. Continuing without log." << std::endl;
    }

    std::cout << R"(   _________    ____  __    ____  _____    __  _______  _   ___________________
  / ____/   |  / __ \/ /   / __ \/ ___/   /  |/  / __ \/ | / /_  __/ ____/ ___/
 / /   / /| | / /_/ / /   / / / /\__ \   / /|_/ / / / /  |/ / / / / __/  \__ \
/ /___/ ___ |/ _, _/ /___/ /_/ /___/ /  / /  / / /_/ / /|  / / / / /___ ___/ / 
\____/_/  |_/_/ |_/_____\/____//____/  /_/  /_/\____/_/ |_/ /_/ /_____//____/  
                                                                               )" << std::endl;
    std::cout << "--- Carlos Montes Iniciando (C++17) ---" << std::endl;

    TimeCalculator timer;
    timer.start();

    // 1. Declare variables to hold read data
    SimulationParameters params(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, true, false, false);
    std::string file_in;
    std::string file_out;

    // 2. Read parameters and file list from "input.txt"
    if (!readInputFile(inputFile, params, file_in, file_out)) {
        std::cerr << "Falló la lectura del archivo de configuración. Terminando programa." << std::endl;
        return 1;
    }
    
    std::cout << "Archivo de entrada encontrado: " << inputFile << std::endl;
    std::cout << "Parametros de la simulacion:" << std::endl;
    std::cout << params << std::endl;

    // 3. Start the simulation loop using the read data
    std::cout << "Estructura inicial tomada de: " << file_in << std::endl;
    std::cout << "Archivo de salida de la simulacion: " << file_out << std::endl;
        
    // Pass the parameters to the main simulation loop
    SimulationLoop(params, file_in.c_str(), file_out.c_str());
    

    timer.stop();
    
    // Display the time taken in the format "hh:mm:ss"
    timer.displayDuration();

    // stop logging and restore streams
    slog::stop();

    return 0;
}
