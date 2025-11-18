#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <iomanip>
#include <fstream> // For file input
#include <sstream> // For string parsing

#include "simulation.h"
#include "input_parser.h"

int main(int argc, char* argv[]){
    if (argc != 3){
        std::cerr << "# Modo de uso: ./mc_simulation -in <input_file_name>" << std::endl;
    } else{
        std::cout << "--- Carlos Montes Iniciando (C++17) ---" << std::endl;
    }

    clock_t start = clock(); // Iniciamos el reloj

    // 1. Declare variables to hold read data
    SimulationParameters params(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false);
    std::string file_in;
    std::string file_out;
    std::string inputFile = argv[2];

    // 2. Read parameters and file list from "input.txt"
    if (!readInputFile(inputFile, params, file_in, file_out)) {
        std::cerr << "Falló la lectura del archivo de configuración. Terminando programa." << std::endl;
        return 1;
    } else{
        std::cout << "Archivo de entrada encontrado: " << inputFile << std::endl;
        std::cout << "Parametros de la simulacion:" << std::endl;
        std::cout << params << std::endl;
    }

    // 3. Start the simulation loop using the read data
    std::cout << "Estructura inicial tomada de: " << file_in << std::endl;
    std::cout << "Archivo de salida de la simulacion: " << file_out << std::endl;
        
    // Pass the parameters to the main simulation loop
    SimulationLoop(params, file_in.c_str(), file_out.c_str());
    

    clock_t end = clock(); // Cortamos el reloj
    double timeTaken = double(end - start) / CLOCKS_PER_SEC;
    
    // Display the time taken in the format "hh:mm:ss"
    int hours = int(timeTaken / 3600);
    int minutes = int((timeTaken - hours * 3600) / 60);
    int seconds = int(timeTaken - hours * 3600 - minutes * 60);
    std::cout << "-------------------------------" << std::endl;
    std::cout << "Tiempo total: "
              << std::setfill('0') << std::setw(2) << hours << ":"
              << std::setfill('0') << std::setw(2) << minutes << ":"
              << std::setfill('0') << std::setw(2) << seconds << std::endl;
    std::cout << "-------------------------------" << std::endl;
    return 0;
}
