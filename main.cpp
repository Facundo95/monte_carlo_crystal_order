#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <iomanip>
#include <fstream> // For file input
#include <sstream> // For string parsing

#include "simulation.h"
#include "input_parser.h"

using namespace std; 

int main(int argc, char* argv[]){
    if (argc != 3){
        cerr << "# Modo de uso: ./mc_simulation -in <input_file_name>" << endl;
    } else{
        cout << "--- Carlos Montes Iniciando (C++17) ---" << endl;
    }

    clock_t start = clock(); // Iniciamos el reloj

    // 1. Declare variables to hold read data
    SimulationParameters params(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false); // Dummy initialization
    vector<string> files;
    string inputFile = argv[2];

    // 2. Read parameters and file list from "input.txt"
    if (!readInputFile(inputFile, params, files)) {
        cerr << "Falló la lectura del archivo de configuración. Terminando programa." << endl;
        return 1;
    } else{
        cout << "Archivo de entrada encontrado: " << inputFile << endl;
        cout << "Parametros de la simulacion:" << endl;
        cout << params << endl;
    }

    // 3. Start the simulation loop using the read data
    for (const auto& file : files){
        cout << "Estructura inicial tomada de: " << file << endl;
        
        // Pass the parameters to the main simulation loop
        SimulationLoop(params, file.c_str());
    }

    clock_t end = clock(); // Cortamos el reloj
    double timeTaken = double(end - start) / CLOCKS_PER_SEC;
    
    // Display the time taken in the format "hh:mm:ss"
    int hours = int(timeTaken / 3600);
    int minutes = int((timeTaken - hours * 3600) / 60);
    int seconds = int(timeTaken - hours * 3600 - minutes * 60);
    cout << "-------------------------------" << endl;
    cout << "Tiempo total: "
              << setfill('0') << setw(2) << hours << ":"
              << setfill('0') << setw(2) << minutes << ":"
              << setfill('0') << setw(2) << seconds << endl;
    cout << "-------------------------------" << endl;
    return 0;
}
