#ifndef FILE_HANDLER_H
#define FILE_HANDLER_H

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
bool OpenLROParametersFile(const char* nombrefile, std::ofstream& output_stream);

/**
 * @brief Constructs the final configuration filename and attempts to open the ofstream.
 * * @param nombrefile The base name of the simulation file.
 * @param Hache The magnetic field value to include in the filename.
 * @param count A counter or index to include in the filename.
 * @param output_stream The std::ofstream object to be initialized and opened.
 * @return bool True if the file was successfully opened, false otherwise.
 */
bool OpenFinalRedFile(const char* nombrefile, float Hache, float TEMPERA, int count, std::ofstream& output_stream);

#endif // FILE_HANDLER_H
