#ifndef FILE_HANDLER_H
#define FILE_HANDLER_H

#include <fstream>
#include <string>

// Declarations for file opening helper functions
bool OpenLROParametersFile(const char* nombrefile, std::ofstream& output_stream);
bool OpenFinalRedFile(const char* nombrefile, float Hache, float TEMPERA, int count, std::ofstream& output_stream);

#endif // FILE_HANDLER_H
