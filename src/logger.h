#pragma once

#include <string>

namespace slog {
    // Start logging to `filename`. If a file with that name exists a new
    // unique filename is chosen (filename or filename(n)) and opened. The
    // existing file is left untouched. Returns true when the log file was
    // opened successfully.
    bool start(const std::string &filename);

    // Stop logging, restore original streams and close file.
    void stop();

    // Return an available filename based on `path`. If `path` does not exist
    // the same string is returned. Otherwise returns `path(n)` where n is
    // the smallest positive integer making the filename available.
    std::string makeUnique(const std::string &path);
}
