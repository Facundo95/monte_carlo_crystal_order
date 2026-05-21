#include "logger.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <filesystem>

// streambuf that forwards output to two streambufs (like tee)
struct dualbuf : public std::streambuf {
    std::streambuf* sb1;
    std::streambuf* sb2;
    dualbuf(std::streambuf* a, std::streambuf* b) : sb1(a), sb2(b) {}
    int overflow(int c) override {
        if (c == EOF) return !EOF;
        const int r1 = sb1->sputc(c);
        const int r2 = sb2->sputc(c);
        return (r1 == EOF || r2 == EOF) ? EOF : c;
    }
    std::streamsize xsputn(const char* s, std::streamsize n) override {
        auto w1 = sb1->sputn(s, n);
        auto w2 = sb2->sputn(s, n);
        return (w1 == n && w2 == n) ? n : std::min(w1, w2);
    }
    int sync() override {
        const int r1 = sb1->pubsync();
        const int r2 = sb2->pubsync();
        return (r1 == 0 && r2 == 0) ? 0 : -1;
    }
};

static std::ofstream g_logFile;
static std::streambuf* g_origCout = nullptr;
static std::streambuf* g_origCerr = nullptr;
static dualbuf* g_dualCout = nullptr;
static dualbuf* g_dualCerr = nullptr;

namespace {
    // no helpers here; use slog::makeUnique for consistent policy
}

namespace slog {

bool start(const std::string &filename) {
    // pick a unique filename (filename or filename(n)) and open it. Do NOT
    // rename or move existing files; leave them untouched.
    std::string chosen = makeUnique(filename);
    g_logFile.open(chosen, std::ios::out);
    if (!g_logFile) return false;

    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    g_logFile << "=== Run: " << std::put_time(std::localtime(&now_c), "%F %T") << " ===\n";

    g_origCout = std::cout.rdbuf();
    g_origCerr = std::cerr.rdbuf();
    g_dualCout = new dualbuf(g_origCout, g_logFile.rdbuf());
    g_dualCerr = new dualbuf(g_origCerr, g_logFile.rdbuf());
    std::cout.rdbuf(g_dualCout);
    std::cerr.rdbuf(g_dualCerr);
    return true;
}

void stop() {
    if (g_origCout) std::cout.rdbuf(g_origCout);
    if (g_origCerr) std::cerr.rdbuf(g_origCerr);
    if (g_dualCout) { delete g_dualCout; g_dualCout = nullptr; }
    if (g_dualCerr) { delete g_dualCerr; g_dualCerr = nullptr; }
    if (g_logFile) {
        g_logFile << "=== End Run ===\n\n";
        g_logFile.flush();
        g_logFile.close();
    }
}

} // namespace slog

// Public utility to compute a unique filename without renaming anything.
std::string slog::makeUnique(const std::string &path) {
    namespace fs = std::filesystem;
    try {
        if (!fs::exists(path)) return path;
        int n = 1;
        std::string newName;
        do {
            auto pos = path.find_last_of('.');
            std::string prefix = path.substr(0, pos);
            std::string suffix = path.substr(pos);
            newName = prefix + "(" + std::to_string(n) + ")" + suffix;
            ++n;
        } while (fs::exists(newName));
        return newName;
    } catch (...) {
        return path; // on error, fall back to original
    }
}
