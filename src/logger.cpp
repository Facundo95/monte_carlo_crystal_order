#include "logger.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <algorithm>

namespace {
    bool fileExists(const std::string &path) {
        std::ifstream file(path.c_str(), std::ios::in);
        return file.good();
    }
}

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
    std::tm now_tm = *std::localtime(&now_c);
    char timestamp[32];
    std::strftime(timestamp, sizeof(timestamp), "%F %T", &now_tm);
    g_logFile << "=== Run: " << timestamp << " ===\n";

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
    if (!fileExists(path)) {
        return path;
    }

    std::string prefix = path;
    std::string suffix;
    const std::string::size_type pos = path.find_last_of('.');

    if (pos != std::string::npos && pos != 0) {
        prefix = path.substr(0, pos);
        suffix = path.substr(pos);
    }

    for (int n = 1; ; ++n) {
        const std::string candidate = prefix + "(" + std::to_string(n) + ")" + suffix;
        if (!fileExists(candidate)) {
            return candidate;
        }
    }
}
