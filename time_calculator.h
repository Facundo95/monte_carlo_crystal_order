#ifndef TIME_CALCULATOR_H
#define TIME_CALCULATOR_H

#include <chrono>
#include <iostream>
#include <iomanip>

class TimeCalculator {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> endTime;
    double durationSeconds;
 
public:
    TimeCalculator() : durationSeconds(0.0) {}

    void start() {
        startTime = std::chrono::high_resolution_clock::now();
    }
 
    void stop() {
        endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = endTime - startTime;
        durationSeconds = duration.count();
    }
    
    double getDurationSeconds() const {
        return durationSeconds;
    }
    
    void displayDuration() const {
        int hours = static_cast<int>(durationSeconds / 3600);
        int minutes = static_cast<int>((durationSeconds - hours * 3600) / 60);
        int seconds = static_cast<int>(durationSeconds - hours * 3600 - minutes * 60);
        std::cout << "Total time: "
        << std::setfill('0') << std::setw(2) << hours << ":"
        << std::setfill('0') << std::setw(2) << minutes << ":"
        << std::setfill('0') << std::setw(2) << seconds << std::endl;
    }
};

#endif // TIME_CALCULATOR_H