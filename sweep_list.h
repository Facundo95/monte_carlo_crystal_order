#ifndef SWEEP_LIST_H
#define SWEEP_LIST_H

#include <vector>
#include <stdexcept>

/**
 * @brief Generic template to create a sweep list for any parameter.
 * @tparam T The numeric type (double, int, etc.)
 * @param start Starting value
 * @param end Ending value
 * @param step Step size (must be positive)
 * @param loop If true, sweep goes from start to end and back. If false, just start to end.
 * @return Vector containing sweep values
 * @throws std::invalid_argument if step <= 0
 */
template<typename T>
std::vector<T> createSweepList(T start, T end, T step, bool loop = false) {
    std::vector<T> list;
    
    if (step <= 0) {
        throw std::invalid_argument("Step must be positive.");
    }
    
    if (start == end) {
        list.push_back(start);
        return list;
    }

    if (!loop) {
        // Non-looping case: go from start to end only
        if (start < end) {
            for (T val = start; val <= end; val += step) {
                list.push_back(val);
            }
        } else {
            for (T val = start; val >= end; val -= step) {
                list.push_back(val);
            }
        }
        return list;
    }
    
    // Looping case: go from start to end, then back to start
    if (start < end) {
        for (T val = start; val <= end; val += step) {
            list.push_back(val);
        }
        for (T val = end - step; val >= start; val -= step) {
            list.push_back(val);
        }
    } else {
        for (T val = start; val >= end; val -= step) {
            list.push_back(val);
        }
        for (T val = end + step; val <= start; val += step) {
            list.push_back(val);
        }
    }
    return list;
}

#endif // SWEEP_LIST_H
