#include "rng.h"

// Constants defined in the original code
#define semilla 0x00001158e460913dL
#define mascara48 0x00001ffffffffffL

// Static variables used by the generator (must be defined once)
static unsigned long ALE32 = semilla;
static union
{
    long tmp;
    float fl;
} xx;

/** @brief Returns a random integer. */
long RanEnt(){
    ALE32 = ALE32 * semilla;
    ALE32 = ALE32 & mascara48;
    return((long) ALE32);
}

/** @brief Returns a random float in [0, 1). */
float Ran0a1(){
    xx.tmp = RanEnt();
    // Masking bits to ensure result is between 0 and 1 (standard LCG trick)
    xx.tmp = (xx.tmp & 0x3fffffffL) | 0x3f800000l;
    return((double)(xx.fl - 1.0));
}

/** @brief Returns a random integer in [1, 8]. */
int RanEnt1a8(){
    // Use the random float generation to map to an integer between 1 and 8
    xx.fl = RanEnt();
    xx.fl = (xx.tmp & 0x00000007L);
    return((int)(xx.fl + 1.0));
}
