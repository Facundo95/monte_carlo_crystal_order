#ifndef RNG_H
#define RNG_H

#include "heisenberg_hamiltonian.h"

/** @brief Returns a random integer. */
long RanEnt();

/** @brief Returns a random float in [0, 1). */
float Ran0a1();

/** @brief Returns a random integer in [1, 8]. */
int RanEnt1a8();

/** @brief Returns a uniformly distributed random unit vector. */
HeisenbergVector randomUnitVector();

/** @brief Returns a three-component standard Gaussian random vector. */
HeisenbergVector randomGaussianVector();

/**
 * @brief Generates a normalized Gaussian trial state around a unit moment.
 * @param moment The current unit vector.
 * @param sigma Gaussian move width.
 */
HeisenbergVector gaussianSpinProposal(const HeisenbergVector& moment,
									  double sigma);

#endif // RNG_H