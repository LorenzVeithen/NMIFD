//
// Created by Lorenz Veithen on 17/10/2024.
//

#include "linAlg.h"

/**
 * Misc functions
 */
void crossProduct(const std::vector<double> &v_A,
                  const std::vector<double> &v_B,
                  std::vector<double> &c_P) {

    c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
    c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
    c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}
