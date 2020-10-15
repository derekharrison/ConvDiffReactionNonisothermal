/*
 * user_types.hpp
 *
 *  Created on: Oct 13, 2020
 *      Author: d-w-h
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

typedef struct physical_params {
    double U;
    double Cao;
    double Cbo;
    double To;
    double Da;
    double rho;
    double Cp;
    double gamma;
    double Ha;
} p_params;

typedef struct grid_params {
    int num_nodes;
    double L;
} g_params;

typedef struct solver_data {
    double* Ca;
    double* Cb;
    double* T;
    double* z_c;
} s_data;

#endif /* USER_TYPES_HPP_ */
