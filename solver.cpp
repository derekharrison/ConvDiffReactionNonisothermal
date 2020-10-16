/*
 * solver.cpp
 *
 *  Created on: Oct 13, 2020
 *      Author: d-w-h
 */

#include <math.h>
#include <stdio.h>
#include "main.hpp"
#include "user_types.hpp"

double dra_dCa(double Ca, double Cb, double T) {
    double dCa = 0.0001;
    return (ra(Ca + dCa, Cb, T) - ra(Ca, Cb, T)) / dCa;
}

double dra_dCb(double Ca, double Cb, double T) {
    double dCb = 0.0001;
    return (ra(Ca, Cb + dCb, T) - ra(Ca, Cb, T)) / dCb;
}

double dra_dT(double Ca, double Cb, double T) {
    double dT = 0.0001;
    return (ra(Ca, Cb, T + dT) - ra(Ca, Cb, T)) / dT;
}

double drb_dCa(double Ca, double Cb, double T) {
    double dCa = 0.0001;
    return (rb(Ca + dCa, Cb, T) - rb(Ca, Cb, T)) / dCa;
}

double drb_dCb(double Ca, double Cb, double T) {
    double dCb = 0.0001;
    return (rb(Ca, Cb + dCb, T) - rb(Ca, Cb, T)) / dCb;
}

double drb_dT(double Ca, double Cb, double T) {
    double dT = 0.0001;
    return (rb(Ca, Cb, T + dT) - rb(Ca, Cb, T)) / dT;
}

void solver(p_params physical_parameters, g_params grid_parameters, s_data* solver_data) {
    double L, U, Cain, Cbin, Tin, rho, Cp, gamma, Cao, Cbo, To, Da, Ha, del_z;
    double error_Ca, error_Cb, error_T;
    int num_nodes, i, outer_it, inner_it, max_outer_it, max_inner_it;

    /* Parameters */
    num_nodes = grid_parameters.num_nodes;
    L = grid_parameters.L;
    U = physical_parameters.U;
    Cao = physical_parameters.Cao;
    Cbo = physical_parameters.Cbo;
    To = physical_parameters.To;
    Da = physical_parameters.Da;
    rho = physical_parameters.rho;
    Cp = physical_parameters.Cp;
    gamma = physical_parameters.gamma;
    Ha = physical_parameters.Ha;

    max_outer_it = 300;
    max_inner_it = 300;

    /* Start simulation */
    del_z = L / num_nodes;
    double* Ca_prev_it = new double[num_nodes];
    double* Cb_prev_it = new double[num_nodes];
    double* T_prev_it = new double[num_nodes];

    /* Initialize Ca and z_c*/
    for(i = 0; i < num_nodes; ++i) {
        solver_data->Ca[i] = 0.0;
        solver_data->Cb[i] = 0.0;
        solver_data->T[i] = 0.0;
        solver_data->z_c[i] = i*del_z + 0.5*del_z;
    }

    /* Gauss Seidel iterations */
    outer_it = 0;
    while(outer_it < max_outer_it) {
        for(i = 0; i < num_nodes; ++i) {
            Ca_prev_it[i] = solver_data->Ca[i];
            Cb_prev_it[i] = solver_data->Cb[i];
            T_prev_it[i] = solver_data->T[i];
        }
        inner_it = 0;
        while(inner_it < max_inner_it) {
            /* Component A */
            //Inlet node
            Cain = (Cao + Da*solver_data->Ca[0]/(U*0.5*del_z)) / (1+Da/(U*0.5*del_z));
            // Left most node
            solver_data->Ca[0] = (Da*Cain/(0.5*del_z) + Da*solver_data->Ca[1]/del_z + U*Cain +
                                  ra(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*del_z -
                                  dra_dCa(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*Ca_prev_it[0]*del_z +
                                  dra_dCb(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*(solver_data->Cb[0] - Cb_prev_it[0])*del_z +
                                  dra_dT(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*(solver_data->T[0] - T_prev_it[0])*del_z) /
                                 (Da/(0.5*del_z) + Da/del_z + U - dra_dCa(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*del_z);
            // Central nodes
            for(i = 1; i < num_nodes - 1; ++i) {
                solver_data->Ca[i] = (Da*solver_data->Ca[i-1]/del_z + Da*solver_data->Ca[i+1]/del_z + U*solver_data->Ca[i-1] +
                                      ra(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*del_z -
                                      dra_dCa(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*Ca_prev_it[i]*del_z +
                                      dra_dCb(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*(solver_data->Cb[i] - Cb_prev_it[i])*del_z +
                                      dra_dT(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*(solver_data->T[i] - T_prev_it[i])*del_z) /
                                     (Da/del_z + Da/del_z + U - dra_dCa(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*del_z);
            }
            // Right most node
            solver_data->Ca[num_nodes-1] = (Da*solver_data->Ca[num_nodes-2]/del_z + U*solver_data->Ca[num_nodes-2] +
                                            ra(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*del_z -
                                            dra_dCa(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*Ca_prev_it[num_nodes-1]*del_z +
                                            dra_dCb(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*(solver_data->Cb[num_nodes-1] - Cb_prev_it[num_nodes-1])*del_z +
                                            dra_dT(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*(solver_data->T[num_nodes-1] - T_prev_it[num_nodes-1])*del_z) /
                                           (Da/del_z + U - dra_dCa(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*del_z);

            /* Component B */
            //Inlet node
            Cbin = (Cbo + Da*solver_data->Cb[0]/(U*0.5*del_z)) / (1+Da/(U*0.5*del_z));
            // Left most node
            solver_data->Cb[0] = (Da*Cbin/(0.5*del_z) + Da*solver_data->Cb[1]/del_z + U*Cbin +
                                  rb(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*del_z -
                                  drb_dCb(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*Cb_prev_it[0]*del_z +
                                  drb_dCa(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*(solver_data->Ca[0] - Ca_prev_it[0])*del_z +
                                  drb_dT(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*(solver_data->T[0] - T_prev_it[0])*del_z) /
                                 (Da/(0.5*del_z) + Da/del_z + U - drb_dCb(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*del_z);
            // Central nodes
            for(i = 1; i < num_nodes - 1; ++i) {
                solver_data->Cb[i] = (Da*solver_data->Cb[i-1]/del_z + Da*solver_data->Cb[i+1]/del_z + U*solver_data->Cb[i-1] +
                                      rb(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*del_z -
                                      drb_dCb(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*Cb_prev_it[i]*del_z +
                                      drb_dCa(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*(solver_data->Ca[i] - Ca_prev_it[i])*del_z +
                                      drb_dT(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*(solver_data->T[i] - T_prev_it[i])*del_z) /
                                     (Da/del_z + Da/del_z + U - drb_dCb(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*del_z);
            }
            // Right most node
            solver_data->Cb[num_nodes-1] = (Da*solver_data->Cb[num_nodes-2]/del_z + U*solver_data->Cb[num_nodes-2] +
                                            rb(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*del_z -
                                            drb_dCb(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*Cb_prev_it[num_nodes-1]*del_z +
                                            drb_dCa(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*(solver_data->Ca[num_nodes-1] - Ca_prev_it[num_nodes-1])*del_z +
                                            drb_dT(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*(solver_data->T[num_nodes-1] - T_prev_it[num_nodes-1])*del_z) /
                                           (Da/del_z + U - drb_dCb(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*del_z);

            /* T */
            //Inlet node
            Tin = (To + gamma*solver_data->T[0]/(U*rho*Cp*0.5*del_z)) / (1+gamma/(U*rho*Cp*0.5*del_z));
            // Left most node
            solver_data->T[0] = (gamma*Tin/(0.5*del_z) + gamma*solver_data->T[1]/del_z + U*rho*Cp*Tin +
                                 ra(Ca_prev_it[0], Cb_prev_it[0], T_prev_it[0])*del_z*Ha) /
                                (gamma/(0.5*del_z) + gamma/del_z + U*rho*Cp);
            // Central nodes
            for(i = 1; i < num_nodes - 1; ++i) {
                solver_data->T[i] = (gamma*solver_data->T[i-1]/del_z + gamma*solver_data->T[i+1]/del_z + U*rho*Cp*solver_data->T[i-1] +
                                     ra(Ca_prev_it[i], Cb_prev_it[i], T_prev_it[i])*del_z*Ha) /
                                    (gamma/del_z + gamma/del_z + U*rho*Cp);
            }
            // Right most node
            solver_data->T[num_nodes-1] = (gamma*solver_data->T[num_nodes-2]/del_z + U*rho*Cp*solver_data->T[num_nodes-2] +
                                           ra(Ca_prev_it[num_nodes-1], Cb_prev_it[num_nodes-1], T_prev_it[num_nodes-1])*del_z*Ha) /
                                          (gamma/del_z + U*rho*Cp);

            ++inner_it;
        }

        ++outer_it;
    }

    /* Check convergence */
    error_Ca = 0.0;
    error_Cb = 0.0;
    error_T = 0.0;
    error_Ca += fabs(1.0 - (-Da*(solver_data->Ca[0] - Cain)/(0.5*del_z) + Da*(solver_data->Ca[1] - solver_data->Ca[0])/del_z + U*Cain - U*solver_data->Ca[0]) / (-ra(solver_data->Ca[0], solver_data->Cb[0], solver_data->T[0])*del_z));
    error_Cb += fabs(1.0 - (-Da*(solver_data->Cb[0] - Cbin)/(0.5*del_z) + Da*(solver_data->Cb[1] - solver_data->Cb[0])/del_z + U*Cbin - U*solver_data->Cb[0]) / (-rb(solver_data->Ca[0], solver_data->Cb[0], solver_data->T[0])*del_z));
    error_T += fabs(1.0 - (-gamma*(solver_data->T[0] - Tin)/(0.5*del_z) + gamma*(solver_data->T[1] - solver_data->T[0])/del_z + U*rho*Cp*Tin - U*rho*Cp*solver_data->T[0]) / (-ra(solver_data->Ca[0], solver_data->Cb[0], solver_data->T[0])*del_z*Ha));
    for(i = 1; i < num_nodes - 1; ++i) {
        error_Ca += fabs(1.0 - (-Da*(solver_data->Ca[i] - solver_data->Ca[i-1])/del_z + Da*(solver_data->Ca[i+1] - solver_data->Ca[i])/del_z + U*solver_data->Ca[i-1] - U*solver_data->Ca[i]) / (-ra(solver_data->Ca[i], solver_data->Cb[i], solver_data->T[i])*del_z));
        error_Cb += fabs(1.0 - (-Da*(solver_data->Cb[i] - solver_data->Cb[i-1])/del_z + Da*(solver_data->Cb[i+1] - solver_data->Cb[i])/del_z + U*solver_data->Cb[i-1] - U*solver_data->Cb[i]) / (-rb(solver_data->Ca[i], solver_data->Cb[i], solver_data->T[i])*del_z));
        error_T += fabs(1.0 - (-gamma*(solver_data->T[i] - solver_data->T[i-1])/del_z + gamma*(solver_data->T[i+1] - solver_data->T[i])/del_z + U*rho*Cp*solver_data->T[i-1] - U*rho*Cp*solver_data->T[i]) / (-ra(solver_data->Ca[i], solver_data->Cb[i], solver_data->T[i])*del_z*Ha));
    }

    error_Ca += fabs(1.0 - (-Da*(solver_data->Ca[num_nodes-1] - solver_data->Ca[num_nodes-2])/(del_z) + U*solver_data->Ca[num_nodes-2] - U*solver_data->Ca[num_nodes-1]) / (-ra(solver_data->Ca[num_nodes-1], solver_data->Cb[num_nodes-1], solver_data->T[num_nodes-1])*del_z));
    error_Cb += fabs(1.0 - (-Da*(solver_data->Cb[num_nodes-1] - solver_data->Cb[num_nodes-2])/(del_z) + U*solver_data->Cb[num_nodes-2] - U*solver_data->Cb[num_nodes-1]) / (-rb(solver_data->Ca[num_nodes-1], solver_data->Cb[num_nodes-1], solver_data->T[num_nodes-1])*del_z));
    error_T +=  fabs(1.0 - (-gamma*(solver_data->T[num_nodes-1] - solver_data->T[num_nodes-2])/(del_z) + U*rho*Cp*solver_data->T[num_nodes-2] - U*rho*Cp*solver_data->T[num_nodes-1]) / (-ra(solver_data->Ca[num_nodes-1], solver_data->Cb[num_nodes-1], solver_data->T[num_nodes-1])*del_z*Ha));

    printf("Convergence Ca: %E\n", error_Ca / num_nodes);
    printf("Convergence Cb: %E\n", error_Cb / num_nodes);
    printf("Convergence T: %E\n", error_T / num_nodes);

    delete [] Ca_prev_it;
    delete [] Cb_prev_it;
    delete [] T_prev_it;
}
