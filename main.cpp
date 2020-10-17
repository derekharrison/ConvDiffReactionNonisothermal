/*
 * main.cpp
 *
 *  Created on: Oct 13, 2020
 *      Author: d-w-h
 *
 *      This code solves the coupled convection diffusion and energy equations:
 *      Da/U*d2Ca/dz2 - dCa/dz + ra/U = 0
 *      Da/U*d2Ca/dz2 - dCa/dz + ra/U = 0
 *      gamma/(U*rho*Cp)*d2T/dz2 - dT/dz + ra*Ha/(U*rho*Cp) = 0
 *      Using the Gauss-Seidel method.
 */

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "solver.hpp"
#include "user_types.hpp"

double ra(double Ca, double Cb, double T) {
    /* Reaction rate law component a */
    double k = 1.32e+22;
    return -k*exp(-14017/T)*Ca*Cb*Cb*Ca*Ca;
}

double rb(double Ca, double Cb, double T) {
    /* Reaction rate law component b */
    double k = 1.32e+22;
    return -k*2*exp(-14017/T)*Ca*Cb*Cb*Ca*Ca;
}

int main(int argc, char* argv[]) {
    p_params physical_parameters;
    g_params grid_parameters;
    s_data solver_data;
    clock_t start_time = clock();

    /* Parameters */
    grid_parameters.num_nodes = 160;      //Number of nodes
    grid_parameters.L = 4.0;              //Length of domain
    physical_parameters.U = 2.0;          //Fluid velocity
    physical_parameters.Cao = 1.1;        //Inlet concentration component a
    physical_parameters.Cbo = 2.9;        //Inlet concentration component b
    physical_parameters.To = 273.0;       //Inlet temperature
    physical_parameters.Da = 1.0;         //Diffusion coefficient
    physical_parameters.rho = 1.0;        //Density
    physical_parameters.Cp = 4.0;         //Heat capacity
    physical_parameters.gamma = 1.0;      //Heat conduction coefficient
    physical_parameters.Ha = -20.0;      //Reaction enthalpy

    /* Allocate data for solver results */
    solver_data.Ca = new double[grid_parameters.num_nodes];
    solver_data.Cb = new double[grid_parameters.num_nodes];
    solver_data.T = new double[grid_parameters.num_nodes];
    solver_data.z_c = new double[grid_parameters.num_nodes];

    /* Execute solver */
    solver(physical_parameters, grid_parameters, &solver_data);

    /* Print data */
    for(int i = 0; i < grid_parameters.num_nodes; ++i) {
        printf("i: %i, z: %f, Ca: %f, Cb: %f, T: %f\n", i, solver_data.z_c[i], solver_data.Ca[i], solver_data.Cb[i], solver_data.T[i]);
    }

    clock_t end_time = clock();
    double execution_time = (double) (end_time - start_time)/CLOCKS_PER_SEC;
    printf("execution time: %f\n", execution_time);

    /* Deallocate data */
    delete [] solver_data.Ca;
    delete [] solver_data.Cb;
    delete [] solver_data.T;
    delete [] solver_data.z_c;

    return 0;
}

