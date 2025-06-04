#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <casadi/casadi.hpp>
#include <rlab/simulation_parameters.h>
#include <rlab/dynamics.h>
#include <string>
#include <random>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace casadi;

class OptimizationSetup {
    // Optimization setup for the controller
    int seed;

    DM Q, R, D, M;
    DM Q_norm, R_norm, D_norm, M_norm;
    DM X_max, U_max, U_rate_max, H_max_dm;
    DM lb_X, ub_X, lb_U, ub_U;
    DM X_f, slack, X_i;

    void setupCostMatrices();
    void setupNormalization();
    void setupConstraints();
    void generateInitialConditions();

public:
    OptimizationSetup(int seed = -1);

    DMDict getCostMatrices() const;
    DMDict getBounds() const;
    DMDict getConditions() const;

};

class Optimizer {

    int N_X, N_U;

    Function f, F;
    DM Q, R, D, M;
    DM X_i, X_f, slack;
    DM lb_X, ub_X, lb_U, ub_U;
    DM U_max, U_rate_max;
    DM X_max;
    Function solver;

    bool implicit, gim_dyn;
    std::string guess;

    Opti opti;
    MX X_all, X, U, dt, H;

public:
    Optimizer(const Dynamics& dyn, const OptimizationSetup& setup, const std::string& guess = "");

    void setupProblem();
    void setInitialGuess();

    DMDict solve(
        bool console_print = false,
        bool plot = false,
        bool save_traj = false,
        bool save_plots = false
    );
};


#endif // CONTROLLER_H