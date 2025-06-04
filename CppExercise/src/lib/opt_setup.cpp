#include <rlab/controller.h>

using namespace casadi;

OptimizationSetup::OptimizationSetup(int seed) : seed(seed) {
    setupCostMatrices();
    setupNormalization();
    setupConstraints();
}

void OptimizationSetup::setupCostMatrices() {
    DM Q_weights = {W_x, W_z, W_theta, W_vx, W_vz, W_wy};
    DM R_weights = {W_T, W_gim};
    DM D_weights = {W_T_rate, W_gim_rate};
    DM M_weights = {W_H};

    Q = diag(Q_weights);
    R = diag(R_weights);
    D = diag(D_weights);
    M = diag(M_weights);
}

void OptimizationSetup::setupNormalization() {
    X_max = DM({x_max, z_max, theta_max, vx_max, vz_max, wy_max});
    DM X_norm = diag(1 / pow(X_max, 2));
    Q_norm = mtimes(Q, X_norm);

    U_max = DM({T_max, gim_max});
    DM U_norm = diag(1 / pow(U_max, 2) / N_steps);
    R_norm = mtimes(R, U_norm);

    U_rate_max = DM({throt_rate, gim_rate});
    DM U_rate_norm = diag(1 / pow(U_rate_max, 2) / N_steps);
    D_norm = mtimes(D, U_rate_norm);

    H_max_dm = DM({H_0});
    DM H_norm = diag(1 / pow(H_max_dm, 2));
    M_norm = mtimes(M, H_norm);
}

void OptimizationSetup::setupConstraints() {
    generateInitialConditions();

    X_f = DM({x_f, h_f, theta_f, vx_f, vz_f, wy_f});
    slack = DM({x_e, h_e, theta_e, vx_e, vz_e, wy_e});

    lb_X = DM({x_lb, z_lb, theta_lb, vx_lb, vz_lb, wy_lb});
    ub_X = DM({x_ub, z_ub, theta_ub, vx_ub, vz_ub, wy_ub});

    lb_U = DM({T_min, gim_min});
    ub_U = DM({T_max, gim_max});
}

void OptimizationSetup::generateInitialConditions() {
    // Random initial conditions from dispersion parameters
    std::default_random_engine gen;
    if (seed >= 0) {
        gen = std::default_random_engine(seed);
    } else {
        std::random_device rd;
        gen = std::default_random_engine(rd());
    }
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    double x_r = x_d * dist(gen);
    double h_r = h_d * dist(gen);
    double theta_r = theta_d * dist(gen);
    double vx_r = vx_d * dist(gen);
    double vz_r = vz_d * dist(gen);
    double wy_r = wy_d * dist(gen);

    X_i = DM({x_i + x_r, h_i + h_r, theta_i + theta_r,
              vx_i + vx_r, vz_i + vz_r, wy_i + wy_r});
}

DMDict OptimizationSetup::getCostMatrices() const {
    return {
        {"Q", Q_norm},
        {"R", R_norm},
        {"D", D_norm},
        {"M", M_norm}
    };
}

DMDict OptimizationSetup::getBounds() const {
    return {
        {"lb_X", lb_X},
        {"ub_X", ub_X},
        {"lb_U", lb_U},
        {"ub_U", ub_U},
        {"X_max", X_max},
        {"U_max", U_max},
        {"U_rate_max", U_rate_max}
    };
}

DMDict OptimizationSetup::getConditions() const {
    return {
        {"X_f", X_f},
        {"slack", slack},
        {"X_i", X_i}
    };
}