// Optimizer.cpp
#include <rlab/controller.h>

using namespace casadi;

Optimizer::Optimizer(const Dynamics& dyn, const OptimizationSetup& setup, const std::string& guess)
    : guess(guess) {

    N_X = dyn.getStateDim();
    N_U = dyn.getInputDim();
    f   = dyn.getContinuousDynamics();
    F   = dyn.getDiscreteDynamics();
    gim_dyn = (N_X == 8); // Check if gimbal dynamics are included

    DMDict cost_matrices = setup.getCostMatrices();
    Q = cost_matrices["Q"];
    R = cost_matrices["R"];
    D = cost_matrices["D"];
    M = cost_matrices["M"];

    DMDict conditions = setup.getConditions();
    X_i = conditions["X_i"];
    X_f = conditions["X_f"];
    slack = conditions["slack"];

    DMDict bounds = setup.getBounds();
    U_max = bounds["U_max"];
    lb_X = bounds["lb_X"];
    ub_X = bounds["ub_X"];
    lb_U = bounds["lb_U"];
    ub_U = bounds["ub_U"];
    X_max = bounds["X_max"];

    opti = Opti();

    DM U_mat = repmat(U_max, 1, N_steps);
    U = U_mat * opti.variable(N_U, N_steps);
    dt = opti.variable(N_steps);
    H = sum1(dt);

    if (gim_dyn) {
        X_max = vertcat(X_max, gim_max, gim_rate);
        ub_X = vertcat(ub_X, gim_max, gim_rate);
        lb_X = vertcat(lb_X, -gim_max, -gim_rate);
    }

    DM X_mat = repmat(X_max, 1, N_steps + 1);
    X_all = X_mat * opti.variable(N_X, N_steps + 1);

    if (gim_dyn)
        X = X_all(Slice(0, N_X - 2), Slice());
    else
        X = X_all;

    setupProblem();
}

void Optimizer::setupProblem() {

    opti.subject_to(opti.bounded(lb_U, U, ub_U));
    opti.subject_to(opti.bounded(lb_X, X_all, ub_X));
    opti.subject_to(opti.bounded(lb_dt, dt, ub_dt));

    opti.subject_to(X(Slice(), 0) == X_i);
    opti.subject_to(opti.bounded(-slack, X(Slice(), N_steps) - X_f, slack));

    MX objective = 0;
    if (gim_dyn) {
        MX T_max_vec = dt * throt_rate;
        MX U_rate;
        objective += bilin(R, U(Slice(), 0)); // Control cost
        for (int k = 0; k < N_steps-1; ++k) {
            objective += bilin(R, U(Slice(), k+1)); // Control cost
            U_rate = (U(Slice(), k + 1) - U(Slice(), k)) / dt(k); // Control rate
            objective += bilin(D, U_rate); // Control rate cost
            opti.subject_to(
            opti.bounded(U(0, k) - T_max_vec(k), U(0, k + 1), T_max_vec(k) + U(0, k))); // Throttle rate constraint
        }
    } else {
        for (int k = 0; k < N_steps; ++k) {
            objective += bilin(R, U(Slice(), k)); // Control cost
        }    
    }

    MX X_kp1 = F(MXDict{{"x0", X_all(Slice(),Slice(0, N_steps))}, {"u", U}, {"p", dt}}).at("xf");
    opti.subject_to(X_all(Slice(),Slice(1, N_steps+1)) == X_kp1); // Enforce the discretized dynamics
    objective += bilin(Q, X(Slice(), N_steps)); // State cost
    objective += bilin(M, H);

    setInitialGuess();

    opti.minimize(objective);
}


void Optimizer::setInitialGuess() {
    try {        
        std::ifstream file(guess);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open initial guess file: " + guess);
        }
        
        std::vector<std::vector<double>> data;
        std::string line;
        
        // Skip header line
        if (!std::getline(file, line)) {
            throw std::runtime_error("Empty CSV file or could not read header");
        }        
        // Read all data lines
        int line_number = 1;
        while (std::getline(file, line)) {
            line_number++;
            
            // Skip empty lines
            if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos) {
                continue;
            }
            
            std::vector<double> row;
            std::stringstream ss(line);
            std::string cell;
            
            // Parse each cell in the row
            while (std::getline(ss, cell, ',')) {
                // Trim whitespace
                cell.erase(0, cell.find_first_not_of(" \t\r\n"));
                cell.erase(cell.find_last_not_of(" \t\r\n") + 1);
                
                // Handle empty cells
                if (cell.empty()) {
                    row.push_back(0.0);
                } else {
                    try {
                        double value = std::stod(cell);
                        row.push_back(value);
                    } catch (const std::exception& e) {
                        std::cerr << "Warning: Could not parse '" << cell 
                                 << "' as number at line " << line_number 
                                 << ". Using 0.0" << std::endl;
                        row.push_back(0.0);
                    }
                }
            }
            
            if (!row.empty()) {
                data.push_back(row);
            }
        }
        
        file.close();
        
        if (data.empty()) {
            throw std::runtime_error("No valid data rows found in CSV");
        }        
        // Determine expected number of columns
        int expected_cols = 1 + N_X + N_U;  // dt + state + control
        
        // Initialize matrices
        DM X_guess = DM::zeros(N_X, N_steps + 1);
        DM U_guess = DM::zeros(N_U, N_steps);
        DM dt_guess = DM::zeros(1, N_steps);
        
        // Process data rows
        int steps_processed = 0;
        for (size_t i = 0; i < data.size() && steps_processed < N_steps; ++i) {
            const auto& row = data[i];
            
            // Check if row has enough columns
            if (row.size() < expected_cols) {
                std::cerr << "Warning: Row " << i << " has only " << row.size() 
                         << " columns, expected " << expected_cols << ". Skipping." << std::endl;
                continue;
            }
            
            // Extract dt (first column)
            dt_guess(0, steps_processed) = row[0];
            
            // Extract state variables (next N_X columns)
            for (int j = 0; j < N_X; ++j) {
                if (1 + j < row.size()) {
                    X_guess(j, steps_processed) = row[1 + j];
                } else {
                    X_guess(j, steps_processed) = 0.0;
                }
            }
            
            // Extract control variables (next N_U columns)
            for (int j = 0; j < N_U; ++j) {
                if (1 + N_X + j < row.size()) {
                    U_guess(j, steps_processed) = row[1 + N_X + j];
                } else {
                    U_guess(j, steps_processed) = 0.0;
                }
            }
            
            steps_processed++;
        }
        
        // Handle final state X[:, N_steps]
        if (steps_processed > 0) {
            // Use the last processed row's state as the final state
            for (int j = 0; j < N_X; ++j) {
                X_guess(j, N_steps) = X_guess(j, steps_processed - 1);
            }
        }
        
        // Apply initial guesses to the optimizer
        opti.set_initial(X_all, X_guess);
        opti.set_initial(U, U_guess);
        opti.set_initial(dt, dt_guess);
        
        std::cout << "Initial guess set successfully" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error reading initial guess: " << e.what() << std::endl;
        std::cerr << "Using default initial values instead." << std::endl;
        
        // Fallback to default values
        opti.set_initial(X_all, DM::zeros(N_X, N_steps + 1));
        opti.set_initial(U, DM::zeros(N_U, N_steps));
        opti.set_initial(dt, DM::ones(1, N_steps) * (H_0 / N_steps));
    }
}


DMDict Optimizer::solve(bool console_print, bool plot, bool save_traj, bool save_plots) {

    Dict p_opts{{"expand", true}};
    Dict s_opts{
        {"print_level", 5}, 
        {"max_iter", 10000}
    };

    opti.solver("ipopt", p_opts, s_opts);

    DM X_opt, U_opt, dt_opt, H_opt;

    try {
        OptiSol sol = opti.solve();
        X_opt = opti.value(X_all);
        U_opt = opti.value(U);
        dt_opt = opti.value(dt);
        H_opt = opti.value(H);
    } catch (...) {
        std::cerr << "Solver failed. Returning debug values." << std::endl;
        X_opt = opti.debug().value(X_all);
        U_opt = opti.debug().value(U);
        dt_opt = opti.debug().value(dt);
        H_opt = opti.debug().value(H);
    }

    DMDict res{
        {"X", X_opt},
        {"U", U_opt},
        {"dt", dt_opt},
        {"H", H_opt}
    };


    if (console_print) {
        std::cout << "\nManeuver duration: " << std::fixed << std::setprecision(2) << H_opt << " s\n";
        std::cout << "\nFinal State:\n";
        std::cout << "x: " << std::fixed << std::setprecision(2) << X_opt(0, N_steps) << "m, ";
        std::cout << "z: " << std::fixed << std::setprecision(2) << X_opt(1, N_steps) << "m, ";
        std::cout << "theta: " << std::fixed << std::setprecision(2) << X_opt(2, N_steps) * RAD << "deg, ";
        std::cout << "vx: " << std::fixed << std::setprecision(2) << X_opt(3, N_steps) << "m/s, ";
        std::cout << "vz: " << std::fixed << std::setprecision(2) << X_opt(4, N_steps) << "m/s, ";
        std::cout << "wy: " << std::fixed << std::setprecision(2) << X_opt(5, N_steps) * RAD << "deg/s\n";
    }

    if (plot) {
        save_traj = true; // Ensure trajectory is saved if plotting is requested
    }

    if (save_traj) {
        std::ofstream out("../output/results.csv");
        if (N_X == 6)
            out << "dt,x,z,theta,vx,vz,wy,T,gimbal\n";
        else
            out << "dt,x,z,theta,vx,vz,wy,gim_real,gim_rate,T,gimbal\n";
        for (int i = 0; i < N_steps; ++i) {
            out << dt_opt(i);
            for (int j = 0; j < N_X; ++j)
                out << "," << X_opt(j, i);
            for (int j = 0; j < N_U; ++j)
                out << "," << U_opt(j, i);
            out << "\n";
        }
        // Write the final state        
        for (int j = 0; j < N_X; ++j)
            out << "," << X_opt(j, N_steps);
        out << "," << ",";
        out.close();
    }

    if (plot) {

        std::string filename = "results.csv";
        std::string title = "Trajectory Optimization Results";
        // Convert booleans to int-like strings
        std::string cmd = "python3 ../src/lib/plotter.py " + filename + " " +
                    std::to_string(N_steps) + " " +
                    std::to_string(N_X) + " \"" +
                    title + "\" " +
                    std::to_string(gim_dyn) + " " +
                    std::to_string(save_plots);
        std::system(cmd.c_str());
    }

    return res;
}
