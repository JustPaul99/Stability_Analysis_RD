function [theta, conf, StabilityMatrix, eigenvalue_function, k_max, max_eigenvalue] = analyze_data(u, v, t, L, tbegin, tend, tstepsize, xbegin, xend, stepsize)
    %% --- Step 1: Define analysis resolution based on user-defined sampling ---
    tteststeps = floor(1 + (tend - tbegin) / tstepsize);     % Number of time samples
    dttest = t / (tteststeps - 1);                           % Effective time step size

    xteststeps = floor(1 + (xend - xbegin) / stepsize);      % Number of spatial samples
    dxtest = L / (xteststeps - 1);                           % Effective space step size

    %% --- Step 2: Extract and center data ---
    % If you prefer to subtract a fixed reference (e.g., steady state), use:
    % U = transpose(u(tbegin:tstepsize:tend, xbegin:stepsize:xend)) - ustar;
    % V = transpose(v(tbegin:tstepsize:tend, xbegin:stepsize:xend)) - vstar;

    % Current implementation subtracts the mean of the sampled data
    U = transpose(u(tbegin:tstepsize:tend, xbegin:stepsize:xend)) - mean(transpose(u(tbegin:tstepsize:tend, xbegin:stepsize:xend)));
    V = transpose(v(tbegin:tstepsize:tend, xbegin:stepsize:xend)) - mean(transpose(v(tbegin:tstepsize:tend, xbegin:stepsize:xend)));

    %% --- Step 3: Initialize derivative matrices ---
    A = zeros(xteststeps, tteststeps, 2);  % Stores time derivatives of U and V
    B = zeros(xteststeps, tteststeps, 4);  % Stores U, V, and second spatial derivatives

    %% --- Step 4: Compute time and spatial derivatives ---
    % Assign values to matrices as defined in your reference paper
    for j = 1:tteststeps-1
        % Time derivatives (forward difference)
        A(:, j, 1) = (U(:, j+1) - U(:, j)) / dttest;
        A(:, j, 2) = (V(:, j+1) - V(:, j)) / dttest;

        % Store U and V values
        B(:, j, 1) = U(:, j);
        B(:, j, 2) = V(:, j);

        % Second-order spatial derivatives (central difference)
        B(1, j, 3) = 2 * (U(2, j) - U(1, j)) / dxtest^2;
        B(xteststeps, j, 3) =  (U(xteststeps-1, j) - U(xteststeps, j)) / dxtest^2;
        B(1, j, 4) = 2 * (V(2, j) - V(1, j)) / dxtest^2;
        B(xteststeps, j, 4) = (V(xteststeps-1, j) - V(xteststeps, j)) / dxtest^2;

        % Interior points: second-order central difference
        B(2:end-1, j, 3) = (U(3:end, j) + U(1:end-2, j) - 2 * U(2:end-1, j)) / dxtest^2;
        B(2:end-1, j, 4) = (V(3:end, j) + V(1:end-2, j) - 2 * V(2:end-1, j)) / dxtest^2;
    end

    %% --- Step 5: Reshape matrices for least squares fitting ---
    AcheckU = reshape(A(:, :, 1), xteststeps * tteststeps, 1);  % Flatten U time derivatives
    AcheckV = reshape(A(:, :, 2), xteststeps * tteststeps, 1);  % Flatten V time derivatives
    Bmatrix = reshape(B, xteststeps * tteststeps, 4);           % Flatten B matrix

    % Construct block matrix for fitting both U and V equations
    Bmatrix2 = sparse([Bmatrix, zeros(length(Bmatrix), 4); zeros(length(Bmatrix), 4), Bmatrix]);
    Acheck = [AcheckU; AcheckV];                                % Combine U and V derivatives

    %% --- Step 6: Define constraints for least squares optimization ---
    ineq = [0, 0, 0, 0, 0, 0, 0, -1];      % Inequality constraint (e.g., damping condition)
    bineq = 0;
    eq = [0, 0, 0, 1, 0, 0, 0, 0;          % Enforce specific parameter values
          0, 0, 1, 0, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0, 1, 0];
    beq = [0, 1, 0];                       % Corresponding equality targets
    lb = -200 * ones(8, 1);               % Lower bounds
    ub = 200 * ones(8, 1);                % Upper bounds
    x0 = [1, 1, 1, 0, 1, 1, 0, 1];        % Initial guess for parameters

    %% --- Step 7: Perform constrained least squares fitting ---
    fun = @(x, xdata) xdata * transpose(x);  % Linear model function
    [theta, resnorm, residual, exitflag, output, lambda, jacobian] = ...
        lsqcurvefit(fun, x0, Bmatrix2, Acheck, lb, ub, ineq, bineq, eq, beq);

    display(theta)  % Display fitted parameters

    %% --- Step 8: Compute confidence intervals (optional) ---
    % If you don't want to compute confidence intervals, you can set:
    % conf = 0;
    conf = nlparci(theta, residual, "Jacobian", jacobian);  % Confidence intervals (ignores variable interactions)
    % display(conf)

    %% --- Step 9: Stability analysis ---
    % Define stability matrix as a function of wavenumber k
    StabilityMatrix = @(k1) [theta(1) - theta(3) * k1.^2, theta(2);
                             theta(5), theta(6) - k1.^2 * theta(8)];

    % Define function to compute dominant eigenvalue
    eigenvalue_function = @(k1) max(eig(StabilityMatrix(k1)));

    % Find wavenumber k that maximizes the dominant eigenvalue
    k_max = fminbnd(@(k1) -eigenvalue_function(k1), 0, 100);
    max_eigenvalue = eigenvalue_function(k_max);

    %% --- Optional: Print stability result ---
    % if max_eigenvalue > 0
    %     fprintf('For k=%.2f we have that λ=%.2f and thus our solution is unstable.\n', k_max, max_eigenvalue);
    % else
    %     fprintf('For k=%.2f we have that λ=%.2f and thus our solution is stable.\n', k_max, max_eigenvalue);
    % end
end
