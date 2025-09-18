function [u, v] = simulate_Klausmeier_RK(delta, h, m, pvec, t, tsteps, L, xsteps, noiseType, correlation_length, NoiseScale)
    %% --- Grid setup ---
    dt = t / (tsteps - 1);               % Time step size
    dx = L / (xsteps - 1);               % Spatial step size
    time = 0:dt:t;                       % Time vector
    place = 0:dx:L;                      % Space vector
    [T, X] = meshgrid(time, place);      % Meshgrid for plotting (optional)

    %% --- Compute steady states u* and v* for intial solution---
    ustar = (2*h*m + pvec(1) + 2*h^2*pvec(1) - sqrt(-4*m^2 - 4*h*m*pvec(1) + pvec(1)^2)) / (2 + 2*h^2);
    vstar = (pvec(1) + sqrt(pvec(1)^2 - 4*m*(m + h*pvec(1)))) / (2*m + 2*h*pvec(1));

    %% --- Initialize solution matrices with steady states ---
    u = ones(tsteps, xsteps) * ustar;
    v = ones(tsteps, xsteps) * vstar;
    epsilon = NoiseScale * sqrt(dt);  % Noise scaling factor

    %% --- Time-stepping loop using RK4 ---
    for i = 1:tsteps-1
        %% --- Generate noise based on type ---
         if strcmp(noiseType, 'Correlated')
            x = -xsteps:1:xsteps;
            filter = exp(-x.^2 / (correlation_length^2 / dx^2));
            filter = filter / sqrt(sum(filter.^2));  % Normalize filter

            uncorrelated_noise_u = epsilon * randn(xsteps, 1);
            uncorrelated_noise_v = epsilon * randn(xsteps, 1);

            Noise = conv(uncorrelated_noise_u, filter, 'same');
            Noise2 = conv(uncorrelated_noise_v, filter, 'same');

        elseif strcmp(noiseType, 'White')
            Noise = epsilon * randn(xsteps, 1);
            Noise2 = epsilon * randn(xsteps, 1);

        else  % Uniform noise
            Noise = epsilon * (rand(xsteps, 1) - 0.5) * 2;
            Noise2 = epsilon * (rand(xsteps, 1) - 0.5) * 2;
        end

        %% --- RK4 integration for each spatial point ---
        for j = 1:xsteps
            % RK4 steps for u
            k1u = dt * du_dt(u(i, :), v(i, :), pvec(i), dx, j, xsteps);
            k2u = dt * du_dt(u(i, :) + 0.5 * k1u, v(i, :) + 0.5 * k1u, pvec(i), dx, j, xsteps);
            k3u = dt * du_dt(u(i, :) + 0.5 * k2u, v(i, :) + 0.5 * k2u, pvec(i), dx, j, xsteps);
            k4u = dt * du_dt(u(i, :) + k3u, v(i, :) + k3u, pvec(i), dx, j, xsteps);

            % RK4 steps for v
            k1v = dt * dv_dt(u(i, :), v(i, :), delta, h, m, dx, j, xsteps);
            k2v = dt * dv_dt(u(i, :) + 0.5 * k1v, v(i, :) + 0.5 * k1v, delta, h, m, dx, j, xsteps);
            k3v = dt * dv_dt(u(i, :) + 0.5 * k2v, v(i, :) + 0.5 * k2v, delta, h, m, dx, j, xsteps);
            k4v = dt * dv_dt(u(i, :) + k3v, v(i, :) + k3v, delta, h, m, dx, j, xsteps);

            % Update solution with RK4 and noise
            u(i+1, j) = u(i, j) + (k1u + 2*k2u + 2*k3u + k4u) / 6 + Noise(j);
            v(i+1, j) = v(i, j) + (k1v + 2*k2v + 2*k3v + k4v) / 6 + Noise2(j);
        end
    end
end

%% --- Function to compute du/dt ---
function dudt = du_dt(u, v, p, dx, j, xsteps)
    if j == 1
        dudt = 2 * (u(2) - u(1)) / dx^2 + p - u(1) - u(1) * v(1)^2;
    elseif j == xsteps
        dudt = 2 * (u(xsteps-1) - u(xsteps)) / dx^2 + p - u(xsteps) - u(xsteps) * v(xsteps)^2;
    else
        dudt = (u(j+1) + u(j-1) - 2 * u(j)) / dx^2 + p - u(j) - u(j) * v(j)^2;
    end
end

%% --- Function to compute dv/dt ---
function dvdt = dv_dt(u, v, delta, h, m, dx, j, xsteps)
    if j == 1
        dvdt = delta * 2 * (v(2) - v(1)) / dx^2 + u(1) * v(1)^2 * (1 - h * v(1)) - m * v(1);
    elseif j == xsteps
        dvdt = delta * 2 * (v(xsteps-1) - v(xsteps)) / dx^2 + u(xsteps) * v(xsteps)^2 * (1 - h * v(xsteps)) - m * v(xsteps);
    else
        dvdt = delta * (v(j+1) + v(j-1) - 2 * v(j)) / dx^2 + u(j) * v(j)^2 * (1 - h * v(j)) - m * v(j);
    end
end
