function [u, v] = simulate_Klausmeier(delta, h, m, pvec, t, tsteps, L, xsteps, noiseType, correlation_length, NoiseScale)
    %%Specific parameter scenarios for tipping and/or Turing can be found
    %%in the paper in table 1

    %% --- Create spatial and temporal grids ---
    dt = t / (tsteps - 1);
    dx = L / (xsteps - 1);
    time = 0:dt:t;
    place = 0:dx:L;
    [T, X] = meshgrid(time, place);

    %% --- Compute steady states u* and v* for intial solution---
    ustar = (2*h*m + pvec(1) + 2*h^2*pvec(1) - sqrt(-4*m^2 - 4*h*m*pvec(1) + pvec(1)^2)) / (2 + 2*h^2);
    vstar = (pvec(1) + sqrt(pvec(1)^2 - 4*m*(m + h*pvec(1)))) / (2*m + 2*h*pvec(1));

    %% --- Initialize solution matrices with steady states ---
    u = ones(tsteps, xsteps) * ustar;
    v = ones(tsteps, xsteps) * vstar;
    epsilon = NoiseScale * sqrt(dt);  % Noise scaling factor

    %spatial correlation filter for noise
    if correlation_length>0
    x= dx*(-xsteps:1:xsteps);
    filter=1/(correlation_length*sqrt(pi))*exp(-x.^2/correlation_length^2);
    end

    %% --- Time-stepping loop ---
    for i = 1:tsteps-1
        %% --- Generate noise based on type ---
        if strcmp(noiseType, 'Correlated')
            uncorrelated_noise_u = epsilon * randn(xsteps, 1);
            uncorrelated_noise_v = epsilon * randn(xsteps, 1);

            Noise = dx*conv(uncorrelated_noise_u,filter, 'same');
            Noise2 = dx*conv(uncorrelated_noise_v,filter, 'same');

        elseif strcmp(noiseType, 'White')
            Noise = epsilon * randn(xsteps, 1);
            Noise2 = epsilon * randn(xsteps, 1);

        else  % Uniform noise
            Noise = epsilon * (rand(xsteps, 1) - 0.5) * 2;
            Noise2 = epsilon * (rand(xsteps, 1) - 0.5) * 2;
        end

        %% --- Boundary conditions (free flow) ---
        % Left boundary
        u(i+1, 1) = u(i, 1) + dt * (2 * (u(i, 2) - u(i, 1)) / dx^2 + pvec(i) - u(i, 1) - u(i, 1) * v(i, 1)^2) + Noise(1);
        v(i+1, 1) = v(i, 1) + dt * (delta * 2 * (v(i, 2) - v(i, 1)) / dx^2 + u(i, 1) * v(i, 1)^2 * (1 - h * v(i, 1)) - m * v(i, 1)) + Noise2(1);

        % Right boundary
        u(i+1, xsteps) = u(i, xsteps) + dt * (2 * (u(i, xsteps-1) - u(i, xsteps)) / dx^2 + pvec(i) - u(i, xsteps) - u(i, xsteps) * v(i, xsteps)^2) + Noise(xsteps);
        v(i+1, xsteps) = v(i, xsteps) + dt * (delta * 2 * (v(i, xsteps-1) - v(i, xsteps)) / dx^2 + u(i, xsteps) * v(i, xsteps)^2 * (1 - h * v(i, xsteps)) - m * v(i, xsteps)) + Noise2(xsteps);

        %% --- Interior points (central difference) ---
        for j = 2:xsteps-1
            u(i+1, j) = u(i, j) + dt * ((u(i, j+1) + u(i, j-1) - 2 * u(i, j)) / dx^2 + pvec(i) - u(i, j) - u(i, j) * v(i, j)^2) + Noise(j);
            v(i+1, j) = v(i, j) + dt * (delta * (v(i, j+1) + v(i, j-1) - 2 * v(i, j)) / dx^2 + u(i, j) * v(i, j)^2 * (1 - h * v(i, j)) - m * v(i, j)) + Noise2(j);
        end
    end

    %% --- Optional: Plot results ---
    % xbegin = 1;
    % xend = xsteps;
    % timesteps = tsteps / 100;
    % figure();
    % mesh(T(xbegin:xend, 1:timesteps:end), X(xbegin:xend, 1:timesteps:end), real(transpose(u(1:timesteps:end, xbegin:xend))));
    % figure();
    % mesh(T(xbegin:xend, 1:timesteps:end), X(xbegin:xend, 1:timesteps:end), real(transpose(v(1:timesteps:end, xbegin:xend))));
end
