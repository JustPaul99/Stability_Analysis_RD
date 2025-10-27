%% --- Stability Analysis Across Multiple Simulations ---
% This script runs multiple simulations of the Klausmeier model using
% simulate_Klausmeier and analyze_data to evaluate fitted parameters and
% stability metrics under varying noise and diffusion conditions.

clear
%Loops can be changed depending on the interest of your simulations.
% Loop over correlation lengths to define noise type
for correlation_length = [0, 0.1, 0.2]
    if correlation_length == 0
        noiseType = 'White';
    else 
        noiseType = 'Correlated';
    end

    % Loop over noise amplitudes
    for NoiseScale = [0.1,1,10]*sqrt(0.1003) %Fixed scaling of NoiseScales due to found error in normalisation in earlier code, scaling is presented to get same noisestrengths as in paper

        % Loop over diffusion coefficients
        for delta = [0.01,0.5]

            simuls = 100;  % Number of simulations, can make smaller for faster runs

            % Define wavenumber range for stability analysis
            if delta < 0.2
                ksize = 10;
            else 
                ksize = 2;
            end
            kpoints=1000;
            k1 = linspace(-ksize, ksize, kpoints);

            % Initialize storage variables
            plot_data = zeros(simuls, length(k1));
            k_max_data = zeros(simuls, 1);
            eigenvalue_data=zeros(simuls,1);
            theta_data = zeros(simuls, 8);
            data_u = cell(simuls, 1);
            data_v = cell(simuls, 1);

            %% --- Grid and time setup ---
            L = 40;              % Domain length
            xsteps = 400;        % Number of spatial steps
            t = 1;               % Total simulation time
            scale = 100;         % Scaling factor for time resolution, above 5 should be enough in most cases to work (improves runtime dramatically)
            tsteps = floor(scale * t * xsteps^2 / L^2);  % Time steps based on CFL-like scaling
            dt = t / (tsteps - 1);

            %% --- Model parameters ---
            h = 0.1;             % Feedback strength
            m = 0.5;             % Mortality rate
            pbegin = 6;          % Initial precipitation
            pend = 6;            % Final precipitation
            pvec = pbegin + (pend - pbegin) * (0:dt:t) / t;  % Time-varying precipitation
            time_change = (pbegin ~= pend);  % Flag for dynamic precipitation

            %% --- Compute steady states at t = 0 ---
            ustar = (2*h*m + pvec(1) + 2*h^2*pvec(1) - sqrt(-4*m^2 - 4*h*m*pvec(1) + pvec(1)^2)) / (2 + 2*h^2);
            vstar = (pvec(1) + sqrt(pvec(1)^2 - 4*m*(m + h*pvec(1)))) / (2*m + 2*h*pvec(1));

            %% --- Define theoretical stability matrix at t = 0 ---
            StabilityMatrix_real = @(k1) [-k1.^2 - 1 - vstar^2, -2*vstar*ustar;
                                          vstar^2*(1 - h*vstar), 2*vstar*ustar - 3*h*ustar*vstar^2 - m - delta*k1.^2];
            eigenvalue_function_real = @(k1) max(real(eig(StabilityMatrix_real(k1))));
            k_max_real = fminbnd(@(k1) -eigenvalue_function_real(k1), 0, 1000);
            max_eigenvalue_real = eigenvalue_function_real(k_max_real);

            %% --- If precipitation changes, compute steady states at t = end ---
            if time_change
                ustar_end = (2*h*m + pvec(end) + 2*h^2*pvec(end) - sqrt(-4*m^2 - 4*h*m*pvec(end) + pvec(end)^2)) / (2 + 2*h^2);
                vstar_end = (pvec(end) + sqrt(pvec(end)^2 - 4*m*(m + h*pvec(end)))) / (2*m + 2*h*pvec(end));

                StabilityMatrix_real_end = @(k1) [-k1.^2 - 1 - vstar_end^2, -2*vstar_end*ustar_end;
                                                  vstar_end^2*(1 - h*vstar_end), 2*vstar_end*ustar_end - 3*h*ustar_end*vstar_end^2 - m - delta*k1.^2];
                eigenvalue_function_real_end = @(k1) max(real(eig(StabilityMatrix_real_end(k1))));
                k_max_real_end = fminbnd(@(k1) -eigenvalue_function_real_end(k1), 0, 10);
                max_eigenvalue_real_end = eigenvalue_function_real_end(k_max_real);
            end

            %% --- Parameters for analyze_data ---
            tbegin = 1;
            tend = tsteps;
            tstepsize = 1;
            xbegin = 1;
            xend = xsteps;
            stepsize = 1;

            %% --- Run simulations in parallel ---
            parfor sim2 = 1:simuls
                rng(sim2) %for reproducability
                disp(['Simulation ', num2str(sim2)]);
                
                % Run the first function to simulate diffusion
                [u, v] = simulate_Klausmeier(delta, h, m, pvec, t, tsteps, L, xsteps, noiseType, correlation_length, NoiseScale);
                
                % Run the second function to analyze data
                [theta, conf, StabilityMatrix, eigenvalue_function, k_max, max_eigenvalue] = analyze_data(u, v, t, L, tbegin, tend, tstepsize, xbegin, xend, stepsize);
               
                data_u{sim2}=u;
                data_v{sim2}=v;
                k_max_data(sim2)=abs(k_max);
                eigenvalue_data(sim2)=max_eigenvalue;
                theta_data(sim2,:)=theta;
                for kpos = 1:kpoints
                    plot_data(sim2, kpos) = real(eigenvalue_function(k1(kpos)));
                end
            end

            %% --- Save results to disk ---
            %Feel free to change the subfolder names if you are changing
            %different parameter settings or you just want a specific
            %folder to save your files to. This one is set to distinguish
            %correlation length and noise amplitude now together with basic
            %setup settings.
            base_dir = pwd;  % Get current working directory
            subfolder_name = sprintf('Results_sim=%d_L=%d_xsteps=%d_Correlationlength=%.2f_NoiseAmplitude=%.2f/', ...
                         simuls, L, xsteps, correlation_length, NoiseScale);
            submap_path = fullfile(base_dir, subfolder_name);  % Construct full path

            if ~exist(submap_path, 'dir')
                mkdir(submap_path);
            end

            filename_suffix = sprintf('delta=%.2f', delta);
            save([submap_path 'theta_data_' filename_suffix '.mat'], 'theta_data');
            save([submap_path 'plot_data_' filename_suffix '.mat'], 'plot_data');
            save([submap_path 'k_max_data_' filename_suffix '.mat'], 'k_max_data');
            save([submap_path 'u_data_' filename_suffix '.mat'], 'data_u');
            save([submap_path 'v_data_' filename_suffix '.mat'], 'data_v');

            %% --- Plot all eigenvalue curves ---
            figure;
            hold on
            for sim = 1:simuls
                plot(k1, plot_data(sim, :), 'LineWidth', 1);
            end
            hold off
            savefig([submap_path 'all_plots_' filename_suffix '.fig']);

        %% --- Statistical summary plot ---
        meanValues = mean(plot_data, 1);
        percentile95 = prctile(plot_data, 95, 1);
        percentile5 = prctile(plot_data, 5, 1);

        figure;
        hold on;

        %% --- Plot theoretical and statistical curves ---
        fplot(eigenvalue_function_real, [-ksize, ksize], 'k-', 'LineWidth', 1);
        if time_change
            fplot(eigenvalue_function_real_end, [-ksize, ksize], 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1);
        end
        plot(k1, meanValues, 'b', 'LineWidth', 2);
        plot(k1, percentile95, 'r--', 'LineWidth', 1.5);
        plot(k1, percentile5, 'g--', 'LineWidth', 1.5);
        plot(k_max_real, max_eigenvalue_real, 'p', 'MarkerSize', 20, 'MarkerFaceColor', 'yellow');

        %% Set consistent Y-limits after plotting
        %You can set your own limits or let the limits be defined by your
        %dispersion relation as is commented out below
        xlim([-ksize, ksize]);
        y_min=-4;
        y_max=0.5;
        %ymin = min([percentile5, percentile95, meanValues]);
        %ymax = max([percentile5, percentile95, meanValues]);
        %padding = 0.05 * (ymax - ymin);  % Small padding
        %y_min=ymin-padding;
        %y_max=ymin-padding;
        ylim([y_min, y_max]);


        %% --- Visualize spread in (k_max, lambda_max) --f-
        scatter(k_max_data, eigenvalue_data, 30, 'filled', 'MarkerFaceAlpha', 0.3);

        data = [k_max_data, eigenvalue_data];
        plotvariation = 3;
        mu = mean(data);
        sigma = std(data);

        if sigma(1) < 1e-5 || sigma(2) < 1e-5 || plotvariation == 1 %rectangular
            rectangle('Position', [mu(1)-sigma(1), mu(2)-sigma(2), 2*sigma(1), 2*sigma(2)], ...
                      'EdgeColor', [1, 0.5, 0], 'LineWidth', 2);
        elseif plotvariation == 2 %polynomial least area fit
            C = cov(data);
            mahal_dist = sqrt(sum(((data - mu) / chol(C)).^2, 2));
            inliers = data(mahal_dist <= 1, :);
            if size(inliers, 1) >= 3
                k = convhull(inliers(:,1), inliers(:,2));
                plot(inliers(k,1), inliers(k,2), 'y-', 'LineWidth', 2);
            end
        else %Elliptic fit via covariance
            C = cov(data);
            [V, D] = eig(C);
            theta = linspace(0, 2*pi, 100);
            circle = [cos(theta); sin(theta)];
            ellipsoid = V * sqrt(D) * circle + mu';
            plot(ellipsoid(1,:), ellipsoid(2,:), '-', 'Color', [1, 0.5, 0], 'LineWidth', 2);
        end


         %% --- Density overlays ---
        [f2, xi2] = ksdensity(eigenvalue_data, 'Bandwidth', 0.05);
        scaled_f2 = 0.25 * f2 / max(f2) * 2 * ksize - ksize;
        fill([scaled_f2, fliplr(scaled_f2)], [xi2, -ksize * ones(size(xi2))], ...
             'magenta', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        plot(scaled_f2, xi2, 'magenta-', 'LineWidth', 2);

        [f, xi] = ksdensity(k_max_data, 'Bandwidth', 0.05);
        scaled_f = 0.25 * f / max(f) * (y_max - y_min) + y_min;
        fill([xi, fliplr(xi)], [scaled_f, (y_min) * ones(size(scaled_f))], ...
             'magenta', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        plot(xi, scaled_f, 'magenta-', 'LineWidth', 2);

        %% --- Axes and annotation setup ---

        yyaxis right;
        set(gca, 'YColor', 'magenta');
        ylim([0, 4]);
        yticks(0:0.5:1);
        
        % Match x-axis limits and ticks on top overlay
        axTop = axes('Position', get(gca, 'Position'), ...
                     'Color', 'none', 'XAxisLocation', 'top', ...
                     'YAxisLocation', 'left', 'YTick', [], ...
                     'XColor', 'magenta', 'YColor', 'none');
        xlim(axTop, [0, 4]);
        xticks(axTop, 0:0.5:1);
            
            %% --- Return to primary Y-axis for labels and legend ---         
            
%% --- Annotate fitted parameters (theta and stability metrics) ---
        theta_mean = mean(theta_data(:, [1, 2, 5, 6, 8]), 1);
        theta_std = std(theta_data(:, [1, 2, 5, 6, 8]), 0, 1);

        labeloptions = 2;  % Choose label format: 1 = full theta with stability parameters, 2 = stability parameters, 3 = raw theta only

        if labeloptions == 1
            average = [theta_mean, mu];
            deviation = [theta_std, std(data)];
            labels = {'a', 'b', 'c', 'd', '\delta', 'k_*', '\lambda_*'};

        elseif labeloptions == 2
            average = mu - [k_max_real, max_eigenvalue_real];
            deviation = std(data);
            labels = {'\Delta{k_*}', '\Delta\lambda_*'};

        else
            average = theta_mean;
            deviation = theta_std;
            labels = {'a', 'b', 'c', 'd', '\delta'};
        end

        legend_lines = strings(1, numel(average));
        for i = 1:numel(average)
            val = average(i);
            unc = deviation(i);

            if unc < 1e-4
                legend_lines(i) = sprintf('%s: %.4f(0)', labels{i}, val);
                continue;
            end

            unc_sig = round(unc, -floor(log10(unc)) + 1);
            digits = max(0, -floor(log10(unc_sig)) + 1);
            val_rounded = round(val, digits);
            unc_digits = round(unc_sig * 10^digits);

            format_str = sprintf('%%s: %%.%df(%%d)', digits);
            legend_lines(i) = sprintf(format_str, labels{i}, val_rounded, unc_digits);
        end

        legend_text = cellstr(strtrim(legend_lines));

        % Dynamically size and position annotation box
        ax = axes('Units', 'normalized', 'Position', [0 0 1 1], 'Visible', 'off');
        temp_text = text(0.5, 0.5, legend_text, 'Units', 'normalized', 'FontSize', 19, 'Visible', 'off');
        extent = get(temp_text, 'Extent');
        text_width = extent(3);
        delete(temp_text); delete(ax);

        annotation('textbox', [0.86 - text_width, 0.82, 0.1, 0.1], ...
            'String', legend_text, 'FitBoxToText', 'on', ...
            'BackgroundColor', [0.678, 0.847, 0.902], ...
            'EdgeColor', 'blue', 'FontSize', 19, 'Margin', 0);

        hold off;
        savefig([submap_path 'statistical_summary' filename_suffix '.fig']);
        end
    end
end
