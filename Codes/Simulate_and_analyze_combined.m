%% --- Stability Analysis Across Multiple Simulations ---
% This script runs multiple simulations of the Klausmeier model using
% simulate_Klausmeier and analyze_data to evaluate fitted parameters and
% stability metrics under varying noise and diffusion conditions.

clear

% Loop over correlation lengths to define noise type
for correlation_length = [0, 0.1, 0.2]
    if correlation_length == 0
        noiseType = 'White';
    else 
        noiseType = 'Correlated';
    end

    % Loop over noise amplitudes
    for NoiseScale = [0.1, 1, 10]

        % Loop over diffusion coefficients
        for delta = [0.01, 0.5]

            simuls = 100;  % Number of simulations

            % Define wavenumber range for stability analysis
            if delta < 0.2
                ksize = 10;
            else 
                ksize = 2;
            end
            k1 = linspace(-ksize, ksize, 1000);

            % Initialize storage variables
            plot_data = zeros(simuls, length(k1));
            k_max_data = zeros(simuls, 1);
            theta_data = zeros(simuls, 8);
            data_u = cell(simuls, 1);
            data_v = cell(simuls, 1);

            %% --- Grid and time setup ---
            L = 40;              % Domain length
            xsteps = 400;        % Number of spatial steps
            t = 1;               % Total simulation time
            scale = 100;         % Scaling factor for time resolution
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
            Paul = true;
            tbegin = 1;
            tend = tsteps;
            tstepsize = 1;
            xbegin = 1;
            xend = xsteps;
            stepsize = 1;

            %% --- Run simulations in parallel ---
            parfor sim = 1:simuls
                rng(sim)  % Ensure reproducibility
                disp(['Simulation ', num2str(sim)]);

                % Simulate data
                [u, v] = simulate_Klausmeier(delta, h, m, pvec, t, tsteps, L, xsteps, noiseType, correlation_length, NoiseScale);

                % Analyze stability
                [theta, conf, StabilityMatrix, eigenvalue_function, k_max, max_eigenvalue] = ...
                    analyze_data(u, v, ustar, vstar, t, tsteps, L, xsteps, Paul, tbegin, tend, tstepsize, xbegin, xend, stepsize);

                % Store results
                data_u{sim} = u;
                data_v{sim} = v;
                k_max_data(sim) = abs(k_max);
                theta_data(sim, :) = theta;

                for kpos = 1:length(k1)
                    plot_data(sim, kpos) = real(eigenvalue_function(k1(kpos)));
                end
            end

            %% --- Save results to disk ---
            base_dir = pwd;  % Get current working directory
            subfolder_name = sprintf('Results_sim=%d_L=%d_xsteps=%d_Correlationlength=%.2f_NoiseAmplitude=%.2f', ...
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
            fplot(eigenvalue_function_real, [-ksize, ksize], 'k-', 'LineWidth', 1);
            if time_change
                fplot(eigenvalue_function_real_end, [-ksize, ksize], 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1);
            end
            plot(k1, meanValues, 'b', 'LineWidth', 2);
            plot(k1, percentile95, 'r--', 'LineWidth', 1.5);
            plot(k1, percentile5, 'g--', 'LineWidth', 1.5);

            % Vertical lines for theoretical k_max
            yl = ylim;
            plot([k_max_real, k_max_real], yl, 'k--', 'LineWidth', 1.5);
            plot(k_max_real, max_eigenvalue_real, 'p', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow');
            if time_change
                plot([k_max_real_end, k_max_real_end], yl, '--', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5);
                plot(k_max_real_end, max_eigenvalue_real_end, 'p', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow');
            end

            %% --- Density plot of k_max values ---
            [f, xi] = ksdensity(k_max_data, 'Bandwidth', 0.05);  % Estimate density of k_max values
            scaled_f = 0.5 * f / max(f) * (y_max - y_min) + (y_min - padding);  % Scale to match eigenvalue axis
            
            % Fill area under density curve
            fill([xi, fliplr(xi)], [scaled_f, (y_min - padding) * ones(size(scaled_f))], ...
                'magenta', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            plot(xi, scaled_f, 'magenta-', 'LineWidth', 2);
            
            % Set x-axis limits
            xlim([-ksize, ksize]);
            
            %% --- Add secondary Y-axis for density ---
            yyaxis right;
            set(gca, 'YColor', 'magenta');
            ylim([0, 2]);  % Adjust as needed
            yticks(0:0.2:1);
            
            %% --- Return to primary Y-axis for labels and legend ---
            yyaxis left;
            xlabel('Wavenumber k₁');
            ylabel('Stability Eigenvalue λ', 'Interpreter', 'latex');
            title('Statistical Summary of Stability Analysis');
            
            % Add legend depending on whether precipitation changes
            if time_change
                lgd = legend('Theoretical λ (start)', 'Theoretical λ (end)', ...
                             'Mean', '95th Percentile', '5th Percentile', ...
                             'Location', 'northwest');
            else
                lgd = legend('Theoretical λ', 'Mean', '95th Percentile', '5th Percentile', ...
                             'Theoretical k_{max}', 'Theoretical λ_{max}', ...
                             'Location', 'northwest');
            end
            set(lgd, 'FontSize', 8);
            
            %% --- Compute mean and uncertainty of selected θ parameters ---
            theta_mean = mean(theta_data(:, [1, 2, 5, 6, 8]), 1);
            theta_std = std(theta_data(:, [1, 2, 5, 6, 8]), 0, 1);
            
            % Format legend text with uncertainties
            labels = {'a', 'b', 'c', 'd', '\delta'};
            legend_lines = strings(1, numel(theta_mean));
            
            for i = 1:numel(theta_mean)
                val = theta_mean(i);
                unc = theta_std(i);
            
                if unc == 0
                    legend_lines(i) = sprintf('%s: %.4f(0)', labels{i}, val);
                    continue;
                end
            
                unc_sig = round(unc, -floor(log10(unc)) + 1);  % Round uncertainty to 1–2 sig digits
                digits = max(0, -floor(log10(unc_sig)) + 1);   % Decimal places needed
            
                val_rounded = round(val, digits);
                unc_digits = round(unc_sig * 10^digits);
            
                format_str = sprintf('%%s: %%.%df(%%d)', digits);
                legend_lines(i) = sprintf(format_str, labels{i}, val_rounded, unc_digits);
            end
            
            legend_text = strtrim(legend_lines);  % Clean up whitespace
            annotation('textbox', [0.14, 0.62, 0.3, 0.3], 'String', legend_text, ...
                'FitBoxToText', 'on', 'BackgroundColor', [0.678, 0.847, 0.902], ...
                'EdgeColor', 'blue', 'FontSize', 19, 'Margin', 0);
            
            hold off;
            savefig([submap_path 'statistical_summary_' filename_suffix '.fig']);
        end
    end
end
