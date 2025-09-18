%% --- Plotting Stability Results from Simulations ---
% This script loads simulation results from subfolders and visualizes
% eigenvalue curves, statistical summaries, and parameter distributions.

close all
clear

% 📁 Define the folder containing your results
% You can use the same logic as in your simulation script:
% base_dir = pwd;
% subfolder_name = sprintf('Results_sim=%d_L=%d_xsteps=%d_Correlationlength=%.2f_NoiseAmplitude=%.2f', ...
%                          simuls, L, xsteps, correlation_length, NoiseScale);
% folder = fullfile(base_dir, subfolder_name);

baseDir = pwd;  % Use current working directory

for delta = [0.01, 0.5]
    % Set y-axis limits based on delta
    if delta == 0.01
        ymin = -2.5; ymax = 0.5;
    else
        ymin = -4; ymax = 0.5;
    end
    delta_str = sprintf('delta=%.2f', delta);

    % Get subfolders
    subs = dir(baseDir);
    subs = subs([subs.isdir] & ~ismember({subs.name}, {'.', '..'}));

    tracker2 = 0;
    legend_combined = cell(length(subs), 1);

    for foldernum = 1:length(subs)
        tracker2 = tracker2 + 1;

        %% --- Model parameters ---
        h = 0.1; m = 0.5; pbegin = 6; pend = 6;
        time_change = (pbegin ~= pend);

        % Compute steady states
        ustar = (2*h*m + pbegin + 2*h^2*pbegin - sqrt(-4*m^2 - 4*h*m*pbegin + pbegin^2)) / (2 + 2*h^2);
        vstar = (pbegin + sqrt(pbegin^2 - 4*m*(m + h*pbegin))) / (2*m + 2*h*pbegin);

        % Define theoretical stability matrix
        StabilityMatrix_real = @(k1) [-k1.^2 - 1 - vstar^2, -2*vstar*ustar;
                                      vstar^2*(1 - h*vstar), 2*vstar*ustar - 3*h*ustar*vstar^2 - m - delta*k1.^2];
        eigenvalue_function_real = @(k1) max(real(eig(StabilityMatrix_real(k1))));
        k_max_real = fminbnd(@(k1) -eigenvalue_function_real(k1), 0, 1000);
        max_eigenvalue_real = eigenvalue_function_real(k_max_real);

        %% --- Load simulation results ---
        folder = fullfile(baseDir, subs(foldernum).name);
        files = dir(fullfile(folder, ['*' delta_str '*.mat']));
        for f = 1:length(files)
            load(fullfile(folder, files(f).name));
        end

        simuls = length(k_max_data);
        ksize = 10;
        k1 = linspace(-ksize, ksize, 1000);
        eigenvalue_data = k_max_data;
        plot_data = zeros(simuls, length(k1));

        %% --- Reconstruct eigenvalue curves from fitted parameters ---
        for num = 1:simuls
            theta = theta_data(num, :);
            StabilityMatrix = @(k1) [theta(1) - theta(3) * k1.^2, theta(2);
                                     theta(5), theta(6) - k1.^2 * theta(8)];
            eigenvalue_function = @(k1) max(eig(StabilityMatrix(k1)));
            k_max = fminbnd(@(k1) -eigenvalue_function(k1), 0, 100);
            max_eigenvalue = eigenvalue_function(k_max);
            eigenvalue_data(num) = max_eigenvalue;

            for kpos = 1:length(k1)
                plot_data(num, kpos) = real(eigenvalue_function(k1(kpos)));
            end
        end

        %% --- Plot all eigenvalue curves ---
        figure;
        hold on
        for sim = 1:simuls
            plot(k1, plot_data(sim, :), 'LineWidth', 1);
        end
        ylim([ymin, ymax]);
        hold off
        savefig(fullfile(folder, ['all_plots_' delta_str '.fig']));

        %% --- Statistical summary plot ---
        meanValues = mean(plot_data, 1);
        percentile95 = prctile(plot_data, 95, 1);
        percentile5 = prctile(plot_data, 5, 1);

        figure;
        hold on;
        xlim([-ksize, ksize]);
        ylim([ymin, ymax]);

        scatter(k_max_data, eigenvalue_data, 30, 'filled', 'MarkerFaceAlpha', 0.3);

        %% --- Visualize spread in (k_max, lambda_max) ---
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

        %% --- Plot theoretical and statistical curves ---
        fplot(eigenvalue_function_real, [-ksize, ksize], 'k-', 'LineWidth', 1);
        if time_change
            fplot(eigenvalue_function_real_end, [-ksize, ksize], 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1);
        end
        plot(k1, meanValues, 'b', 'LineWidth', 2);
        plot(k1, percentile95, 'r--', 'LineWidth', 1.5);
        plot(k1, percentile5, 'g--', 'LineWidth', 1.5);
        plot(k_max_real, max_eigenvalue_real, 'p', 'MarkerSize', 20, 'MarkerFaceColor', 'yellow');

        %% --- Density overlays ---
        [f2, xi2] = ksdensity(eigenvalue_data, 'Bandwidth', 0.05);
        scaled_f2 = 0.25 * f2 / max(f2) * 2 * ksize - ksize;
        fill([scaled_f2, fliplr(scaled_f2)], [xi2, -ksize * ones(size(xi2))], ...
             'magenta', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        plot(scaled_f2, xi2, 'magenta-', 'LineWidth', 2);

        [f, xi] = ksdensity(k_max_data, 'Bandwidth', 0.05);
        scaled_f = 0.25 * f / max(f) * (ymax - ymin) + ymin;
        fill([xi, fliplr(xi)], [scaled_f, ymin * ones(size(scaled_f))], ...
             'magenta', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        plot(xi, scaled_f, 'magenta-', 'LineWidth', 2);

        %% --- Axes and annotation setup ---
        xlim([-ksize, ksize]);
        ylim([ymin, ymax]);

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
        savefig(fullfile(folder, ['statistical_summary_' delta_str '.fig']));
    end
end
