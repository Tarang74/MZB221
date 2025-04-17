%% Phase Line (dx/dt vs. x(t))
clear; close all; clc;

dxdt = @(x) -(x + 1) .* x .* (x - 1).^2;
x = linspace(-2, 2, 2001); x(end) = [];

% Plot
figure
plot(x, dxdt(x))

% Export data
writetable(array2table([x.' dxdt(x).'], 'VariableNames', {'x', 'g'}), 'phase_line_analysis.csv', 'LineEnding', '\n')

%% Solution Curves (x(t) vs. t)
clear; close all; clc;

dxdt = @(t, x) -(x + 1) .* x .* (x - 1).^2;
t = linspace(0, 10, 2001); t(end) = [];

% Select initial values
x0s = [-1.01, -1, -2/3, -1/3, 0, 1/3, 2/3, 1, 1.99];

% Save data [t, x1, x2, ...]
export_data = nan([length(t), length(x0s) + 1]);
export_data(:, 1) = t.';

% Solve ODE numerically
for i = 1:length(x0s)
    [~, x] = ode45(dxdt, t, x0s(i));

    % Up to length(x) in cases of finite-time blow up
    export_data(1:length(x), i+1) = x;
end

% Plot
figure; hold on
for i = 1:length(x0s)
    plot(t, export_data(:, i+1), LineWidth=2)
end
ylim([-2, 2])

% Export data
writetable(array2table(export_data, 'VariableNames', ['t', arrayfun(@(k) sprintf("x%d", k), 1:length(x0s))]), 'solution_curve_analysis.csv', 'LineEnding', '\n')

%% Bifurcation Analysis
%% Phase Line Analysis (dx/dt vs. x(t))
clear; close all; clc;

dxdt = @(x, lambda) lambda .* x + x.^2 - x.^3;
x = linspace(-0.5, 1.5, 2001); x(end) = [];

lambdas = [-1/4 - 0.25,  ... % lambda < -1/4
           -1/4,         ... % lambda = -1/4
           -1/4 + 0.125, ... % -1/4 < lambda < 0
           0,            ... % lambda = 0
           0.25];            % lambda > 0

% Plot
figure; hold on
for lambda = lambdas
    plot(x, dxdt(x, lambda))
end
grid on
xlim([-1/2, 3/2])
ylim([-1/2, 1/2])

% Export data
export_data = cell2mat(arrayfun(@(lambda) dxdt(x, lambda).', lambdas, 'UniformOutput', false));
writetable(array2table([x.' export_data], 'VariableNames', ['x', arrayfun(@(k) sprintf("g%d", k), 1:length(lambdas))]), 'parametrised_phase_line.csv', 'LineEnding', '\n')

%% Solution Curves (x(t) vs. t)
clear; close all; clc;

dxdt = @(t, x, lambda) lambda .* x + x.^2 - x.^3;
t = linspace(0, 50, 2001); t(end) = [];

lambdas = [-1/4 - 0.25,  ... % lambda < -1/4
           -1/4,         ... % lambda = -1/4
           -1/4 + 0.125, ... % -1/4 < lambda < 0
           0,            ... % lambda = 0
           0.25];            % lambda > 0

% Equilibrium points:
xe1 = 0;
xe2 = @(lambda) 1/2 - 1/2 * sqrt(1 + 4*lambda);
xe3 = @(lambda) 1/2 + 1/2 * sqrt(1 + 4*lambda);

% Select initial values
x0cell = {[-0.5, xe1, 0.75, 1.5], ...                                      % lambda < -1/4
          [-0.5, xe1, ...                                                  % lambda = -1/4
           xe1 + 2/3*(xe3(lambdas(2)) - xe1), ...
           xe3(lambdas(2)), 1, 1.5], ...
          [-0.5, xe1, ...                                                  % -1/4 < lambda < 0
           xe1 + 2/3*(xe2(lambdas(3)) - xe1), ...
           xe2(lambdas(3)), ...
           xe2(lambdas(3)) + 1/3*(xe3(lambdas(3)) - xe2(lambdas(3))), ...
           xe2(lambdas(3)) + 2/3*(xe3(lambdas(3)) - xe2(lambdas(3))), ...
           xe3(lambdas(3)), 1.5], ...
          [-0.5, xe1, ...                                                  % lambda = 0
           xe1 + 1/3*(xe3(lambdas(4)) - xe1), ...
           xe1 + 2/3*(xe3(lambdas(4)) - xe1), ...
           xe3(lambdas(4)), 1.5], ...
          [-0.5, xe2(lambdas(5)), ...                                      % lambda > 0
           xe2(lambdas(5)) + (xe1 - xe2(lambdas(5)))/2, ...
           xe1, ...
           xe1 + 1/3*(xe3(lambdas(5)) - xe1), ...
           xe1 + 2/3*(xe3(lambdas(5)) - xe1), ...
           xe3(lambdas(5)), 1.5]};

% Plot
% for j = 1:length(lambdas)
%     lambda = lambdas(j);
%     x0s = x0cell{j};
% 
%     figure; hold on
%     for i = 1:length(x0s)
%         [~, x] = ode45(@(t, x) dxdt(t, x, lambda), t, x0s(i));
%         plot(t, x, LineWidth=2)
%     end
%     ylim([-0.5, 1.5])
% end

%% Bifurcation Diagram (x(t) vs. lambda)
clear; close all; clc;

dxdt = @(x, lambda) lambda .* x + x.^2 - x.^3;
dxdtdx = @(x, lambda) lambda + 2*x - 3*x.^2; % d/dx [dx/dt]

% Graph may be multivalued so consider the entire space of points
l = linspace(-1, 1, 2001); l(end) = []; x = linspace(-2, 2, 2001); x(end) = [];
[L, X] = meshgrid(l, x);

% Find regions of stability
% Stable: d/dx [dx/dt] < 0
stable = dxdt(X, L);
stable(dxdtdx(X, L) > 0) = nan;
% Unstable: d/dx [dx/dt] > 0
unstable = dxdt(X, L);
unstable(dxdtdx(X, L) < 0) = nan;

% Find contours
stable_segments = extract_contour_segments(contourc(l, x, stable, [0, 0]));
unstable_segments = extract_contour_segments(contourc(l, x, unstable, [0, 0]));

% Plot
figure
hold on
for i = 1:length(stable_segments) % Stable
    plot(stable_segments{i}(1, :), stable_segments{i}(2, :), 'k', LineWidth=1.5)
end
for i = 1:length(unstable_segments) % Unstable
    plot(unstable_segments{i}(1, :), unstable_segments{i}(2, :), 'k--', LineWidth=1.5)
end

% Combine segments
stable_all = [];
for k = 1:length(stable_segments)
    seg = stable_segments{k}';
    stable_all = [stable_all; seg; nan nan];
end
unstable_all = [];
for k = 1:length(unstable_segments)
    seg = unstable_segments{k}';
    unstable_all = [unstable_all; seg; nan nan];
end

figure
hold on
plot(stable_all(:, 1), stable_all(:, 2), 'k', LineWidth=1.5)
plot(unstable_all(:, 1), unstable_all(:, 2), 'k--', LineWidth=1.5)

% Export data
writetable(array2table(stable_all, 'VariableNames', {'x', 'y'}), 'bifurcation_stable.csv', 'LineEnding', '\n')
writetable(array2table(unstable_all, 'VariableNames', {'x', 'y'}), 'bifurcation_unstable.csv', 'LineEnding', '\n')

function [segments] = extract_contour_segments(C)
    segments = {};
    i = 1;
    while i < size(C, 2)
        n = C(2, i);
        segments{end+1} = C(:, i+1:i+n);
        i = i + n + 1;
    end
end
